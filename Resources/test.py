"""
Minimal reproduction of Antescofo's StateChain score parsing.

From the paper (Section 2.3 + 3.1) and debug output, the structure is:

  Score Event (macro)
       |
  State (semi-Markov macrostate)
       |
       +-- Microstate 1  (Markov, same pitch)
       +-- Microstate 2  (Markov, same pitch)
       +-- Microstate 3  (Markov, same pitch)
       ...

The StateChain sees ONE logical State per score event.
Each State owns r Markov microstates internally.
During inference, the forward variable iterates over occupancy u,
accumulating microstate emissions as the implicit duration prior
(negative binomial, eq. 2 in paper), while the semi-Markov
survival distribution (eq. 19) provides the adaptive tempo-driven
duration model on top.

Score:
    BPM 80
    NOTE C4 2   -> type=0, pitch=6000 (C4 ~= 261.6 Hz, stored as int mHz or cents)
    NOTE D4 2   -> type=0, pitch=6200
    NOTE E4 2   -> type=0, pitch=6400

Debug shows:
    chain size = 5  (initial silence + 3 notes + terminal empty state)
    type=1  -> silence/rest,  1 sub-state,  pitch=0
    type=0  -> normal note,   r sub-states, all same pitch
    type=2  -> chord,         1 combined + N individual sub-states
"""

import math
from dataclasses import dataclass, field
from typing import List, Optional


# ---------------------------------------------------------------------------
# Microstate: one node in the internal Markov chain of a macrostate
# ---------------------------------------------------------------------------

@dataclass
class Microstate:
    """
    A single Markov microstate inside a semi-Markov macrostate.
    All microstates within the same macrostate share the same pitch(es),
    so they have identical observation probabilities.
    From Figure 2 in the paper: r microstates with transition prob q
    (self-loop / next) and p (exit). With p=0 this gives the negative
    binomial occupancy distribution.
    """
    pitches: List[int]          # pitch list in mHz (0 = silence)
    q: float = 0.9             # next-microstate / self-loop probability
    p: float = 0.0             # exit probability (0 = no early exit)

    def __post_init__(self):
        print(f"State::State: add new event ( {' '.join(str(p) for p in self.pitches)} )")


# ---------------------------------------------------------------------------
# State: one semi-Markov macrostate, owning r Markov microstates
# ---------------------------------------------------------------------------

class State:
    """
    Semi-Markov macrostate corresponding to one score event.

    event_type:
        0 -> normal timed note    (semi-Markov, r microstates)
        1 -> silence / rest       (semi-Markov, 1 microstate, pitch=0)
        2 -> chord / polyphonic   (semi-Markov, 1 combined + N individual)

    The number of microstates r is derived from the note duration and
    the BPM (minimum frame occupancy). Here we hardcode r=3 as seen
    in the debug output for quarter/half notes at 80 BPM.

    The survival distribution d_j(u) = exp(-(t - t_{n-1}) / (psi_{n-1} * l_j))
    (equation 19) is computed externally by the tempo agent and passed in
    during inference. The microstates provide the implicit duration floor.
    """

    # r=3 microstates is what Antescofo uses for normal notes at typical BPM.
    # For the last note in the score an extra terminal microstate is added,
    # which is why E4 shows 4 sub-states in the debug output.
    DEFAULT_R = 3

    def __init__(self,
                 event_type: int,
                 pitches: List[int],
                 duration_beats: float,
                 index: int,
                 is_last: bool = False):
        self.event_type = event_type
        self.pitches = pitches          # all pitches of this event
        self.duration_beats = duration_beats
        self.index = index
        self.is_last = is_last

        # forward variable (set during inference)
        self.forward: float = 0.0
        self.forward_last: float = 0.0

        # build internal microstates
        self.microstates: List[Microstate] = self._build_microstates()

    def _build_microstates(self) -> List[Microstate]:
        microstates = []

        if self.event_type == 1:
            # Silence: single microstate, pitch=0
            microstates.append(Microstate(pitches=[0]))

        elif self.event_type == 0:
            # Normal timed note: r identical microstates, same pitch
            # The last note in the score gets r+1 to allow transition
            # to the terminal state (this is what produces the 4th sub-state
            # for E4 in the debug output).
            r = self.DEFAULT_R + (1 if self.is_last else 0)
            for _ in range(r):
                microstates.append(Microstate(pitches=self.pitches))

        elif self.event_type == 2:
            # Chord: one combined state + one per individual pitch
            # Combined state uses max observation over all pitches (TRILL logic)
            microstates.append(Microstate(pitches=self.pitches))       # combined
            for p in self.pitches:
                microstates.append(Microstate(pitches=[p]))            # individual

        return microstates

    def observation(self, audio_pitch_hz: float) -> float:
        """
        Simplified KL-divergence based observation probability (eq. 6).
        Here we just use a Gaussian over pitch distance as a stand-in.
        In Antescofo this is computed from FFT frames vs pitch templates.
        """
        if not self.pitches or self.pitches == [0]:
            # Silence state: high obs when energy is low
            return 0.1

        # Use the best-matching pitch (max observation, as in TRILL class)
        best = max(
            self._pitch_obs(p, audio_pitch_hz)
            for p in self.pitches if p > 0
        )
        return best

    @staticmethod
    def _pitch_obs(state_pitch_mhz: int, audio_hz: float) -> float:
        """
        Gaussian on log-frequency distance. state_pitch in milli-Hz.
        """
        state_hz = state_pitch_mhz / 1000.0
        if state_hz <= 0 or audio_hz <= 0:
            return 0.0
        # distance in semitones
        dist_semitones = abs(12.0 * math.log2(audio_hz / state_hz))
        sigma = 0.5  # half-tone variance, as mentioned in paper Section 6
        return math.exp(-0.5 * (dist_semitones / sigma) ** 2)

    def survival(self, u: int, tempo_psi: float) -> float:
        """
        Survival distribution d_j(u) from equation 19:
            d_j(t - t_{n-1}) = exp(-(t - t_{n-1}) / (psi_{n-1} * l_j))
        u      : discrete time steps elapsed in this state
        tempo_psi : current tempo in seconds/beat (from tempo agent)
        """
        expected_duration = tempo_psi * self.duration_beats   # seconds
        # convert u (frame count) to seconds assuming ~23ms frames
        frame_seconds = 0.023
        t_elapsed = u * frame_seconds
        return math.exp(-t_elapsed / expected_duration) if expected_duration > 0 else 0.0

    def __repr__(self):
        return (f"State(type={self.event_type}, pitches={self.pitches}, "
                f"dur={self.duration_beats}b, microstates={len(self.microstates)})")


# ---------------------------------------------------------------------------
# Terminal state: empty pitch list, as seen in debug: state ( )
# ---------------------------------------------------------------------------

class TerminalState(State):
    def __init__(self, index: int):
        # bypass normal __init__ to avoid microstate creation
        self.event_type = -1
        self.pitches = []
        self.duration_beats = 0.0
        self.index = index
        self.is_last = True
        self.forward = 0.0
        self.forward_last = 0.0
        self.microstates = []
        print("State::State: add new event (  )")   # empty, as in debug

    def observation(self, audio_pitch_hz: float) -> float:
        return 0.0  # terminal absorbs nothing

    def __repr__(self):
        return "State(terminal, pitches=[])"


# ---------------------------------------------------------------------------
# StateChain: builds and holds the full left-right state sequence
# ---------------------------------------------------------------------------

class StateChain:
    """
    Builds the Hidden Hybrid Markov/semi-Markov chain from score events.

    The chain is strictly left-right (no backward transitions), which
    allows the real-time forward-only inference of equation (3).

    Structure for BPM 80 / NOTE C4 2 / NOTE D4 2 / NOTE E4 2:

        State 0: silence   type=1  pitches=[0]       microstates=1
        State 1: C4        type=0  pitches=[6000]    microstates=3
        State 2: D4        type=0  pitches=[6200]    microstates=3
        State 3: E4        type=0  pitches=[6400]    microstates=4  (last)
        State 4: terminal  type=-1 pitches=[]        microstates=0

    This matches exactly:
        chain size = 5
        score with 3 events (silence is implicit)
    """

    def __init__(self, bpm: float = 80.0):
        self.bpm = bpm
        self.tempo_psi = 60.0 / bpm   # seconds per beat
        self.states: List[State] = []
        self._add_initial_silence()

    def _add_initial_silence(self):
        """
        Antescofo always prepends an implicit silence state (type=1).
        This is the state ( 0 ) seen at the top of every inference trace.
        """
        s = State(event_type=1, pitches=[0], duration_beats=0.0, index=0)
        self.states.append(s)

    def newStateFomEvent(self,
                         event_type: int,
                         pitches: List[int],
                         duration_beats: float,
                         is_last: bool = False):
        """
        Mirrors StateChain::newStateFomEvent() in Antescofo.
        Called once per score event during parsing.
        """
        print(f"\nStateChain::newStateFomEvent(): event->type() = {event_type}")
        index = len(self.states)
        s = State(event_type=event_type,
                  pitches=pitches,
                  duration_beats=duration_beats,
                  index=index,
                  is_last=is_last)
        self.states.append(s)

    def finalize(self):
        """
        After all events are added, append the terminal empty state.
        Also triggers the decode window reset seen in the debug output.
        """
        terminal = TerminalState(index=len(self.states))
        self.states.append(terminal)

    def print_chain(self):
        print("\n--- StateChain contents ---")
        for s in self.states:
            print(f"  {s}")

    def infer(self, audio_pitch_hz: float):
        """
        Simplified forward inference (equation 3/4).
        This is a skeleton — the real version integrates the survival
        distribution and tempo agent. Here we just show the structure.
        """
        print(f"\nStateChain::Infere(): BEGIN: size of N_t vector = {len(self.states)}")

        # Propagate forward variable left-to-right
        for i, state in enumerate(self.states):
            obs = state.observation(audio_pitch_hz)

            if i == 0:
                # Initial state keeps its prior
                new_forward = state.forward_last * obs
            else:
                prev = self.states[i - 1]
                # Semi-Markov forward: max over occupancy u (simplified to 1 step here)
                # Full version integrates survival d_j(u) as in eq. 3
                state_in = prev.forward_last
                surv = state.survival(u=1, tempo_psi=self.tempo_psi)
                new_forward = obs * state_in * surv

            state.forward = new_forward
            pitch_str = ' '.join(str(p) for p in state.pitches)
            print(f"StateChain::Infere: state ( {pitch_str} ) - "
                  f"ForwardLast = {state.forward_last:.6g}, "
                  f"obs = {obs:.6f}, "
                  f"Forward = {new_forward:.6g}")

        # Decode: argmax forward
        best = max(self.states, key=lambda s: s.forward)
        pitch_str = ' '.join(str(p) for p in best.pitches)
        print(f"\nStateChain::Decode(): BEST forward found: "
              f"- sub-state pitch list ( {pitch_str} ) "
              f"- state pitch list: ( {pitch_str} )")
        print(f"AnteAlign::Decode(): most likely event = index() {best.index}, "
              f"pitch list ( {pitch_str} ), type = {best.event_type}")

        # Update forward_last for next frame
        for state in self.states:
            state.forward_last = state.forward


# ---------------------------------------------------------------------------
# Parse the score and build the chain
# ---------------------------------------------------------------------------

def parse_score():
    """
    BPM 80
    NOTE C4 2    -> type=0, pitch=6000 mHz (C4 = 261.6 Hz -> ~6000 in Antescofo units)
    NOTE D4 2    -> type=0, pitch=6200 mHz
    NOTE E4 2    -> type=0, pitch=6400 mHz  (last note)
    """
    score = [
        # (event_type, pitches,  duration_beats, is_last)
        (0, [6000], 2.0, False),   # C4
        (0, [6200], 2.0, False),   # D4
        (0, [6400], 2.0, True),    # E4  <- last, gets extra microstate
    ]

    chain = StateChain(bpm=80.0)

    for i, (etype, pitches, dur, is_last) in enumerate(score):
        chain.newStateFomEvent(etype, pitches, dur, is_last)

    chain.finalize()
    chain.print_chain()
    return chain


if __name__ == "__main__":
    print("=" * 60)
    print("SCORE PARSING")
    print("=" * 60)
    chain = parse_score()

    print("\n" + "=" * 60)
    print("INFERENCE (playing C4 = 261.6 Hz)")
    print("=" * 60)

    # Initialize: silence state gets prior belief 1.0
    chain.states[0].forward_last = 1.0

    # Simulate a few frames of C4
    for frame in range(3):
        print(f"\n--- Frame {frame} ---")
        chain.infer(audio_pitch_hz=261.6)   # C4
