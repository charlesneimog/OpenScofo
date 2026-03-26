<CsoundSynthesizer>
<CsOptions>
--opcode-lib=/home/neimog/Documents/Git/OpenScofo/build/Sources/Wrappers/CSound/OScofoCSound.so
</CsOptions>

<CsInstruments>
sr = 48000
ksmps = 512
nchnls = 1
0dbfs = 1

instr 1
    ; audio input from file or real time audio input
    a1 diskin2 "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.wav", 1, 0, 0

    ; kEvent is the event index of the score, first note, second note, third note, etc (rest do not count as event).
    ; kBPM is the current kBPM that the player is play a given score.
    ; kTrig is 1 when new event index is detected, 0 otherwise.
    kEvent, kBPM, kTrig, OpenScofoScore a1, "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.txt"

    ; first input parameters is the audio channel
    ; third parameter is the score path for "score" mode or descriptor name for "description" mode.

    if (kTrig == 1) then
        ; here you can do something
        printk2 kEvent
    endif

    out a1
endin

</CsInstruments>

<CsScore>
i1 0 60
</CsScore>
</CsoundSynthesizer>
<bsbPanel>
 <label>Widgets</label>
 <objectName/>
 <x>100</x>
 <y>100</y>
 <width>320</width>
 <height>240</height>
 <visible>true</visible>
 <uuid/>
 <bgcolor mode="background">
  <r>240</r>
  <g>240</g>
  <b>240</b>
 </bgcolor>
</bsbPanel>
<bsbPresets>
</bsbPresets>
