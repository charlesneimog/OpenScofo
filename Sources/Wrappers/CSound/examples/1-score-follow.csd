<CsoundSynthesizer>

<CsInstruments>
sr = 48000
ksmps = 64
nchnls = 1
0dbfs = 1

instr 1
    ; audio input from file or real time audio input
    a1 diskin2 "/home/neimog/Documents/Git/OpenScofo/Resources/Pieces/miniatura1.wav", 1, 0, 0

    ; kEvent is the event index of the score, first note, second note, third note, etc (rest do not count as event).
    ; kBPM is the current kBPM that the player is play a given score.
    ; kTrig is 1 when new event index is detected, 0 otherwise.
    kEvent, kBPM, kTrig OpenScofoScore a1, "/home/neimog/Documents/Git/OpenScofo/Tests/miniaturas/Extras/miniatura1-csound.scofo", 2048, 512
        
    ; first input (a1) is the audio signal (recording or realtime audio from the MIC)
    ; second is the Score Text (check https://charlesneimog.github.io/OpenScofo/score/intro/ for Documentation)
    ; FFT Size (recommended is 2048)
    ; Hop Size (recommended is 512)


		; imprime apenas quando kTrig for 1
		printf "Event: %03d | BPM: %.2f\n", kTrig, kEvent, kBPM
		
   out a1
endin

instr 2
    iFreq = p4
    if iFreq <= 0 then
        iFreq = 440
    endif
    iAmp = p5
    if iAmp <= 0 then
        iAmp = 0.20
    endif
    aEnv linseg 0, 0.01, iAmp, p3 - 0.02, iAmp, 0.01, 0
    aTone poscil aEnv, iFreq
    //prints "OpenScofo scheduled instr 2: freq=%f amp=%f\\n", iFreq, iAmp
    out aTone
endin

instr namedPing
    iFreq = p4
    if iFreq <= 0 then
        iFreq = 880
    endif
    iAmp = p5
    if iAmp <= 0 then
        iAmp = 0.12
    endif
    aEnv linseg 0, 0.01, iAmp, p3 - 0.02, iAmp, 0.01, 0
    aTone poscil aEnv, iFreq
    //prints "OpenScofo scheduled namedPing: freq=%f amp=%f\\n", iFreq, iAmp
    out aTone
endin

</CsInstruments>

<CsScore>
i1 0 30
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
