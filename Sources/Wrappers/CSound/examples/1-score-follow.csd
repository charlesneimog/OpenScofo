<CsoundSynthesizer>
<CsOptions>
-d -odac --opcode-lib=/home/neimog/Documents/Git/OpenScofo/build/Sources/Wrappers/CSound/OScofoCSound.so
</CsOptions>

<CsInstruments>

sr = 48000
ksmps = 512
nchnls = 1
0dbfs = 1


; This intrument load the score and do something
instr 1
    ; audio input from file or real time audio input
    a1 diskin2 "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.wav", 1, 0, 0
    kEvent, kBPM, kTrig, kDesc OpenScofo a1, "score", "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.txt"
    ; print when event updates
   if (kTrig == 1) then
   		printf "On event index: %d | BPM: %.2f\n", changed(kEvent), kEvent, kBPM
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
