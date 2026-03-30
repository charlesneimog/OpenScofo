<CsoundSynthesizer>
<CsOptions>
-odac --opcode-lib=/home/neimog/Documents/Git/OpenScofo/build/Sources/Wrappers/CSound/OpenScofo.so -b64
</CsOptions> 

<CsInstruments>

sr = 48000
ksmps = 512
nchnls = 1
0dbfs = 1


; This intrument load the score and do something
instr 1
    ; audio input from file or real time audio input
    a1 soundin "./assets/bwv-1013.wav", 1, 0, 0
    kEvent, kBPM, kTrig OpenScofoScore a1, "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.txt"
    ; print when event updates
   if (kTrig == 1) then
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
