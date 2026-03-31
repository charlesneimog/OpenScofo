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
    a1 diskin2 "./assets/1.wav", 1, 0, 1

k1 downsamp a1
printk2 k1

endin

</CsInstruments>

<CsScore>
i1 0 60
</CsScore>
</CsoundSynthesizer>
