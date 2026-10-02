---
icon: custom/csound
tags:
  - Host Integration
  - Csound
---

# Csound

In Csound, `sendto` schedules an instrument.

```openscofo
NOTE C4 1
    sendto 2 [0 1 0.7]
```

```csound
<CsoundSynthesizer>

<CsInstruments>
sr = 48000
ksmps = 64
nchnls = 1
0dbfs = 1

instr 1
    a1 diskin2 "./miniatura1.wav", 1, 0, 0
    kEvent, kBPM, kTrig OpenScofoScore a1, "./miniatura1-csound.scofo", 2048, 512
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
    out aTone
endin

</CsInstruments>

<CsScore>
i1 0 30
</CsScore>
</CsoundSynthesizer>
```
