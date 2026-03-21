OpenScofo : UGen {
    *ar { arg sampleRate = 48000.0, fftSize = 2048.0, hopSize = 512.0, in = 0.0;
        ^this.multiNew('audio', in, sampleRate, fftSize, hopSize)
    }

    *kr { arg sampleRate = 48000.0, fftSize = 2048.0, hopSize = 512.0, in = 0.0;
        ^this.multiNew('control', in, sampleRate, fftSize, hopSize)
    }

    *findUnitIndex { arg synth;
        var desc, children, openScofoUGen;

        if(synth.isNil) { ^0 };

        desc = SynthDescLib.global.at(synth.defName.asSymbol);
        if(desc.isNil or: { desc.def.isNil }) { ^0 };

        children = desc.def.children;
        if(children.isNil) { ^0 };

        openScofoUGen = children.detect({ |ugen|
            ugen.class == OpenScofo
        });

        if(openScofoUGen.isNil) { ^0 };
        ^openScofoUGen.synthIndex ? 0;
    }

    *cmd { arg synth, command ...args;
        var unitIndex;
        if(synth.isNil) {
            Error("OpenScofo requires a Synth instance").throw;
        };
        unitIndex = this.findUnitIndex(synth);
        synth.server.sendMsg("/u_cmd", synth.nodeID, unitIndex, command.asString, *args);
    }

    *parseScore { arg synth, path;
        ^this.cmd(synth, "parseScore", path.asString);
    }

    *getCurrentEvent { arg synth;
        ^this.cmd(synth, "getCurrentEvent");
    }

    *setFollowScore { arg synth, enabled = true;
        var value = if(enabled, 1, 0);
        ^this.cmd(synth, "setFollowScore", value);
    }

    *setEventNotifications { arg synth, enabled = true;
        var value = if(enabled, 1, 0);
        ^this.cmd(synth, "setEventNotifications", value);
    }

    *getDescriptor { arg synth, descriptorId;
        ^this.cmd(synth, "getDescriptor", descriptorId.asString);
    }

    *loadOnnxModel { arg synth, modelPath, descriptorIds;
        var csv = descriptorIds.asArray.collect(_.asString).join(",");
        ^this.cmd(synth, "loadOnnxModel", modelPath.asString, csv);
    }
}
