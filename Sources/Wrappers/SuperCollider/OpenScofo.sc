OpenScofoUGen : UGen {
    *ar { arg in = 0.0, sampleRate = 48000.0;
        ^this.multiNew('audio', in, sampleRate)
    }

    *kr { arg in = 0.0, sampleRate = 48000.0;
        ^this.multiNew('control', in, sampleRate)
    }
}

OpenScofo {
    var <server, <synth, <scorePath, <sampleRate, <target, <addAction;
    var <inBus, <outBus, <unitIndex, <eventOSCDef;
    var eventAction, synthDefName, eventOSCDefName, actionOSCDefs;

    *new { arg scorePath, inBus = 0, outBus = 0, sampleRate, server, target, addAction = \addToTail,
        eventNotifications = true, eventAction;
        ^super.new.init(scorePath, inBus, outBus, sampleRate, server, target, addAction, eventNotifications, eventAction)
    }

    findUnitIndex {
        var desc, children, openScofoUGen;
        if(synth.isNil) { ^0 };

        desc = SynthDescLib.global.at(synth.defName.asSymbol);
        if(desc.isNil or: { desc.def.isNil }) { ^0 };

        children = desc.def.children;
        if(children.isNil) { ^0 };

        openScofoUGen = children.detect({ |ugen| ugen.class == OpenScofoUGen });
        if(openScofoUGen.isNil) { ^0 };
        ^openScofoUGen.synthIndex ? 0;
    }

    init { arg scorePathArg, inBusArg = 0, outBusArg = 0, sampleRateArg, serverArg, targetArg,
        addActionArg = \addToTail, eventNotifications = true, eventActionArg;
        var sr;

        server = serverArg ? Server.default;
        scorePath = scorePathArg;
        inBus = inBusArg;
        outBus = outBusArg;
        sampleRate = sampleRateArg ? server.sampleRate;
        target = targetArg ? server.defaultGroup;
        addAction = addActionArg;
        eventAction = eventActionArg;
        actionOSCDefs = IdentityDictionary.new;

        sr = sampleRate;
        synthDefName = ("OpenScofo_" ++ this.identityHash).asSymbol;
        eventOSCDefName = ("OpenScofo_event_" ++ this.identityHash).asSymbol;

        SynthDef(synthDefName, { arg inBus = 0, outBus = 0, sampleRate = sr;
            var in, tracker;
            in = In.ar(inBus, 1);
            tracker = OpenScofoUGen.ar(in, sampleRate);
            Out.ar(outBus, Silent.ar(2));
        }).add;

        server.sync;
        synth = Synth(synthDefName, [
            \inBus, inBus,
            \outBus, outBus,
            \sampleRate, sampleRate
        ], target: target, addAction: addAction);
        server.sync;

        unitIndex = this.findUnitIndex;

        if(scorePath.notNil) {
            this.loadScore(scorePath);
        };
        this.setEventNotifications(eventNotifications);
        this.registerEventHandler(eventAction);
        ^this
    }

    cmd { arg command, args = #[];
        if(synth.isNil) {
            Error("OpenScofo instance has no Synth").throw;
        };
        synth.server.sendMsg("/u_cmd", synth.nodeID, unitIndex, command.asString, *args.asArray);
    }

    loadScore { arg path;
        scorePath = path;
        ^this.cmd("loadScore", [path.asString]);
    }

    parseScore { arg path;
        ^this.loadScore(path);
    }

    getCurrentEvent {
        ^this.cmd("getCurrentEvent");
    }

    setFollowScore { arg enabled = true;
        ^this.cmd("setFollowScore", [if(enabled, 1, 0)]);
    }

    setEventNotifications { arg enabled = true;
        ^this.cmd("setEventNotifications", [if(enabled, 1, 0)]);
    }

    getDescriptor { arg descriptorId;
        ^this.cmd("getDescriptor", [descriptorId.asString]);
    }

    loadOnnxModel { arg modelPath, descriptorIds;
        var csv = descriptorIds.asArray.collect(_.asString).join(",");
        ^this.cmd("loadOnnxModel", [modelPath.asString, csv]);
    }

    registerActionReceiver { arg receiver;
        ^this.cmd("registerActionReceiver", [receiver.asString]);
    }

    unregisterActionReceiver { arg receiver;
        ^this.cmd("unregisterActionReceiver", [receiver.asString]);
    }

    listen { arg receiver, action;
        var receiverString, oscName, address;
        receiverString = receiver.asString;
        oscName = ("OpenScofo_action_" ++ this.identityHash ++ "_" ++ receiverString).asSymbol;
        address = ("/OpenScofo/" ++ receiverString).asSymbol;

        this.unlisten(receiverString);
        actionOSCDefs[receiverString] = OSCdef(oscName, { |msg|
            action.value(msg[3..], msg);
        }, address, server.addr);
        this.registerActionReceiver(receiverString);
    }

    unlisten { arg receiver;
        var receiverString, oscdef;
        receiverString = receiver.asString;
        oscdef = actionOSCDefs[receiverString];
        if(oscdef.notNil) {
            oscdef.free;
            actionOSCDefs.removeAt(receiverString);
            this.unregisterActionReceiver(receiverString);
        };
    }

    registerEventHandler { arg action;
        eventAction = action;
        if(eventOSCDef.notNil) {
            eventOSCDef.free;
            eventOSCDef = nil;
        };
        if(action.notNil) {
            eventOSCDef = OSCdef(eventOSCDefName, { |msg|
                action.value(msg[3].asInteger, msg);
            }, '/oscofo/currentEvent', server.addr);
        };
    }

    free {
        if(eventOSCDef.notNil) {
            eventOSCDef.free;
            eventOSCDef = nil;
        };
        actionOSCDefs.keysValuesDo({ |receiver, oscdef|
            oscdef.free;
            this.unregisterActionReceiver(receiver);
        });
        actionOSCDefs.clear;
        if(synth.notNil) {
            synth.free;
            synth = nil;
        };
    }
}
