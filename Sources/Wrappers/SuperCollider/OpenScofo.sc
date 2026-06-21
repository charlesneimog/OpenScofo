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
    var <inBus, <unitIndex, <eventOSCDef, <namespace;
    var eventAction, synthDefName, eventOSCDefName, actionOSCDefs;

    *new { arg scorePath, inBus = 0, sampleRate, server, target, addAction = \addToTail,
        eventNotifications = true, eventAction, namespace = "openscofo";
        ^super.new.init(scorePath, inBus, sampleRate, server, target, addAction, eventNotifications, eventAction, namespace)
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

    init { arg scorePathArg, inBusArg = 0, sampleRateArg, serverArg, targetArg,
        addActionArg = \addToTail, eventNotifications = true, eventActionArg, namespaceArg = "openscofo";
        var sr;

        server = serverArg ? Server.default;
        scorePath = scorePathArg;
        inBus = inBusArg;
        sampleRate = sampleRateArg ? server.sampleRate;
        target = targetArg ? server.defaultGroup;
        addAction = addActionArg;
        namespace = this.normalizeNamespace(namespaceArg);
        eventAction = eventActionArg;
        actionOSCDefs = IdentityDictionary.new;

        sr = sampleRate;
        synthDefName = ("OpenScofo_" ++ this.identityHash).asSymbol;
        eventOSCDefName = ("OpenScofo_event_" ++ this.identityHash).asSymbol;

        SynthDef(synthDefName, { arg inBus = 0, sampleRate = sr;
            var in;
            in = In.ar(inBus, 1);
            OpenScofoUGen.ar(in, sampleRate);
        }).add;

        server.sync;
        synth = Synth(synthDefName, [
            \inBus, inBus,
            \sampleRate, sampleRate
        ], target: target, addAction: addAction);
        server.sync;

        unitIndex = this.findUnitIndex;

        this.setNamespace(namespace);
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

    normalizeNamespace { arg namespaceArg;
        var value = namespaceArg.asString;

        while({ value.notEmpty and: { value[0] == $/ } }) {
            value = if(value.size > 1, { value.copyRange(1, value.size - 1) }, { "" });
        };
        while({ value.notEmpty and: { value[value.size - 1] == $/ } }) {
            value = if(value.size > 1, { value.copyRange(0, value.size - 2) }, { "" });
        };

        if(value.isEmpty) {
            ^"openscofo"
        };
        ^value
    }

    actionAddress { arg receiver;
        ^("/" ++ namespace ++ "/" ++ receiver.asString).asSymbol
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

    setNamespace { arg namespaceArg = "openscofo";
        namespace = this.normalizeNamespace(namespaceArg);
        ^this.cmd("setNamespace", [namespace]);
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

    decodeActionPayload { arg payload;
        var magic = -900001.0, args, argCount, index = 2, type, size;

        if(payload.isNil or: { payload.size == 0 }) {
            ^[]
        };

        if(payload[0] != magic) {
            ^payload
        };

        argCount = payload[1].asInteger;
        args = Array.new(argCount);

        argCount.do({
            type = payload[index].asInteger;
            index = index + 1;

            if(type == 0) {
                args = args.add(payload[index]);
                index = index + 1;
            } {
                if(type == 1) {
                    size = payload[index].asInteger;
                    index = index + 1;
                    args = args.add(String.fill(size, { |i|
                        payload[index + i].asInteger.asAscii
                    }));
                    index = index + size;
                };
            };
        });

        ^args
    }

    listen { arg receiver, action;
        var receiverString, oscName, address;
        receiverString = receiver.asString;
        oscName = ("OpenScofo_action_" ++ this.identityHash ++ "_" ++ receiverString).asSymbol;
        address = this.actionAddress(receiverString);

        this.unlisten(receiverString);
        actionOSCDefs[receiverString] = OSCdef(oscName, { |msg|
            action.value(this.decodeActionPayload(msg[3..]), msg);
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
            }, this.actionAddress("currentEvent"), server.addr);
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
