(() => {
    const base = new URL("../", document.currentScript.src);
    const asset = (path) => new URL(path, base).href;
    let runtime;
    let pd;
    let mounted;
    let playing = false;
    let selected = 0;

    async function isolate() {
        if (window.crossOriginIsolated) return true;
        if (!window.isSecureContext || !navigator.serviceWorker) {
            throw new Error("This player needs HTTPS and service worker support.");
        }
        const key = `openscofo-isolation:${base.pathname}`;
        if (sessionStorage.getItem(key)) {
            sessionStorage.removeItem(key);
            throw new Error("Audio could not be enabled. Allow service workers and reload the page.");
        }
        const registration = await navigator.serviceWorker.register(asset("pd4web-isolation.js"));
        await navigator.serviceWorker.ready;
        if (!navigator.serviceWorker.controller) {
            await new Promise((resolve) => {
                navigator.serviceWorker.addEventListener("controllerchange", resolve, { once: true });
                if (navigator.serviceWorker.controller) resolve();
            });
        }
        // Credentialless mode also permits the documentation's public CDN assets.
        registration.active.postMessage({ type: "coepCredentialless", value: true });
        sessionStorage.setItem(key, "1");
        window.location.reload();
        return false;
    }

    async function loadRuntime() {
        if (!await isolate()) return null;
        sessionStorage.removeItem(`openscofo-isolation:${base.pathname}`);
        await new Promise((resolve, reject) => {
            const script = document.createElement("script");
            script.src = asset("pieces/pd4web.js");
            script.onload = resolve;
            script.onerror = () => reject(new Error("Could not download the audio player."));
            document.head.append(script);
        });
        const module = await Pd4WebModule({ locateFile: (file) => asset(`pieces/${file}`) });
        pd = new module.Pd4Web();
        pd.openPatch("index.pd", {
            renderGui: false,
            hideInitError: true,
            channelCountIn: 0,
            channelCountOut: 2,
            sampleRate: 48000,
            requestMidi: false,
        });
        // The patch sends floats with [s scoreevent], rather than selector messages.
        pd.onFloatReceived("scoreevent", (_receiver, event) => {
            if (!mounted || !playing) return;
            const score = mounted.score;
            if (event === 0) {
                score.querySelectorAll("[data-followed]").forEach((note) => note.removeAttribute("data-followed"));
                return;
            }
            if (!Number.isInteger(event) || event < 0) return;
            const notes = score.querySelectorAll(`[openscofo-event-number="${event}"]`);
            notes.forEach((note) => note.setAttribute("data-followed", ""));
            if (notes.length) {
                const bounds = score.getBoundingClientRect();
                const note = notes[0].getBoundingClientRect();
                // Scroll only the score, leaving the surrounding documentation in place.
                if (note.top < bounds.top || note.bottom > bounds.bottom) {
                    score.scrollTo({ top: score.scrollTop + note.top - bounds.top - bounds.height / 2, behavior: "smooth" });
                }
            }
        });
        return pd;
    }

    function updateButtons(view) {
        view.buttons.forEach((button) => {
            const active = playing && Number(button.dataset.piece) === selected;
            button.setAttribute("aria-pressed", String(active));
            button.textContent = `${active ? "Pause" : "Play"} Miniatura ${button.dataset.piece}`;
        });
    }

    function showPiece(view, scores, piece) {
        view.score.innerHTML = scores[piece - 1];
        view.score.scrollTop = 0;
        view.code.forEach((panel) => {
            panel.hidden = Number(panel.dataset.demoCode) !== piece;
            panel.scrollTop = 0;
            panel.scrollLeft = 0;
        });
    }

    async function mount(root) {
        const view = {
            root,
            score: root.querySelector("[data-demo-score]"),
            code: root.querySelectorAll("[data-demo-code]"),
            status: root.querySelector("[data-demo-status]"),
            buttons: root.querySelectorAll("[data-piece]"),
        };
        mounted = view;
        try {
            const [engine, ...scores] = await Promise.all([
                runtime ??= loadRuntime(),
                ...[1, 2].map(async (piece) => {
                    const response = await fetch(asset(`pieces/score${piece}.svg`));
                    if (!response.ok) throw new Error(`Could not load score ${piece}.`);
                    return response.text();
                }),
            ]);
            if (!engine || mounted !== view || !root.isConnected) return;
            selected = 0;
            showPiece(view, scores, 1);
            view.status.textContent = "Choose a piece. Click its button again to pause.";
            view.buttons.forEach((button) => {
                button.disabled = false;
                button.addEventListener("click", () => {
                    try {
                        const piece = Number(button.dataset.piece);
                        if (selected === piece) {
                            pd.toggleAudio();
                            playing = !playing;
                        } else {
                            // toggleAudio initializes Web Audio on the first user click.
                            if (!playing) pd.toggleAudio();
                            playing = true;
                            selected = piece;
                            showPiece(view, scores, piece);
                            pd.sendFloat("score", piece);
                            pd.sendFloat("audio", piece);
                        }
                        updateButtons(view);
                        view.status.textContent = `${playing ? "Playing" : "Paused"}: Miniatura ${piece}`;
                    } catch (error) {
                        view.status.textContent = "Could not start audio. Reload the page and try again.";
                        console.error(error);
                    }
                });
            });
        } catch (error) {
            if (mounted === view) view.status.textContent = `Player unavailable. ${error.message}`;
            console.error(error);
        }
    }

    function initialize() {
        const root = document.querySelector("[data-score-demo]");
        if (mounted?.root === root) return;
        if (playing) {
            pd.toggleAudio();
            playing = false;
        }
        mounted = null;
        if (root) mount(root);
    }

    // Support the documentation theme's instant navigation, using one audio engine.
    if (typeof document$ !== "undefined") document$.subscribe(initialize);
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", initialize, { once: true });
    } else {
        initialize();
    }
})();
