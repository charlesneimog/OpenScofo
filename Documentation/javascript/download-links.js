const assets_list = ["PureData", "Max", "Vamp", "SuperCollider", "CSound", "Python", "JavaScript", "C++"];

// GLOBAL variable to store release data
window.latestData = null;

// ─────────────────────────────────────
async function getLatestTagWithDetails() {
    if (latestData) return latestData; // return cached data if available

    const repo = "charlesneimog/OpenScofo";

    const tagRes = await fetch(`https://api.github.com/repos/${repo}/tags?per_page=1`);
    if (!tagRes.ok) throw new Error("Failed to fetch tags");

    const tags = await tagRes.json();
    if (tags.length === 0) throw new Error("No tags found");

    const tag = tags[0];

    const relRes = await fetch(`https://api.github.com/repos/${repo}/releases/tags/${tag.name}`);
    let release = null;
    if (relRes.ok) release = await relRes.json();

    latestData = {
        tag: tag.name,
        assets:
            release?.assets?.map((a) => ({
                name: a.name,
                url: a.browser_download_url,
            })) ?? [],
    };

    return latestData;
}

// ─────────────────────────────────────
function renderReleaseTable(data, assetName) {
    const container = document.querySelector(`release[interface="${assetName}"]`);
    if (!container) return;

    if (assetName == "JavaScript") {
        assetName = "Emscripten";
    }
    const filteredAssets = data.assets.filter((a) => a.name.includes(assetName));
    if (filteredAssets.length === 0) return;

    const table = document.createElement("table");
    table.style.width = "50%";
    table.style.borderCollapse = "collapse";
    table.style.textAlign = "center";
    table.style.tableLayout = "fixed";

    const thead = document.createElement("thead");
    thead.innerHTML = `
        <tr>
            <th style="width: 30%;">Version</th>
            <th style="width: 30%;">OS</th>
            <th style="width: 40%;">Download</th>
        </tr>
    `;
    table.appendChild(thead);

    const tbody = document.createElement("tbody");
    filteredAssets.forEach((a) => {
        const parts = a.name.replace(".zip", "").split("-");
        const version = data.tag;
        var os = parts[parts.length - 1] || "unknown";
        if (assetName == "Emscripten") {
            os = "WASM";
        }

        const tr = document.createElement("tr");

        const tdVersion = document.createElement("td");
        tdVersion.innerHTML = `<code>` + version + `</code>`;

        const tdOS = document.createElement("td");
        tdOS.textContent = os;

        const tdDownload = document.createElement("td");
        const link = document.createElement("a");
        link.href = a.url;
        link.textContent = "Download";
        link.target = "_blank";
        tdDownload.appendChild(link);

        tr.appendChild(tdVersion);
        tr.appendChild(tdOS);
        tr.appendChild(tdDownload);

        tbody.appendChild(tr);
    });
    table.appendChild(tbody);

    table.querySelectorAll("th, td").forEach((cell) => {
        cell.style.border = "1px solid #ccc";
        cell.style.padding = "8px 12px";
    });

    // wrap table in a centered div
    const wrapper = document.createElement("div");
    wrapper.style.display = "flex";
    wrapper.style.justifyContent = "center";
    wrapper.style.margin = "20px 0";
    wrapper.appendChild(table);

    table.querySelectorAll("th").forEach((th) => {
        th.style.width = "50%";
    });

    // replace original container with the wrapper
    container.replaceWith(wrapper);
}

// ─────────────────────────────────────
// Attach focus listeners to all release elements
window.lastOnFocus = null;

// ─────────────────────────────────────
document.addEventListener("click", async (event) => {
    assets_list.forEach((name) => {
        setTimeout(async () => {
            const el = document.querySelector(`release[interface="${name}"]`);
            if (el !== null) {
                if (window.lastOnFocus !== name) {
                    window.lastOnFocus = name;
                    if (window.latestData === null) {
                        window.latestData = await getLatestTagWithDetails();
                    }
                    renderReleaseTable(window.latestData, name);
                }
            }
        }, 50); // 50ms de delay
    });
});

// ─────────────────────────────────────
document.addEventListener("DOMContentLoaded", async () => {
    assets_list.forEach((name) => {
        setTimeout(async () => {
            const el = document.querySelector(`release[interface="${name}"]`);
            if (el !== null) {
                if (window.lastOnFocus !== name) {
                    window.lastOnFocus = name;
                    if (window.latestData === null) {
                        window.latestData = await getLatestTagWithDetails();
                    }
                    renderReleaseTable(window.latestData, name);
                }
            }
        }, 50); // 50ms delay
    });
});
