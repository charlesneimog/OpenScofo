const assets_list = ["PureData", "Max", "Vamp", "SuperCollider", "CSound", "Python", "JavaScript", "C++"];

// GLOBAL variable to store release data
window.latestData = null;

// ─────────────────────────────────────
function getIcons(os_name) {
    let iconSrc = "";
    if (os_name == "Windows") {
        iconSrc =
            "https://raw.githubusercontent.com/squidfunk/mkdocs-material/master/material/templates/.icons/fontawesome/brands/windows.svg";
    } else if (os_name == "Linux") {
        iconSrc =
            "https://raw.githubusercontent.com/squidfunk/mkdocs-material/master/material/templates/.icons/simple/linux.svg";
    } else if (os_name == "macOS") {
        iconSrc =
            "https://raw.githubusercontent.com/squidfunk/mkdocs-material/master/material/templates/.icons/simple/apple.svg";
    } else if (os_name == "WASM") {
        iconSrc =
            "https://raw.githubusercontent.com/squidfunk/mkdocs-material/master/material/templates/.icons/simple/webassembly.svg";
    }

    return `<div style="display: flex; align-items: center; gap: 0.5em;">
                <span class="twemoji"><img src="${iconSrc}" ></span>
                <p style="margin: 0;">${os_name}</p>
            </div>`;
}
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

    let filteredAssets;
    if (assetName === "All") {
        filteredAssets = data.assets; // no filter
    } else if (assetName === "JavaScript") {
        filteredAssets = data.assets.filter((a) => a.name.includes("Emscripten"));
    } else {
        filteredAssets = data.assets.filter((a) => a.name.includes(assetName));
    }

    if (filteredAssets.length === 0) return;

    const table = document.createElement("table");
    table.style.width = "80%";
    table.style.borderCollapse = "collapse";
    table.style.textAlign = "center";
    table.style.tableLayout = "fixed";

    const thead = document.createElement("thead");

    if (assetName === "All") {
        thead.innerHTML = `
        <tr>
            <th style="width: 20%;">Enviroment</th>
            <th style="width: 5%;">Version</th>
            <th style="width: 25%;">OS</th>
            <th style="width: 20%;">Download</th>
        </tr>
      `;
    } else {
        thead.innerHTML = `
        <tr>
            <th style="width: 30%;">Version</th>
            <th style="width: 30%;">OS</th>
            <th style="width: 40%;">Download</th>
        </tr>
      `;
    }
    table.appendChild(thead);

    const tbody = document.createElement("tbody");
    filteredAssets.forEach((a) => {
        const parts = a.name.replace(".zip", "").split("-");
        const version = data.tag;
        let os = parts[parts.length - 1] || "unknown";
        if (assetName === "JavaScript" || assetName === "Emscripten") os = "WASM";
        if (parts[0] === "Emscripten") {
            os = "WASM";
        }

        const tr = document.createElement("tr");
        const tdVersion = document.createElement("td");
        tdVersion.innerHTML = `<code>${version}</code>`;

        const tdOS = document.createElement("td");
        // tdOS.textContent = os;
        tdOS.innerHTML = getIcons(os);

        const tdDownload = document.createElement("td");
        const link = document.createElement("a");
        link.href = a.url;
        link.textContent = "Download";
        link.target = "_blank";
        tdDownload.appendChild(link);

        if (assetName == "All") {
            const tdSystem = document.createElement("td");
            tdSystem.innerHTML = parts[0];
            tr.appendChild(tdSystem);
        }

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

    table.querySelectorAll("th").forEach((th) => {
        th.style.width = "50%";
    });

    const wrapper = document.createElement("div");
    wrapper.style.display = "flex";
    wrapper.style.justifyContent = "center";
    wrapper.style.margin = "20px 0";
    wrapper.appendChild(table);

    container.replaceWith(wrapper);
}

// ─────────────────────────────────────
// Attach focus listeners to all release elements
window.lastOnFocus = null;

// ─────────────────────────────────────
document.addEventListener("click", async (event) => {
    assets_list.concat("All").forEach((name) => {
        setTimeout(async () => {
            const el = document.querySelector(`release[interface="${name}"]`);
            if (el && window.lastOnFocus !== name) {
                window.lastOnFocus = name;
                if (!window.latestData) window.latestData = await getLatestTagWithDetails();
                renderReleaseTable(window.latestData, name);
            }
        }, 250);
    });
});

// ─────────────────────────────────────
document.addEventListener("DOMContentLoaded", async () => {
    assets_list.concat("All").forEach((name) => {
        setTimeout(async () => {
            const el = document.querySelector(`release[interface="${name}"]`);
            if (el && window.lastOnFocus !== name) {
                window.lastOnFocus = name;
                if (!window.latestData) window.latestData = await getLatestTagWithDetails();
                renderReleaseTable(window.latestData, name);
            }
        }, 250);
    });
});
