async function getLatestDownloads() {
  const res = await fetch(
    "https://api.github.com/repos/charlesneimog/OpenScofo/releases/latest"
  );
  const data = await res.json();
  console.log(data);

  // all uploaded assets
  const downloads = data.assets.map(a => ({
    name: a.name,
    url: a.browser_download_url
  }));

  return downloads;
}

getLatestDownloads().then(console.log);
