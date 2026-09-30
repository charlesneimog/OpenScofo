---
hide:
  - navigation
  - toc
tags:
  - Overview
---

<style>
  .md-typeset h1,
  .md-content__button {
    display: none;
  }
</style>

# OpenScofo

<p align="center" markdown>
  ![OpenScofo logo](./assets/logo.svg#only-light){ width="15%" }
  ![OpenScofo logo](./assets/logo-dark.svg#only-dark){ width="15%" }
</p>

<h4 latex="false" align="center"><i>Score following for contemporary music.</i></h4>

---

OpenScofo exists so live electronics can stay aligned with a human performer. It follows a notated score in real time and triggers computer actions at musical events. Basically, based on a tradicional notated score, you create an electronic score as showed below.

---

<div class="score-demo" data-score-demo markdown="1">
  <h2 align="center">Listen and follow the score</h2>
  <p>This demo plays a microphone recording of Cassia Carrascoza’s performance. OpenScofo processes the recording live as it plays and follows the score in real time. Everything runs in your browser using <a href="https://charlesneimog.github.io/pd4web/">pd4web</a>, which runs Pure Data on the web.</p>

  --- 

  <div class="score-demo__controls" role="group" aria-label="Play a piece">
    <button type="button" class="md-button" data-piece="1" aria-pressed="false" disabled>Play Miniatura 1</button>
    <button type="button" class="md-button" data-piece="2" aria-pressed="false" disabled>Play Miniatura 2</button>
  </div>

  <p data-demo-status role="status">Loading the interactive player…</p>

  <div class="score-demo__layout" markdown="1">

  <div class="score-demo__panel">

  <p><strong>Music score</strong></p>
  <div class="score-demo__score" data-demo-score role="region" aria-label="Selected music score" tabindex="0"></div>
</div>

<div class="score-demo__panel" markdown="1">
  <p><strong>OpenScofo code</strong></p>
  <div class="score-demo__code" data-demo-code="1" role="region" aria-label="Miniatura 1 OpenScofo code" tabindex="0" markdown="1">

  ```openscofo
  --8<-- "Documentation/pieces/miniatura1.scofo"
  ```
</div>

<div class="score-demo__code" data-demo-code="2" role="region" aria-label="Miniatura 2 OpenScofo code" tabindex="0" hidden markdown="1">
  ```openscofo
  --8<-- "Documentation/pieces/miniatura2.scofo"
  ```
</div>
</div>
</div>
<noscript>Enable JavaScript to play the pieces and follow the score.</noscript>
</div>

---
