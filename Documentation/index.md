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

<div class="grid cards" markdown>
- __Music Score__

    ![Minimal score](./assets/events/minimal-score.png)

- __OpenScofo Score__

    ```openscofo
    NOTE C4 1
    NOTE D4 1
        sendto activated_computer_processing [1]
    ```
</div>

--- 

