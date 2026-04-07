<p align="center">
  <h1 align="center">OpenScofo</h1>
  <h4 align="center">OpenScofo: OpenScore Follower</h4>
</p>

`OpenScofo` is an open-source project designed to provide score following capabilities for contemporary music applications. Originally developed as a PureData (Pd) object, OpenScofo has now been expanded into a versatile C++ library that integrates seamlessly with multiple environments, including a Max Object, Python package, and others comming. Currently under development, OpenScofo is already functional and serves as a valuable tool for researchers and musicians.

## Download

* Check the [Releases](https://github.com/charlesneimog/OpenScofo/releases/latest) Page.

## Examples

* [Cânticos de Silício](https://charlesneimog.github.io/Canticos-de-Silicio-I/) by Charles K. Neimog;
* [BWV-1013](https://charlesneimog.github.io/pd4web/tests/OpenScofo) by Bach;

## Goal

The aim of *OpenScofo* is to provide a straightforward and accessible tool for real-time score following. By keeping the software lightweight, it can run seamlessly on the web through the [pd4web](https://charlesneimog.github.io/pd4web/) platform, thanks to the ability to use PureData directly in web browsers. With _pd4web_ and _OpenScofo_ will be possible to use the software in rehearsals with just a single click, eliminating the need for external libraries, compatibility issues, or complex installations -- ultimately facilitating the sharing and performance of contemporary music.

## Collaboration and Contribution

I invite composers, researchers and developers to contribute to the *OpenScofo* project. Not just with code, but with theory, math, etc. I am trying to make OpenScofo acessible via a Python implementation, to test it should be easy. By sharing the source code, I am trying to provide access to the theories and mathematical formulas that drive the software, all this come from the amazing research work of Arshia Cont and Philippe Cuvillier at IRCAM. 

## Technical Foundations

*OpenScofo* uses several concepts developed by many researches (with focus on the research of Cont and Cuvillier).

* **Pitch Comparison**: Utilizes the Kullback-Leibler (KL) Divergence method for pitch comparison as presented by Christopher Raphael (2006), Arshia Cont in 2008 and 2010.
* **Rhythm Synchronization**: Integrates theories of rhythm synchronization developed by Edward Large and Mari Riess Jones (1999) and Edward Large and Caroline Palmer (2002), as presented for Cont (2010).
* **Forward Algorithm**: For now, *OpenScofo* uses the equation presented by Arshia Cont (2010) and developed by Yann Guédon (2005).
* **Score Language**: Based on the `scofo` (by Miller Puckette) and `antescofo~` (by Arshia Cont, Philippe Cuvillier, and others) language.

## Building

### Requirements

* On Windows, you need `mingw64`.

#### Optional

* treesitter (`npm install tree-sitter`) (If you want to change/update score syntax).
* nanobind: To build Python package.
* PureData (optional): To build the Pd Object.

#### Build Options

* `OSCOFO_BUILD_ALL` (ON/OFF): Build all OpenScofo modules, including Python, Pd, and Max. Default: ON.
* `OSCOFO_BUILD_PD_OBJECT` (ON/OFF): Build the Pure Data (Pd) object. Default: OFF.
* `OSCOFO_BUILD_PY_MODULE` (ON/OFF): Build the Python module. Default: OFF.
* `OSCOFO_BUILD_MAX_OBJECT` (ON/OFF): Build the Max object. Default: OFF.
* `OSCOFO_BUILD_CSOUND_PLUGIN` (ON/OFF): Build the Csound plugin. Default: OFF.
* `OSCOFO_BUILD_TESTS` (ON/OFF): Build test suite. Default: OFF.
* `OSCOFO_BUILD_WITH_LUA` (ON/OFF): Build Lua module embedded in OpenScofo. Default: ON.
* `UPDATE_OSCOFO_LANGUAGE` (ON/OFF): Update OpenScofo language grammar (`grammar.js`). Default: ON.


#### Building Steps

``` bash
git clone https://github.com/charlesneimog/OpenScofo
cmake . -B build
cmake --build build
```

## GPLv3 Licensing Notice

OpenScofo is distributed under the GNU General Public License version 3 (GPLv3).

The GPLv3 is a free and open source license. Any software that incorporates, links to (including dynamically), or forms a combined work with OpenScofo and is distributed must also be licensed under the GPLv3, with its complete corresponding source code made available.

Non-compliance automatically terminates the license rights granted under the GPLv3 and may constitute copyright infringement.

Full license text: [https://www.gnu.org/licenses/gpl-3.0.html](https://www.gnu.org/licenses/gpl-3.0.html)
