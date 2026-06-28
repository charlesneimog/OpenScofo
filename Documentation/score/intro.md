# Introduction

`OpenScofo` uses a textual score where the composer must define the musical `EVENTS` and its consequences, what is called `ACTIONS`. In addition to that, `OpenScofo` has also some `CONFIG` keywords and the possibility to use `LUA` language to create interactive `ACTIONS`. In this section is explained all these concepts.

---
<div class="grid cards" style="font-weigth:bold" markdown>

-   :fontawesome-solid-1: [__How to Configure the Score?__](config.md)

-   :fontawesome-solid-2: [__How to Set New Events?__](events.md)

</div>

<div class="grid cards" style="font-weigth:bold" markdown>
-   :fontawesome-solid-3: [__How to add Actions?__](actions.md)

-   :fontawesome-solid-4: [__How to use Lua?__](lua.md)

</div>

---


!!! tip "Language Parse for VSCode and Neovim"
    There is a language parse for VSCode (search for OpenScofo [extension](https://marketplace.visualstudio.com/items?itemName=charlesneimog.openscofo-language-parse){:target = "_blank"}) and for [Neovim](https://github.com/charlesneimog/OpenScofo/tree/main/Sources/Language/nvim){:target + "_blank"}.

After Configuration, the Score will look like this:

<p align="center">
    <img style="width: 80%; border-radius: 5px" src="https://charlesneimog.github.io/OpenScofo/assets/oscofo-code.png">
</p>

## Filename, Extension, and Syntax Highlighting

An `OpenScofo` score is a plain text file. Although any text file extension can be used, the recommended extension is `.scofo`.

To improve readability and editing, support for `.scofo` files is already available on multiple platforms:

- **VS Code**: Install the [OpenScofo Language Parse extension](https://marketplace.visualstudio.com/items?itemName=charlesneimog.openscofo-language-parse){:target="_blank"}.
- **Neovim**: See the [Neovim configuration](https://github.com/charlesneimog/OpenScofo/tree/main/Sources/Language/nvim){:target="_blank"}.
- **Web browser**: Use the [OpenScofo Online Editor](https://charlesneimog.github.io/OpenScofo/Editor/){:target="_blank"} for editing and experimenting with `.scofo` files directly in your browser.

Support for additional editors is under development.
