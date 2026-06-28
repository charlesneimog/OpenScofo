# OpenScofo for Neovim

Neovim support for `.scofo` files lives in the main OpenScofo repository. It
provides filetype detection, Tree-sitter parser registration, syntax highlights,
and Lua injections inside `LUA { ... }` blocks.

## Installation

Using `lazy.nvim`:

```lua
{
  "charlesneimog/OpenScofo",
  dependencies = { "nvim-treesitter/nvim-treesitter" },
  init = function(plugin)
    vim.opt.rtp:append(plugin.dir .. "/Sources/Language/nvim")
  end,
  config = function()
    require("openscofo").setup()
  end,
}
```

Then install the parser:

```vim
:TSInstall openscofo
```

Open a `.scofo` file and Neovim should use the `openscofo` filetype.

## Development

The parser grammar is stored in `Sources/Language/grammar.js`. The generated
parser source expected by `nvim-treesitter` is stored in
`Sources/Language/src/parser.c`.

When the grammar changes, regenerate the parser and commit the updated generated
files:

```sh
cd Sources/Language
tree-sitter generate
```

The query files used by Neovim are:

- `Sources/Language/nvim/queries/openscofo/highlights.scm`
- `Sources/Language/nvim/queries/openscofo/injections.scm`
