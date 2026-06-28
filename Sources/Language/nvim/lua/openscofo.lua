local M = {}

local parser_url = "https://github.com/charlesneimog/OpenScofo"

local function register_filetype()
  vim.filetype.add({
    extension = {
      scofo = "openscofo",
    },
  })
end

local function register_treesitter_parser()
  local ok, parsers = pcall(require, "nvim-treesitter.parsers")
  if not ok then
    return
  end

  local parser_config = parsers.get_parser_configs()

  parser_config.openscofo = {
    install_info = {
      url = parser_url,
      files = { "src/parser.c" },
      branch = "main",
      location = "Sources/Language",
      requires_generate_from_grammar = false,
    },
    filetype = "openscofo",
  }
end

function M.setup()
  register_filetype()
  register_treesitter_parser()
end

M.setup()

return M
