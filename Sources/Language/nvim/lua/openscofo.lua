local M = {}

local parser_url = "https://github.com/charlesneimog/OpenScofo"

local parser_config = {
  install_info = {
    url = parser_url,
    branch = "main",
    location = "Sources/Language",
    queries = "Sources/Language/nvim/queries/openscofo",

    -- Used by older nvim-treesitter releases.
    files = { "src/parser.c" },
    requires_generate_from_grammar = false,
  },
  filetype = "openscofo",
  tier = 3,
}

local function package_root()
  local source = debug.getinfo(1, "S").source:sub(2)
  return vim.fn.fnamemodify(source, ":p:h:h:h:h:h")
end

local function register_filetype()
  vim.filetype.add({
    extension = {
      scofo = "openscofo",
    },
  })

  if vim.treesitter and vim.treesitter.language and vim.treesitter.language.register then
    vim.treesitter.language.register("openscofo", { "openscofo" })
  end
end

local function register_treesitter_parser()
  local ok, parsers = pcall(require, "nvim-treesitter.parsers")
  if not ok then
    return
  end

  if type(parsers.get_parser_configs) == "function" then
    parsers.get_parser_configs().openscofo = parser_config
    return
  end

  local config = vim.deepcopy(parser_config)
  config.install_info.path = package_root()
  parsers.openscofo = config
end

local function register_treesitter_autocmd()
  local group = vim.api.nvim_create_augroup("OpenScofoTreesitter", { clear = true })

  vim.api.nvim_create_autocmd("User", {
    group = group,
    pattern = "TSUpdate",
    callback = register_treesitter_parser,
  })

  vim.api.nvim_create_autocmd("FileType", {
    group = group,
    pattern = "openscofo",
    callback = function(args)
      pcall(vim.treesitter.start, args.buf, "openscofo")
    end,
  })
end

local function parser_installed()
  local ok, treesitter = pcall(require, "nvim-treesitter")
  if ok and type(treesitter.get_installed) == "function" then
    return vim.list_contains(treesitter.get_installed("parsers"), "openscofo")
  end

  return #vim.api.nvim_get_runtime_file("parser/openscofo.*", true) > 0
end

local function install_parser()
  if parser_installed() then
    return
  end

  local ok, treesitter = pcall(require, "nvim-treesitter")
  if ok and type(treesitter.install) == "function" then
    treesitter.install({ "openscofo" })
    return
  end

  if vim.fn.exists(":TSInstall") == 2 then
    vim.cmd("silent! TSInstall openscofo")
  end
end

function M.setup(opts)
  opts = opts or {}

  register_filetype()
  register_treesitter_parser()
  register_treesitter_autocmd()

  if opts.auto_install ~= false then
    vim.schedule(install_parser)
  end
end

return M
