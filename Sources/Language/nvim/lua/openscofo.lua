local M = {}

-- ─────────────────────────────────────
local module_file = debug.getinfo(1, "S").source:sub(2)
local language_root = vim.fs.dirname(vim.fs.dirname(vim.fs.dirname(module_file)))

-- ─────────────────────────────────────
local function ensure_parser()
	local parser = vim.fs.joinpath(language_root, "parser.so")
	local source = vim.fs.joinpath(language_root, "src", "parser.c")

	local parser_stat = vim.uv.fs_stat(parser)
	local source_stat = vim.uv.fs_stat(source)

	if parser_stat and source_stat and parser_stat.mtime.sec >= source_stat.mtime.sec then
		return parser
	end

	local temporary = parser .. ".tmp"
	local result = vim.system({
		"cc",
		"-shared",
		"-fPIC",
		"-Isrc",
		"-o",
		temporary,
		"src/parser.c",
	}, {
		cwd = language_root,
		text = true,
	}):wait()

	if result.code ~= 0 then
		error("OpenScofo: could not compile parser:\n" .. (result.stderr or ""))
	end

	local ok, err = vim.uv.fs_rename(temporary, parser)
	assert(ok, "OpenScofo: could not install parser: " .. (err or ""))

	return parser
end

-- ─────────────────────────────────────
function M.setup()
	local parser_path = ensure_parser()

	local ok, err = vim.treesitter.language.add("openscofo", {
		path = parser_path,
	})
	assert(ok, err)

	vim.treesitter.language.register("openscofo", "openscofo")

	vim.filetype.add({
		extension = {
			scofo = "openscofo",
		},
	})

	local group = vim.api.nvim_create_augroup("OpenScofoTreesitter", { clear = true })

	vim.api.nvim_create_autocmd("FileType", {
		group = group,
		pattern = "openscofo",
		callback = function(ev)
			vim.treesitter.start(ev.buf, "openscofo")
		end,
	})
end

return M
