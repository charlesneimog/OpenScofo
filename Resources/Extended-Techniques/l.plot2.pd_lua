local lplot = pd.Class:new():register("l.plot2")

-- ─────────────────────────────────────
function lplot:initialize(_, args)
	self.inlets = 1
	self.args = args

	self.tick_ms = 50
	self.tick_count = 0
	self.dx_per_tick = 2
	self.background = { 225, 225, 225 }

	if #args == 2 then
		self.width = args[1]
		self.height = args[2]
	else
		self.width = 400
		self.height = 200
	end

	self.draws = {}
	self.color = { 255, 0, 0 }
	self.background = { 245, 245, 245 }
	self:set_size(self.width, self.height)

	self.draws = {}

	return true
end

-- ─────────────────────────────────────
function lplot:postinitialize()
	self.clock = pd.Clock:new():register(self, "tick")
	self.clock:delay(self.tick_ms)
end

-- ─────────────────────────────────────
function lplot:tick()
	self.tick_count = self.tick_count + 1
	self:repaint()
	if self.clock then
		self.clock:delay(self.tick_ms)
	end
end
-- ─────────────────────────────────────
function lplot:in_1(sel, args)
	if #args > 1 then
		self:error("Just possible to render float")
		return
	end

	local val = tonumber(args[1]) or 0
	val = math.max(0, math.min(1, val))
	local isbang = val == "bang"

	-- Ensure selector exists
	if not self.draws[sel] then
		self.draws[sel] = {
			color = { math.random(0, 255), math.random(0, 255), math.random(0, 255) },
			points = {},
			count_since_last_draw = 0,
		}
	end

	local draw = self.draws[sel]

	if isbang then
		-- Determine the maximum length among all other draws
		local max_len = 0
		for _, d in pairs(self.draws) do
			if #d.points > max_len then
				max_len = #d.points
			end
		end

		-- Fill missing points with 0 so the bang will appear at the end
		while #draw.points < max_len do
			table.insert(draw.points, 0)
		end

		-- Append 1 for the bang
		table.insert(draw.points, 1)
	else
		-- Number input: append value normally
		local numval = tonumber(val) or 0
		table.insert(draw.points, numval)
	end

	draw.count_since_last_draw = draw.count_since_last_draw + 1

	-- Discard oldest if exceeding max points
	local max_points = math.floor(self.width / self.dx_per_tick)
	if #draw.points > max_points then
		table.remove(draw.points, 1)
	end

	-- Align all arrays to the same length
	local max_len = 0
	for _, d in pairs(self.draws) do
		if #d.points > max_len then
			max_len = #d.points
		end
	end

	for _, d in pairs(self.draws) do
		local last_val = d.points[#d.points] or 0
		while #d.points < max_len do
			table.insert(d.points, last_val)
		end
	end
end

-- ─────────────────────────────────────
local function hsv_to_rgb(h, s, v)
	local c = v * s
	local x = c * (1 - math.abs((h / 60) % 2 - 1))
	local m = v - c
	local r, g, b
	if h < 60 then
		r, g, b = c, x, 0
	elseif h < 120 then
		r, g, b = x, c, 0
	elseif h < 180 then
		r, g, b = 0, c, x
	elseif h < 240 then
		r, g, b = 0, x, c
	elseif h < 300 then
		r, g, b = x, 0, c
	else
		r, g, b = c, 0, x
	end
	return { math.floor((r + m) * 255), math.floor((g + m) * 255), math.floor((b + m) * 255) }
end

-- ─────────────────────────────────────
function lplot:string_to_color(sel)
	local hash = 0
	for i = 1, #sel do
		hash = (hash * 31 + sel:byte(i)) % 360 -- keep hash in 0..359
	end
	local hue = hash
	local sat = 0.8 -- fixed saturation
	local val = 0.9 -- fixed brightness
	return hsv_to_rgb(hue, sat, val)
end

-- ─────────────────────────────────────
function lplot:paint(g)
	g:set_color(table.unpack(self.background))
	g:fill_all()

	local margin = 5
	local line_height = 8 -- adjust font size/spacing as needed
	local y_text = margin

	for sel, draw in pairs(self.draws) do
		local col = draw.color
		g:set_color(table.unpack(col))

		local name = sel:sub(1, 8) -- truncate if longer
		name = name .. string.rep(" ", 8 - #name) -- pad if shorter
		local text_x = margin
		g:draw_text(name, text_x, y_text, 200, 3) -- width nil, fontsize 10
		y_text = y_text + line_height

		-- Fill missing value if nothing added
		if draw.count_since_last_draw == 0 then
			table.insert(draw.points, 0)
		end

		local x = self.width - margin
		local prev_x, prev_y = nil, nil
		draw.count_since_last_draw = 0

		for i = #draw.points, 1, -1 do
			local y = margin + (1 - draw.points[i]) * (self.height - 2 * margin)

			if prev_x then
				g:draw_line(x, y, prev_x, prev_y, 1)
			end

			prev_x, prev_y = x, y
			x = x - self.dx_per_tick
			if x < margin then
				break
			end
		end
	end
end

-- ─────────────────────────────────────
function lplot:in_1_reload()
	if self.clock then
		pcall(function()
			self.clock:destruct()
		end)
		self.clock = nil
	end
	self:dofilex(self._scriptname)
	self:initialize("l.plot", self.args)
	self:postinitialize()
end
