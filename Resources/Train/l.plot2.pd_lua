local lplot = pd.Class:new():register("l.plot2")

-- ─────────────────────────────────────
function lplot:initialize(_, args)
	self.inlets = 1
	self.args = args

	self.tick_ms = 50
	self.tick_count = 0
	self.dx_per_tick = 2

	if #args == 2 then
		self.width = args[1]
		self.height = args[2]
	else
		self.width = 400
		self.height = 200
	end

	self.background = { 245, 245, 245 }
	self:set_size(self.width, self.height)

	self.draws = {}

	-- color system
	self.color_index = 0
	self.BASE_COUNT = 6
	self.base_colors = {}
	self.current_max_category_name = ""
	self.current_max_category_value = 0

	for i = 0, self.BASE_COUNT - 1 do
		local hue = (i / self.BASE_COUNT) * 360
		local val = (i % 2 == 0) and 0.9 or 0.7
		self.base_colors[i + 1] = self:hsv_to_rgb(hue, 0.8, val)
	end

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
function lplot:properties(p)
	local i = 0
	p:new_frame("Colors", #self.draws)

	for k, v in pairs(self.draws) do
		local method = "update_color_" .. i
		self[method] = function(obj, color)
			obj.draws[k].color = color[1]
		end
		p:add_color(k, method, v.color)
		i = i + 1
	end
end

-- ─────────────────────────────────────
function lplot:get_next_color()
	self.color_index = self.color_index + 1

	if self.color_index <= self.BASE_COUNT then
		return self.base_colors[self.color_index]
	end

	-- generate near color
	local base_idx = ((self.color_index - 1) % self.BASE_COUNT) + 1
	local base_hue = ((base_idx - 1) / self.BASE_COUNT) * 360

	local sat = 0.6 + 0.3 * math.random()
	local val = 0.7 + 0.25 * math.random()

	return self:hsv_to_rgb(base_hue, sat, val)
end

-- ─────────────────────────────────────
function lplot:in_1(sel, args)
	if #args > 1 then
		self:error("Just possible to render float")
		return
	end

	local raw = args[1]
	local isbang = raw == "bang"

	if raw > self.current_max_category_value then
		self.current_max_category_name = sel
		self.current_max_category_value = raw
	end

	local val = tonumber(raw) or 0
	val = math.max(0, math.min(1, val))

	-- Ensure selector exists
	if not self.draws[sel] then
		self.draws[sel] = {
			color = self:get_next_color(),
			points = {},
			count_since_last_draw = 0,
		}
	end

	local draw = self.draws[sel]

	if isbang then
		local max_len = 0
		for _, d in pairs(self.draws) do
			if #d.points > max_len then
				max_len = #d.points
			end
		end

		while #draw.points < max_len do
			table.insert(draw.points, 0)
		end

		table.insert(draw.points, 1)
	else
		table.insert(draw.points, val)
	end

	draw.count_since_last_draw = draw.count_since_last_draw + 1

	local max_points = math.floor(self.width / self.dx_per_tick)
	if #draw.points > max_points then
		table.remove(draw.points, 1)
	end

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
function lplot:hsv_to_rgb(h, s, v)
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

	return {
		math.floor((r + m) * 255),
		math.floor((g + m) * 255),
		math.floor((b + m) * 255),
	}
end

-- ─────────────────────────────────────
function lplot:paint(g)
	g:set_color(table.unpack(self.background))
	g:fill_all()

	local margin = 5
	local line_height = 8
	local y_text = margin

	local w, h = self:get_size()

	for sel, draw in pairs(self.draws) do
		local col = draw.color
		g:set_color(table.unpack(col))

		local name = sel:sub(1, 8)
		name = name .. string.rep(" ", 8 - #name)
		g:draw_text(name, margin, y_text, 200, 3)
		y_text = y_text + line_height

		if sel == self.current_max_category_name then
			g:draw_text(name, w - (9 * margin), margin, 200, 8)
			self.current_max_category_value = 0
		end

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
