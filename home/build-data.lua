-- Pre-render script (see _quarto.yml) that writes home/gallery.json for the
-- gallery home page. Runs with Quarto's bundled Pandoc, so the site still
-- renders without R. Inputs:
--   _quarto.yml                   sidebar sections and their order
--   examples/*/index.qmd          front matter
--   home/examples.yml             curated metadata (see the comments there)
--   examples/*/model.R, module-fx.R   scanned for EpiModel calls, ERGM terms,
--                                 control.net() settings, and custom modules
--   home/epimodel-functions.tsv   EpiModel exports and their reference pages
--   renv.lock                     pinned EpiModel version

local stringify = pandoc.utils.stringify
local ptype = pandoc.utils.type

local function read_file(path)
  local f = io.open(path, "r")
  if not f then return nil end
  local s = f:read("a")
  f:close()
  return s
end

local function write_if_changed(path, text)
  if read_file(path) == text then return end
  local f = assert(io.open(path, "w"))
  f:write(text)
  f:close()
end

local function front_matter(text)
  local yaml = text:match("^%-%-%-\n(.-)\n%-%-%-")
  if not yaml then return {} end
  return pandoc.read("---\n" .. yaml .. "\n---\n", "markdown").meta
end

-- Metadata values arrive as Pandoc Inlines, Lists, or booleans.
local function str(v)
  if v == nil then return nil end
  return stringify(v)
end

local function list(v)
  local out = {}
  if v == nil then return out end
  if ptype(v) == "List" then
    for _, x in ipairs(v) do out[#out + 1] = stringify(x) end
  else
    out[1] = stringify(v)
  end
  return out
end

local function inline_html(v)
  if v == nil then return nil end
  local inl = ptype(v) == "Inlines" and v or pandoc.Inlines(stringify(v))
  local html = pandoc.write(pandoc.Pandoc({ pandoc.Plain(inl) }), "html",
    { wrap_text = "none" })
  return (html:gsub("\n+$", ""))
end

local function states(v)
  local out = {}
  if v == nil then return out end
  for _, s in ipairs(v) do
    local parts = list(s)
    out[#out + 1] = { label = parts[1], role = parts[2] or "i",
                      name = parts[3] or parts[1] }
  end
  return out
end

local function set_of(t)
  local s = {}
  for _, x in ipairs(t) do s[x] = true end
  return s
end

local function sorted_keys(t)
  local k = {}
  for x in pairs(t) do k[#k + 1] = x end
  table.sort(k)
  return k
end

-- ---- sidebar sections, in order ---------------------------------------------
local sections, section_of, sidebar_pos = {}, {}, {}
do
  local current, n = nil, 0
  for line in (read_file("_quarto.yml") or ""):gmatch("[^\n]+") do
    local sec = line:match("^%s*%-%s*section:%s*(.-)%s*$")
    if sec then
      current = sec:gsub('^"(.*)"$', "%1")
      sections[#sections + 1] = current
    end
    local ex = line:match("^%s*%-%s*examples/([^/]+)/index%.qmd%s*$")
    if ex and current then
      n = n + 1
      section_of[ex] = current
      sidebar_pos[ex] = n
    end
  end
end

-- ---- EpiModel exports ---------------------------------------------------------
local exports = {}
for line in (read_file("home/epimodel-functions.tsv") or ""):gmatch("[^\n]+") do
  if not line:match("^#") then
    local fn, topic = line:match("^([^\t]+)\t?(.*)$")
    if fn then exports[fn] = (topic ~= "" and topic or fn) end
  end
end

local epimodel_version = nil
do
  local lock = read_file("renv.lock")
  if lock then
    local ok, d = pcall(pandoc.json.decode, lock, false)
    if ok and d and d.Packages and d.Packages.EpiModel then
      epimodel_version = d.Packages.EpiModel.Version
    end
  end
end

-- ---- code scanning ------------------------------------------------------------
local ERGM_TERMS = set_of({
  "edges", "concurrent", "concurrentties", "degree", "degrange", "nodefactor",
  "nodematch", "nodemix", "nodecov", "absdiff", "absdiffcross", "isolates",
  "gwesp", "gwdegree", "triangle", "mutual", "offset", "sociality", "b1factor",
  "b2factor", "Persist", "Form", "Diss"
})

local SETTINGS = {
  "tergmLite", "resimulate.network", "cumulative.edgelist", "truncate.el.cuml",
  "save.cumulative.edgelist", "save.transmat", "save.nwstats", "nwstats.formula",
  "module.order", "attr.rules", "set.control.tergm", "set.control.ergm",
  "save.run", "epi.by", "save.other"
}

-- Drops comments, leaving strings alone.
local function strip_comments(src)
  local out = {}
  for line in (src .. "\n"):gmatch("([^\n]*)\n") do
    local q, cut, i = nil, #line + 1, 1
    while i <= #line do
      local c = line:sub(i, i)
      if q then
        if c == "\\" then i = i + 1 elseif c == q then q = nil end
      elseif c == '"' or c == "'" or c == "`" then
        q = c
      elseif c == "#" then
        cut = i
        break
      end
      i = i + 1
    end
    out[#out + 1] = line:sub(1, cut - 1)
  end
  return table.concat(out, "\n")
end

-- Right-hand sides of formulas: from each ~ to the end of the expression.
local function formula_terms(code, found)
  local pos = 1
  while true do
    local s = code:find("~", pos, true)
    if not s then break end
    local depth, i, n = 0, s + 1, #code
    while i <= n do
      local c = code:sub(i, i)
      if c == "(" then
        depth = depth + 1
      elseif c == ")" then
        if depth == 0 then break end
        depth = depth - 1
      elseif c == "," and depth == 0 then
        break
      elseif c == "\n" and depth == 0 then
        local before = code:sub(s + 1, i - 1):match("([^%s])%s*$")
        if before ~= "+" and before ~= "*" and before ~= "~" then break end
      end
      i = i + 1
    end
    for w in code:sub(s + 1, i - 1):gmatch("[%a][%w%._]*") do
      if ERGM_TERMS[w] then found[w] = true end
    end
    pos = i
  end
end

local function scan_code(dir)
  local model = strip_comments(read_file(dir .. "/model.R") or "")
  local mods = strip_comments(read_file(dir .. "/module-fx.R") or "")
  local code = model .. "\n" .. mods
  local fns, terms, settings, modules = {}, {}, {}, {}
  for name in code:gmatch("([%a%.][%w%._]*)%s*%(") do
    if exports[name] then fns[name] = true end
  end
  formula_terms(code, terms)
  for _, s in ipairs(SETTINGS) do
    local pat = "[^%w%._]" .. s:gsub("%.", "%%.") .. "%s*=[^=]"
    if ("\n" .. code):find(pat) then settings[#settings + 1] = s end
  end
  local lines = 0
  for line in (mods .. "\n"):gmatch("([^\n]*)\n") do
    local m = line:match("^([%a%.][%w%._]*)%s*<%-%s*function%s*%(%s*dat%s*,")
    if m then modules[#modules + 1] = m end
    if line:match("%S") then lines = lines + 1 end
  end
  return sorted_keys(fns), sorted_keys(terms), settings, modules, lines
end

-- ---- examples -------------------------------------------------------------------
local curated = {}
do
  local text = read_file("home/examples.yml") or ""
  local meta = pandoc.read("---\n" .. text .. "\n---\n", "markdown").meta
  if meta.examples then
    for k, v in pairs(meta.examples) do curated[k] = v end
  end
end

local LEVELS = { beginner = "Beginner", intermediate = "Intermediate", advanced = "Advanced" }

local examples, used_fns = {}, {}
local dirs = pandoc.system.list_directory("examples")
table.sort(dirs)
for _, name in ipairs(dirs) do
  local dir = "examples/" .. name
  local qmd = read_file(dir .. "/index.qmd")
  if qmd then
    local fm = front_matter(qmd)
    local c = curated[name] or {}
    if not curated[name] then
      io.stderr:write("home/build-data.lua: no entry for " .. name ..
        " in home/examples.yml; the home page shows it without filters\n")
    end
    local level, cats = nil, {}
    for _, cat in ipairs(list(fm.categories)) do
      if LEVELS[cat:lower()] then level = level or LEVELS[cat:lower()]
      else cats[#cats + 1] = cat end
    end
    local fns, terms, settings, modules, lines = scan_code(dir)
    for _, f in ipairs(fns) do used_fns[f] = exports[f] end
    examples[#examples + 1] = {
      id = name,
      href = dir .. "/index.html",
      title = str(fm.title) or name,
      short = str(c.short) or str(fm.title) or name,
      summary = inline_html(c.summary or fm.description) or "",
      description = str(fm.description) or "",
      author = str(fm.author) or "",
      date = str(fm.date) or "",
      order = tonumber(str(fm.order)) or 999,
      section = section_of[name] or "More examples",
      sidebar = sidebar_pos[name] or 999,
      level = level or "Intermediate",
      categories = cats,
      states = states(c.states),
      side = states(c.side),
      loop = c.loop == true,
      model = list(c.model),
      network = list(c.network),
      disease = str(c.disease) or "Other",
      interventions = list(c.interventions),
      techniques = list(c.techniques),
      builds_on = list(c["builds-on"]),
      functions = fns,
      terms = terms,
      settings = settings,
      modules = modules,
      module_lines = lines
    }
  end
end

table.sort(examples, function(a, b)
  if a.order ~= b.order then return a.order < b.order end
  return a.sidebar < b.sidebar
end)

local function has(t, x) for _, y in ipairs(t) do if y == x then return true end end end
for _, e in ipairs(examples) do
  if not has(sections, e.section) then sections[#sections + 1] = e.section end
end

-- Lua tables cannot tell an empty array from an empty object; give empty
-- lists a sentinel the JSON encoder writes as [].
local function arrays(t)
  for k, v in pairs(t) do
    if type(v) == "table" then
      if next(v) == nil then t[k] = pandoc.List({}) else arrays(v) end
    end
  end
  return t
end

local data = {
  epimodel = epimodel_version,
  sections = sections,
  reference = "https://epimodel.github.io/EpiModel/reference/",
  topics = used_fns,
  examples = examples
}
write_if_changed("home/gallery.json", pandoc.json.encode(arrays(data)) .. "\n")
