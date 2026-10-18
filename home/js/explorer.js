// "Find an example": search, facet filters, and three views of the same filtered
// set (cards, the learning path, and an index of the EpiModel code each example
// uses). Filter state lives in the URL, so a filtered view can be linked.

import { el, svg, stateStrip, monthYear, daysSince, onWidth } from "./common.js";
import { glyph, networkLabel } from "./glyph.js";

const FACETS = [
  { key: "level", label: "Level", get: (e) => [e.level], order: ["Beginner", "Intermediate", "Advanced"] },
  { key: "network", label: "Network", get: (e) => e.network, order: ["One layer", "Multiple layers", "Households", "Observed network"] },
  { key: "model", label: "Compartments", get: (e) => e.model, order: ["SI", "SIS", "SIR", "SIRS", "SEI", "SEIR", "SEIRS", "Multi-stage"] },
  { key: "disease", label: "Disease area", get: (e) => [e.disease] },
  { key: "intervention", label: "Intervention", get: (e) => e.interventions },
  { key: "technique", label: "Technique", get: (e) => e.techniques }
];
const VIEWS = [
  { value: "cards", label: "Cards" },
  { value: "path", label: "Learning path" },
  { value: "code", label: "Code index" }
];
const SORTS = [
  { value: "path", label: "Suggested order" },
  { value: "new", label: "Newest first" },
  { value: "az", label: "A to Z" }
];
const LEVEL_DOTS = { Beginner: 1, Intermediate: 2, Advanced: 3 };

// ---- state and matching -----------------------------------------------------------------
function readURL() {
  const p = new URLSearchParams(location.search);
  // Quarto's site search reads and removes ?q=, so the explorer uses ?search=.
  const s = { q: p.get("search") || "", sort: p.get("sort") || "path", view: p.get("view") || "cards", f: {} };
  for (const f of FACETS) s.f[f.key] = new Set((p.get(f.key) || "").split("|").filter(Boolean));
  if (!VIEWS.some((v) => v.value === s.view)) s.view = "cards";
  if (!SORTS.some((v) => v.value === s.sort)) s.sort = "path";
  return s;
}
function writeURL(s) {
  const p = new URLSearchParams();
  if (s.q) p.set("search", s.q);
  for (const f of FACETS) if (s.f[f.key].size) p.set(f.key, [...s.f[f.key]].join("|"));
  if (s.sort !== "path") p.set("sort", s.sort);
  if (s.view !== "cards") p.set("view", s.view);
  const qs = p.toString();
  history.replaceState(null, "", location.pathname + (qs ? "?" + qs : "") + location.hash);
}

const tokens = (q) => q.toLowerCase().replace(/[()]/g, " ").split(/\s+/).filter(Boolean);

function haystack(e) {
  return [e.title, e.short, e.summary.replace(/<[^>]+>/g, ""), e.description, e.author, e.section, e.level,
    e.disease, ...e.categories, ...e.model, ...e.network, ...e.interventions, ...e.techniques,
    ...e.states.map((s) => s.name), ...e.functions, ...e.terms, ...e.settings, ...e.modules].join(" ").toLowerCase();
}

// Names in the example's code that a search term picks out, for the card's "uses" line.
function codeHits(e, toks) {
  if (!toks.length) return [];
  const names = [...e.functions.map((n) => n + "()"), ...e.terms.map((n) => n + " term"), ...e.settings.map((n) => n + " ="), ...e.modules.map((n) => n + "()")];
  return names.filter((n) => toks.some((t) => t.length >= 3 && n.toLowerCase().includes(t)));
}

function passes(e, s, skip) {
  for (const f of FACETS) {
    if (f.key === skip || !s.f[f.key].size) continue;
    if (!f.get(e).some((v) => s.f[f.key].has(v))) return false;
  }
  const toks = tokens(s.q);
  return toks.every((t) => e._hay.includes(t));
}

function sorted(list, sort) {
  const a = list.slice();
  if (sort === "new") a.sort((x, y) => (y.date || "").localeCompare(x.date || "") || x.order - y.order);
  else if (sort === "az") a.sort((x, y) => x.title.localeCompare(y.title));
  else a.sort((x, y) => x.order - y.order);
  return a;
}

function levelDots(level) {
  const n = LEVEL_DOTS[level] || 2;
  return el("span", { class: "dots", "aria-hidden": "true" }, [1, 2, 3].map((i) => el("i", { class: i <= n ? "on" : "" })));
}

// ---- cards --------------------------------------------------------------------------------
function card(e) {
  const g = glyph(e);
  const isNew = daysSince(e.date) < 90;
  const art = el("div", { class: "card-art" }, g.node,
    el("span", { class: "lvl" }, levelDots(e.level), e.level),
    isNew ? el("span", { class: "new", text: "New" }) : null,
    el("span", { class: "net", text: networkLabel(e) }));
  const tags = [...e.interventions, ...e.techniques];
  const shown = tags.slice(0, 4);
  const match = el("div", { class: "card-match", hidden: true });
  const a = el("a", { class: "card-ex", href: e.href },
    art,
    el("div", { class: "card-body" },
      el("div", { class: "card-kicker", text: e.section }),
      el("h3", { text: e.title }),
      stateStrip(e),
      el("p", { class: "card-sum", html: e.summary }),
      match,
      el("div", { class: "card-tags" }, shown.map((t) => el("span", { class: "tag", text: t })),
        tags.length > shown.length ? el("span", { class: "tag more", text: `+${tags.length - shown.length}`, title: tags.slice(4).join(", ") }) : null),
      el("div", { class: "card-foot" },
        el("span", { title: e.modules.length ? "Custom modules: " + e.modules.join(", ") : "", text: `${e.modules.length} custom module${e.modules.length === 1 ? "" : "s"} · ${e.module_lines} lines` }),
        el("span", { text: [e.author, monthYear(e.date)].filter(Boolean).join(" · ") }))));
  a.addEventListener("mouseenter", g.play);
  a.addEventListener("focus", g.play);
  a.addEventListener("mouseleave", g.stop);
  a.addEventListener("blur", g.stop);
  a._match = match;
  return a;
}

function setMatch(c, hits) {
  const m = c._match;
  m.replaceChildren();
  m.hidden = !hits.length;
  if (hits.length) m.append(el("span", { text: "Uses " }), ...hits.slice(0, 4).map((h) => el("code", { text: h })), hits.length > 4 ? ` and ${hits.length - 4} more` : "");
}

// ---- learning path ------------------------------------------------------------------------
function pathView(G, byId) {
  const box = el("div", { class: "path" });
  const detail = el("div", { class: "path-detail", "aria-live": "polite" });
  const legend = el("div", { class: "path-legend" },
    el("span", {}, el("i", { class: "sw up" }), "Read first"),
    el("span", {}, el("i", { class: "sw down" }), "Builds on it"),
    el("span", {}, el("i", { class: "sw dim" }), "Outside the current filters"));
  const scroller = el("div", { class: "path-scroll" });
  box.append(legend, scroller, detail);

  // Depth = longest chain of prerequisites; each depth is a column.
  const ids = G.examples.map((e) => e.id);
  const pre = new Map(ids.map((id) => [id, byId.get(id).builds_on.filter((p) => byId.has(p))]));
  const post = new Map(ids.map((id) => [id, []]));
  for (const [id, ps] of pre) for (const p of ps) post.get(p).push(id);
  const depth = new Map();
  const d = (id, seen = new Set()) => {
    if (depth.has(id)) return depth.get(id);
    if (seen.has(id)) return 0;
    seen.add(id);
    const v = pre.get(id).length ? 1 + Math.max(...pre.get(id).map((p) => d(p, seen))) : 0;
    depth.set(id, v);
    return v;
  };
  ids.forEach((id) => d(id));
  const nCol = Math.max(...depth.values()) + 1;
  const cols = Array.from({ length: nCol }, () => []);
  for (const e of G.examples) cols[depth.get(e.id)].push(e.id);
  // Order each column by the mean row of its prerequisites, two sweeps.
  const rowOf = new Map();
  cols[0].forEach((id, i) => rowOf.set(id, i));
  for (let sweep = 0; sweep < 2; sweep++) {
    for (let c = 1; c < nCol; c++) {
      const key = (id) => { const ps = pre.get(id); return ps.length ? ps.reduce((s, p) => s + (rowOf.get(p) ?? 0), 0) / ps.length : 0; };
      cols[c].sort((a, b) => key(a) - key(b) || byId.get(a).order - byId.get(b).order);
      cols[c].forEach((id, i) => rowOf.set(id, i));
    }
  }
  const anc = (id, acc = new Set()) => { for (const p of pre.get(id)) if (!acc.has(p)) { acc.add(p); anc(p, acc); } return acc; };
  const desc = (id, acc = new Set()) => { for (const p of post.get(id)) if (!acc.has(p)) { acc.add(p); desc(p, acc); } return acc; };

  let nodes = new Map(), links = [], active = new Set(ids), focus = null;
  const COLHEAD = ["Start here", "Next", "Then", "Further", "Further still"];

  function draw(width) {
    const W = Math.max(760, width), colW = W / nCol, boxW = Math.min(190, colW - 34), boxH = 44, rowH = 58, top = 34;
    const maxRows = Math.max(...cols.map((c) => c.length));
    const H = top + maxRows * rowH + 6;
    const root = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, class: "path-svg", role: "group", "aria-label": "Learning path: examples arranged by what to read first" });
    const defs = svg("defs", {});
    for (const k of ["base", "up", "down"]) {
      defs.append(svg("marker", { id: "arr-" + k, viewBox: "0 0 8 8", refX: 7, refY: 4, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" },
        svg("path", { d: "M0,0 L8,4 L0,8 z", class: "arrow " + k })));
    }
    root.append(defs);
    const gL = svg("g", { class: "links" }), gN = svg("g", { class: "nodes" });
    root.append(gL, gN);
    const pos = new Map();
    cols.forEach((c, ci) => {
      const x = ci * colW + (colW - boxW) / 2;
      const off = ((maxRows - c.length) * rowH) / 2;
      c.forEach((id, ri) => pos.set(id, { x, y: top + off + ri * rowH }));
      root.append(svg("text", { x: ci * colW + colW / 2, y: 16, class: "col-head", "text-anchor": "middle" }, `${ci + 1}. ${COLHEAD[ci] || "Further"}`));
    });
    links = [];
    for (const [id, ps] of pre) for (const p of ps) {
      const a = pos.get(p), b = pos.get(id);
      const x1 = a.x + boxW, y1 = a.y + boxH / 2, x2 = b.x - 2, y2 = b.y + boxH / 2, mx = (x1 + x2) / 2;
      const path = svg("path", { d: `M${x1},${y1} C${mx},${y1} ${mx},${y2} ${x2},${y2}`, class: "link", "marker-end": "url(#arr-base)" });
      gL.append(path);
      links.push({ from: p, to: id, path });
    }
    nodes = new Map();
    for (const e of G.examples) {
      const { x, y } = pos.get(e.id);
      const g = svg("a", { href: e.href, class: "pnode", "aria-label": `${e.title}. ${e.level}.` });
      const label = svg("text", { x: x + 11, class: "pn-title" });
      g.append(svg("rect", { x, y, width: boxW, height: boxH, rx: 8 }), label);
      gN.append(g);
      fitLabel(label, e.short, boxW - 20, y + boxH / 2);
      g.addEventListener("mouseenter", () => setFocus(e.id));
      g.addEventListener("focus", () => setFocus(e.id));
      g.addEventListener("mouseleave", () => setFocus(null));
      g.addEventListener("blur", () => setFocus(null));
      nodes.set(e.id, g);
    }
    scroller.replaceChildren(root);
    for (const f of pending) f();
    pending = [];
    paint();
  }

  // Short names wrap to a second line when they do not fit the box. Text can
  // only be measured once the SVG is in the document, so this runs after draw().
  let pending = [];
  function fitLabel(t, text, maxW, cy) {
    pending.push(() => {
      t.textContent = text;
      if (t.getComputedTextLength() <= maxW) { t.setAttribute("y", cy + 4.5); return; }
      const words = text.split(" ");
      let a = words.slice(0, -1), b = words.slice(-1);
      t.textContent = a.join(" ");
      while (a.length > 1 && t.getComputedTextLength() > maxW) { b.unshift(a.pop()); t.textContent = a.join(" "); }
      t.textContent = "";
      t.append(svg("tspan", { x: t.getAttribute("x"), y: cy - 3 }, a.join(" ")), svg("tspan", { x: t.getAttribute("x"), y: cy + 12 }, b.join(" ")));
    });
  }

  function setFocus(id) { focus = id; paint(); }

  function paint() {
    const up = focus ? anc(focus) : new Set(), down = focus ? desc(focus) : new Set();
    for (const [id, g] of nodes) {
      g.classList.toggle("dim", !active.has(id));
      g.classList.toggle("focus", id === focus);
      g.classList.toggle("up", up.has(id));
      g.classList.toggle("down", down.has(id));
      g.classList.toggle("fade", !!focus && id !== focus && !up.has(id) && !down.has(id));
    }
    for (const l of links) {
      const isUp = focus && (l.to === focus || up.has(l.to)) && (up.has(l.from));
      const isDown = focus && (l.from === focus || down.has(l.from)) && down.has(l.to);
      l.path.setAttribute("class", "link" + (isUp ? " up" : isDown ? " down" : focus ? " fade" : ""));
      l.path.setAttribute("marker-end", `url(#arr-${isUp ? "up" : isDown ? "down" : "base"})`);
    }
    detail.replaceChildren();
    if (!focus) {
      detail.append(el("p", { class: "muted", text: "Hover over or tab to an example to see what to read before it and what builds on it. Each arrow points from an example to one that reuses its ideas." }));
      return;
    }
    const e = byId.get(focus);
    const names = (set) => [...set].sort((a, b) => depth.get(a) - depth.get(b) || byId.get(a).order - byId.get(b).order).map((id) => byId.get(id).short);
    detail.append(el("div", { class: "pd-title" }, e.title, el("span", { class: "pd-lvl" }, levelDots(e.level), `${e.level} · ${e.section}`)),
      el("p", { class: "pd-sum", html: e.summary }),
      el("div", { class: "pd-rel" },
        el("span", { class: "k up" }, "Read first: "), up.size ? names(up).join(", ") : "nothing; this is a starting point",
        el("br"),
        el("span", { class: "k down" }, "Builds on it: "), down.size ? names(down).join(", ") : "no other examples yet"));
  }

  onWidth(scroller, draw);
  return {
    node: box,
    update(list) { active = new Set(list.map((e) => e.id)); paint(); }
  };
}

// ---- code index -----------------------------------------------------------------------------
function codeView(G) {
  const box = el("div", { class: "code-index" });
  const ref = (fn) => G.reference + (G.topics[fn] || fn) + ".html";
  const GROUPS = [
    { key: "functions", label: "EpiModel functions", fmt: (n) => n + "()", link: ref },
    { key: "terms", label: "ERGM terms in formulas", fmt: (n) => n },
    { key: "settings", label: "Simulation settings", fmt: (n) => n }
  ];

  function update(list, toks) {
    box.replaceChildren();
    if (!list.length) return;
    box.append(el("p", { class: "muted", text: "Read from each example's model.R and module-fx.R when the site is built. A dot marks an example that calls the function, uses the ERGM term in a formation, dissolution, or diagnostic formula, or sets the option in control.net(), netest(), or netdx(). Function names link to the EpiModel reference." }));
    const table = el("table", { class: "cx" });
    const head = el("tr", {}, el("th", { class: "corner", scope: "col" }, el("span", { class: "sr-only", text: "Name" })),
      list.map((e, j) => el("th", { scope: "col", class: "colh", "data-col": j }, el("a", { href: e.href, title: e.title }, el("span", { text: e.short })))),
      el("th", { scope: "col", class: "cnt", text: "n" }));
    table.append(el("thead", {}, head));
    const body = el("tbody");
    const hit = (name) => toks.some((t) => t.length >= 3 && name.toLowerCase().includes(t));
    for (const g of GROUPS) {
      const counts = new Map();
      for (const e of list) for (const n of e[g.key]) counts.set(n, (counts.get(n) || 0) + 1);
      const everyone = [...counts].filter(([, c]) => c === list.length && list.length > 2).map(([n]) => n).sort();
      const rows = [...counts].filter(([n]) => !everyone.includes(n)).sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
      body.append(el("tr", { class: "grp" }, el("th", { colspan: list.length + 2, scope: "colgroup" }, g.label,
        everyone.length ? el("span", { class: "every" }, ` · in every example shown: `, ...everyone.map((n, i) => [i ? " " : "", g.link ? el("a", { href: g.link(n), target: "_blank", rel: "noopener" }, el("code", { text: g.fmt(n) })) : el("code", { text: g.fmt(n) })])) : null)));
      for (const [n, c] of rows) {
        const label = g.link ? el("a", { href: g.link(n), target: "_blank", rel: "noopener" }, el("code", { text: g.fmt(n) })) : el("code", { text: g.fmt(n) });
        body.append(el("tr", { class: hit(n) ? "hit" : "" }, el("th", { scope: "row" }, label),
          list.map((e, j) => el("td", { "data-col": j, title: e[g.key].includes(n) ? `${g.fmt(n)} in ${e.short}` : "" }, e[g.key].includes(n) ? el("i", { class: "dot" }) : null)),
          el("td", { class: "cnt", text: String(c) })));
      }
    }
    table.append(body);
    table.addEventListener("mouseover", (ev) => {
      const c = ev.target.closest("[data-col]");
      table.querySelectorAll(".colhi").forEach((n) => n.classList.remove("colhi"));
      if (c) table.querySelectorAll(`[data-col="${c.dataset.col}"]`).forEach((n) => n.classList.add("colhi"));
    });
    table.addEventListener("mouseleave", () => table.querySelectorAll(".colhi").forEach((n) => n.classList.remove("colhi")));
    box.append(el("div", { class: "cx-scroll" }, table));
  }
  return { node: box, update };
}

// ---- the explorer -----------------------------------------------------------------------------
export function explorer(G) {
  const S = readURL();
  const byId = new Map(G.examples.map((e) => [e.id, e]));
  for (const e of G.examples) e._hay = haystack(e);

  // Facet values, in a fixed order where one exists, otherwise most common first.
  const values = {};
  for (const f of FACETS) {
    const c = new Map();
    for (const e of G.examples) for (const v of f.get(e)) c.set(v, (c.get(v) || 0) + 1);
    const ord = f.order || [];
    values[f.key] = [...c.keys()].sort((a, b) => {
      const ia = ord.indexOf(a), ib = ord.indexOf(b);
      if (ia >= 0 || ib >= 0) return (ia < 0 ? 99 : ia) - (ib < 0 ? 99 : ib);
      return c.get(b) - c.get(a) || a.localeCompare(b);
    });
  }

  // Search and toolbar.
  const fnNames = [...new Set(G.examples.flatMap((e) => e.functions))].sort();
  const list = el("datalist", { id: "gx-suggest" }, fnNames.map((n) => el("option", { value: n })));
  const search = el("input", { type: "search", id: "gx-q", placeholder: "Search examples, diseases, or functions like get_partners", value: S.q, list: "gx-suggest", autocomplete: "off", "aria-label": "Search examples" });
  const sortSel = el("select", { "aria-label": "Sort cards" }, SORTS.map((o) => el("option", { value: o.value, text: o.label })));
  sortSel.value = S.sort;
  const viewSeg = el("div", { class: "seg", role: "group", "aria-label": "View" });
  const viewBtns = VIEWS.map((v) => {
    const b = el("button", { type: "button", text: v.label, "aria-pressed": String(v.value === S.view) });
    b.addEventListener("click", () => { S.view = v.value; render(); });
    viewSeg.append(b);
    return [v.value, b];
  });
  const sortWrap = el("label", { class: "sort" }, el("span", { text: "Sort" }), sortSel);

  // Facet rail.
  const chips = new Map();
  const rail = el("div", { class: "rail-inner" });
  for (const f of FACETS) {
    const group = el("fieldset", { class: "facet" }, el("legend", { text: f.label }));
    const wrap = el("div", { class: "chips" });
    for (const v of values[f.key]) {
      const n = el("span", { class: "n" });
      const b = el("button", { type: "button", class: "chip", "aria-pressed": "false" }, el("span", { text: v }), n);
      b.addEventListener("click", () => {
        const set = S.f[f.key];
        set.has(v) ? set.delete(v) : set.add(v);
        render();
      });
      chips.set(f.key + "\u0000" + v, { b, n, f, v });
      wrap.append(b);
    }
    group.append(wrap);
    rail.append(group);
  }
  const railSummary = el("summary", { text: "Filters" });
  const railBox = el("details", { class: "rail" }, railSummary, rail);
  const wide = window.matchMedia("(min-width: 960px)");
  const syncRail = () => { railBox.open = wide.matches || railBox.dataset.userOpen === "1"; };
  railSummary.addEventListener("click", () => { railBox.dataset.userOpen = railBox.open ? "0" : "1"; });
  wide.addEventListener("change", syncRail);
  syncRail();

  // Results.
  const count = el("div", { class: "count", "aria-live": "polite" });
  const active = el("div", { class: "active-filters" });
  const cards = new Map(G.examples.map((e) => [e.id, card(e)]));
  const grid = el("div", { class: "cards" });
  const empty = el("div", { class: "empty", hidden: true },
    el("p", { text: "No example matches all of these filters." }),
    el("p", {}, "Try removing one, or ", el("a", { href: "https://github.com/EpiModel/EpiModel-Gallery/issues", text: "ask for an example" }), " that covers what you need."));
  const path = pathView(G, byId);
  const code = codeView(G);
  const panes = { cards: el("div", { class: "pane" }, grid), path: el("div", { class: "pane" }, path.node), code: el("div", { class: "pane" }, code.node) };

  const toolbar = el("div", { class: "toolbar" },
    el("div", { class: "searchbox" }, search, list),
    el("div", { class: "tools" }, sortWrap, viewSeg));
  const results = el("div", { class: "results" }, el("div", { class: "result-head" }, count, active), empty, panes.cards, panes.path, panes.code);
  const root = el("section", { class: "explorer", id: "explore", "aria-label": "Find an example" },
    el("h2", { class: "ex-h", text: "Find an example" }),
    toolbar,
    el("div", { class: "ex-grid" }, railBox, results));

  let timer = null;
  search.addEventListener("input", () => { clearTimeout(timer); timer = setTimeout(() => { S.q = search.value.trim(); render(); }, 120); });
  sortSel.addEventListener("change", () => { S.sort = sortSel.value; render(); });
  document.addEventListener("keydown", (ev) => {
    if (ev.key === "/" && !/^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement?.tagName || "")) { ev.preventDefault(); search.focus(); }
  });

  function render() {
    const toks = tokens(S.q);
    const shown = sorted(G.examples.filter((e) => passes(e, S)), S.sort);
    const n = shown.length, total = G.examples.length;
    count.textContent = n === total && !toks.length ? `${total} examples` : `${n} of ${total} examples`;

    // Chip counts: what each value would give, holding the other filters fixed.
    for (const { b, n: cnt, f, v } of chips.values()) {
      const k = G.examples.filter((e) => passes(e, S, f.key) && f.get(e).includes(v)).length;
      const on = S.f[f.key].has(v);
      cnt.textContent = String(k);
      b.setAttribute("aria-pressed", String(on));
      b.disabled = !on && k === 0;
    }
    const nActive = FACETS.reduce((s, f) => s + S.f[f.key].size, 0);
    railSummary.textContent = nActive ? `Filters (${nActive} on)` : "Filters";

    active.replaceChildren();
    for (const f of FACETS) for (const v of S.f[f.key]) {
      const x = el("button", { type: "button", class: "pill", "aria-label": `Remove filter ${v}` }, el("span", { text: v }), el("span", { "aria-hidden": "true", text: "×" }));
      x.addEventListener("click", () => { S.f[f.key].delete(v); render(); });
      active.append(x);
    }
    if (nActive || S.q) {
      const clr = el("button", { type: "button", class: "clear", text: "Clear all" });
      clr.addEventListener("click", () => { for (const f of FACETS) S.f[f.key].clear(); S.q = ""; search.value = ""; render(); });
      active.append(clr);
    }

    for (const [v, b] of viewBtns) b.setAttribute("aria-pressed", String(v === S.view));
    for (const [k, p] of Object.entries(panes)) p.hidden = k !== S.view;
    sortWrap.hidden = S.view !== "cards";
    empty.hidden = n > 0;

    if (S.view === "cards") {
      grid.replaceChildren(...shown.map((e) => { const c = cards.get(e.id); setMatch(c, codeHits(e, toks)); return c; }));
    } else if (S.view === "path") {
      path.update(shown);
    } else {
      code.update(sorted(shown, "path"), toks);
    }
    writeURL(S);
  }

  render();
  root.showPath = () => { S.view = "path"; render(); };
  return root;
}
