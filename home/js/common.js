// Shared helpers for the gallery home page. Colors are CSS variables defined in
// home/home.css; text from gallery.json is inserted with textContent except the
// card summaries, which build-data.lua renders from the repository's own Markdown.

// Compartment roles (see home/examples.yml). The same colors fill the network
// glyphs, the compartment chips, and the hero animation. Infectious nodes are
// also drawn larger with a halo, so I and R never rely on red and green alone.
export const ROLE = {
  s: { color: "var(--st-s)", hex: "#2a78d6", name: "Susceptible" },
  e: { color: "var(--st-e)", hex: "#eda100", name: "Latent" },
  i: { color: "var(--st-i)", hex: "#e34948", name: "Infectious" },
  r: { color: "var(--st-r)", hex: "#1baf7a", name: "Recovered or immune" },
  v: { color: "var(--st-v)", hex: "#6250d6", name: "Vaccine-immune" },
  t: { color: "var(--st-t)", hex: "#d55181", name: "On treatment" },
  x: { color: "var(--st-x)", hex: "#9aa3ab", name: "Removed" }
};
export const roleOf = (r) => ROLE[r] || ROLE.i;

export function el(tag, attrs = {}, ...kids) {
  const n = document.createElement(tag);
  setAttrs(n, attrs);
  for (const c of kids.flat()) if (c != null && c !== false) n.append(c.nodeType ? c : document.createTextNode(String(c)));
  return n;
}

const SVGNS = "http://www.w3.org/2000/svg";
export function svg(tag, attrs = {}, ...kids) {
  const n = document.createElementNS(SVGNS, tag);
  setAttrs(n, attrs);
  for (const c of kids.flat()) if (c != null && c !== false) n.append(c.nodeType ? c : document.createTextNode(String(c)));
  return n;
}

function setAttrs(n, attrs) {
  for (const [k, v] of Object.entries(attrs)) {
    if (v == null || v === false) continue;
    if (k === "class") n.setAttribute("class", v);
    else if (k === "text") n.textContent = v;
    else if (k === "html") n.innerHTML = v;
    else if (k === "style" && typeof v === "object") Object.assign(n.style, v);
    else if (k.startsWith("on") && typeof v === "function") n.addEventListener(k.slice(2), v);
    else n.setAttribute(k, v === true ? "" : v);
  }
}

// Seeded PRNG (mulberry32) so every card draws the same network on every visit.
export function rng(seed) {
  let a = typeof seed === "string" ? hash(seed) : seed >>> 0;
  return () => {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
export function hash(s) {
  let h = 2166136261;
  for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
  return h >>> 0;
}

export const reducedMotion = () =>
  typeof window !== "undefined" && window.matchMedia &&
  window.matchMedia("(prefers-reduced-motion: reduce)").matches;

const MONTHS = ["January", "February", "March", "April", "May", "June", "July",
  "August", "September", "October", "November", "December"];
export function monthYear(iso) {
  const m = /^(\d{4})-(\d{2})/.exec(iso || "");
  return m ? `${MONTHS[+m[2] - 1]} ${m[1]}` : "";
}
export function daysSince(iso) {
  const t = Date.parse(iso);
  return isNaN(t) ? Infinity : (Date.now() - t) / 864e5;
}

// The compartment strip on a card: S -> E -> I -> R, with a loop mark when the
// last state returns to susceptible and side states after a separator.
export function stateStrip(ex, { compact = false } = {}) {
  const wrap = el("div", { class: "strip" + (compact ? " compact" : ""), role: "img",
    "aria-label": "Compartments: " + ex.states.map((s) => s.name).join(", then ") +
      (ex.loop ? ", then back to susceptible" : "") +
      (ex.side.length ? "; also " + ex.side.map((s) => s.name).join(", ") : "") });
  ex.states.forEach((s, i) => {
    if (i > 0) wrap.append(el("span", { class: "arr", "aria-hidden": "true", text: "→" }));
    wrap.append(chip(s));
  });
  if (ex.loop) wrap.append(el("span", { class: "arr loop", "aria-hidden": "true", title: "Returns to susceptible", text: "↺" }));
  if (ex.side.length) {
    wrap.append(el("span", { class: "sep", "aria-hidden": "true" }));
    for (const s of ex.side) wrap.append(chip(s));
  }
  return wrap;
}
function chip(s) {
  const c = el("span", { class: "st", title: s.name, "aria-hidden": "true", text: s.label });
  c.style.setProperty("--c", roleOf(s.role).color);
  return c;
}

// Width of a container, re-reported on resize.
export function onWidth(node, fn) {
  let last = -1;
  const ro = new ResizeObserver((entries) => {
    const w = Math.round(entries[0].contentRect.width);
    // Deferred a frame so a resize inside fn does not re-enter the observer.
    if (w > 0 && Math.abs(w - last) > 2) { last = w; requestAnimationFrame(() => fn(w)); }
  });
  ro.observe(node);
  return () => ro.disconnect();
}
