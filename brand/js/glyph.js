// Network art: a small seeded network in a given structure (one layer, two
// layers, households plus community ties, or an observed network), colored by a
// toy epidemic that runs through a given list of compartments. Gallery cards use
// it for each example; www.epimodel.org uses it for its resource cards. The art is
// a fixed snapshot; play() replays the outbreak.

import { svg, rng, roleOf, reducedMotion } from "./brand.js";

const W = 320, H = 120;

export function networkType(ex) {
  if (ex.network.includes("Observed network")) return "observed";
  if (ex.network.includes("Households")) return "households";
  if (ex.network.includes("Multiple layers")) return "multilayer";
  return "single";
}

export function networkLabel(ex) {
  const t = networkType(ex);
  if (t === "households") return "Households + community";
  if (t === "observed") return "Observed contacts";
  if (t === "multilayer") return "Two partnership layers";
  return "One contact layer";
}

const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);

function jitteredGrid(R, cols, rows, jitter = 0.32) {
  const padX = 14, padY = 13;
  const cw = (W - 2 * padX) / (cols - 0.5), ch = (H - 2 * padY) / (rows - 1);
  const pts = [];
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      const x = padX + c * cw + (r % 2 ? cw / 2 : 0) + (R() - 0.5) * 2 * jitter * cw;
      const y = padY + r * ch + (R() - 0.5) * 2 * jitter * ch;
      pts.push({ x: clamp(x, 8, W - 8), y: clamp(y, 8, H - 8) });
    }
  }
  return pts;
}

function addEdge(edges, seen, a, b, kind) {
  if (a === b) return false;
  const k = a < b ? a * 1000 + b : b * 1000 + a;
  if (seen.has(k)) return false;
  seen.add(k);
  edges.push({ a, b, kind });
  return true;
}

function nearest(nodes, i, cands) {
  return cands.filter((j) => j !== i).sort((p, q) => dist(nodes[i], nodes[p]) - dist(nodes[i], nodes[q]));
}

export function buildNetwork(ex) {
  const R = rng("net:" + ex.id);
  const type = networkType(ex);
  const edges = [], seen = new Set();
  let nodes;

  if (type === "households") {
    nodes = [];
    const centers = jitteredGrid(R, 6, 2, 0.18);
    centers.forEach((c, h) => {
      const size = 2 + Math.floor(R() * 4);
      const rad = 6 + size * 2, rot = R() * Math.PI * 2;
      const first = nodes.length;
      for (let m = 0; m < size; m++) {
        const a = rot + (m / size) * Math.PI * 2;
        nodes.push({ x: clamp(c.x + rad * Math.cos(a), 6, W - 6), y: clamp(c.y + rad * Math.sin(a), 6, H - 6), hh: h });
      }
      for (let i = first; i < nodes.length; i++)
        for (let j = i + 1; j < nodes.length; j++) addEdge(edges, seen, i, j, "house");
    });
    const all = nodes.map((_, i) => i);
    let tries = 0;
    while (edges.filter((e) => e.kind === "comm").length < Math.round(nodes.length * 0.6) && tries++ < 500) {
      const i = Math.floor(R() * nodes.length);
      const cands = nearest(nodes, i, all).filter((j) => nodes[j].hh !== nodes[i].hh).slice(0, 6);
      if (cands.length) addEdge(edges, seen, i, cands[Math.floor(R() * cands.length)], "comm");
    }
  } else if (type === "observed") {
    // 24 animals in three overlapping home-range clusters; most contacts are local.
    const centers = [{ x: 62, y: 58 }, { x: 160, y: 66 }, { x: 258, y: 54 }];
    nodes = [];
    for (let i = 0; i < 24; i++) {
      const c = centers[i % 3];
      const a = R() * Math.PI * 2, r = 6 + Math.sqrt(R()) * 40;
      nodes.push({ x: clamp(c.x + r * Math.cos(a) * 1.15, 8, W - 8), y: clamp(c.y + r * Math.sin(a) * 0.95, 8, H - 8) });
    }
    for (let i = 0; i < nodes.length; i++)
      for (let j = i + 1; j < nodes.length; j++)
        if (dist(nodes[i], nodes[j]) < 80 && R() < 0.9 * Math.exp(-dist(nodes[i], nodes[j]) / 28)) addEdge(edges, seen, i, j, "obs");
  } else {
    nodes = jitteredGrid(R, 10, 4);
    const all = nodes.map((_, i) => i);
    if (type === "multilayer") {
      const order = all.slice().sort(() => R() - 0.5), used = new Set();
      for (const i of order) {
        if (used.has(i) || R() < 0.3) continue;
        const j = nearest(nodes, i, all).find((k) => !used.has(k));
        if (j != null && dist(nodes[i], nodes[j]) < 60) { used.add(i); used.add(j); addEdge(edges, seen, i, j, "main"); }
      }
      let tries = 0;
      while (edges.filter((e) => e.kind === "casual").length < Math.round(nodes.length * 0.5) && tries++ < 500) {
        const i = Math.floor(R() * nodes.length);
        const c = nearest(nodes, i, all).slice(1, 6);
        addEdge(edges, seen, i, c[Math.floor(R() * c.length)], "casual");
      }
    } else {
      for (let i = 0; i < nodes.length; i++) {
        const k = [0, 1, 1, 1, 2, 2][Math.floor(R() * 6)];
        const c = nearest(nodes, i, all).slice(0, 4);
        for (let m = 0; m < k; m++) {
          const j = c[Math.floor(R() * c.length)];
          if (dist(nodes[i], nodes[j]) < 75) addEdge(edges, seen, i, j, "tie");
        }
      }
      for (let m = 0; m < 2; m++) addEdge(edges, seen, Math.floor(R() * nodes.length), Math.floor(R() * nodes.length), "tie");
    }
  }
  return { type, nodes, edges };
}

// ---- the toy epidemic --------------------------------------------------------------
// State is an index into ex.states; side states are stored as -1, -2, ...
const ADVANCE = { s: 0, e: 0.3, i: 0.17, r: 0.05, v: 0, t: 0.3, x: 0 };
const TRANSMIT = { house: 0.5, comm: 0.22, main: 0.35, casual: 0.22, obs: 0.18, tie: 0.32 };

function makeSim(ex, net) {
  const chain = ex.states.length ? ex.states : [{ role: "s" }, { role: "i" }];
  const firstInf = Math.max(1, chain.findIndex((s) => s.role === "i"));
  const n = net.nodes.length;
  const lastRole = chain[chain.length - 1].role;

  function init(seed) {
    const R = rng(seed);
    const st = new Int8Array(n);
    if (ex.side.length) for (let i = 0; i < n; i++) if (R() < 0.16) st[i] = -1;
    const seeds = net.type === "single" || net.type === "multilayer" ? 2 : 1;
    for (let s = 0; s < seeds; s++) {
      let i, guard = 0;
      do { i = Math.floor(R() * n); } while (st[i] !== 0 && guard++ < 50);
      st[i] = firstInf;
    }
    return { st, R, t: 0 };
  }

  function step(S) {
    const { st, R } = S, next = st.slice();
    for (const e of net.edges) {
      const ra = st[e.a] > 0 ? chain[st[e.a]].role : null;
      const rb = st[e.b] > 0 ? chain[st[e.b]].role : null;
      const p = TRANSMIT[e.kind] ?? 0.25;
      if (ra === "i" && st[e.b] === 0 && next[e.b] === 0 && R() < p) next[e.b] = 1;
      if (rb === "i" && st[e.a] === 0 && next[e.a] === 0 && R() < p) next[e.a] = 1;
    }
    for (let i = 0; i < n; i++) {
      if (st[i] <= 0 || next[i] !== st[i]) continue;
      const role = chain[st[i]].role;
      if (R() < (ADVANCE[role] ?? 0.2)) {
        if (st[i] < chain.length - 1) next[i] = st[i] + 1;
        else if (ex.loop && lastRole !== "x") next[i] = 0;
      }
    }
    S.st = next;
    S.t++;
    return S;
  }

  // The card's still image: the first moment about 45% of nodes have been infected.
  function snapshot() {
    const S = init("sim:" + ex.id);
    const ever = new Uint8Array(n);
    for (let t = 0; t < 80; t++) {
      step(S);
      let k = 0;
      for (let i = 0; i < n; i++) { if (S.st[i] > 0) ever[i] = 1; k += ever[i]; }
      if (k >= 0.45 * n) break;
    }
    return S.st;
  }

  return { init, step, snapshot, chain };
}

export function glyph(ex) {
  const net = buildNetwork(ex);
  const sim = makeSim(ex, net);
  const root = svg("svg", { viewBox: `0 0 ${W} ${H}`, class: "glyph", "aria-hidden": "true", focusable: "false" });
  const gE = svg("g", { class: "g-edges" }), gN = svg("g", { class: "g-nodes" });
  root.append(gE, gN);
  for (const e of net.edges) {
    const a = net.nodes[e.a], b = net.nodes[e.b];
    gE.append(svg("line", { x1: a.x.toFixed(1), y1: a.y.toFixed(1), x2: b.x.toFixed(1), y2: b.y.toFixed(1), class: "e e-" + e.kind }));
  }
  const dots = net.nodes.map((p) => {
    const halo = svg("circle", { cx: p.x.toFixed(1), cy: p.y.toFixed(1), r: 8, class: "halo" });
    const dot = svg("circle", { cx: p.x.toFixed(1), cy: p.y.toFixed(1), r: 3.6, class: "n" });
    gN.append(halo, dot);
    return { halo, dot };
  });

  function paint(st) {
    st.forEach((s, i) => {
      const role = s < 0 ? (ex.side[-s - 1] || { role: "v" }).role : sim.chain[s].role;
      const { dot, halo } = dots[i];
      dot.style.fill = roleOf(role).color;
      dot.setAttribute("r", role === "i" ? 4.6 : 3.6);
      dot.classList.toggle("dead", role === "x");
      halo.style.fill = roleOf(role).color;
      halo.classList.toggle("on", role === "i");
    });
  }

  const still = sim.snapshot();
  paint(still);

  let timer = null, S = null, quiet = 0;
  function play() {
    if (timer || reducedMotion()) return;
    S = sim.init("sim:" + ex.id);
    paint(S.st);
    quiet = 0;
    timer = setInterval(() => {
      const before = S.st.join();
      sim.step(S);
      paint(S.st);
      quiet = S.st.join() === before ? quiet + 1 : 0;
      if (S.t > 80 || quiet > 12) { clearInterval(timer); timer = null; }
    }, 120);
  }
  function stop() {
    if (timer) { clearInterval(timer); timer = null; }
    paint(still);
  }
  return { node: root, play, stop, type: net.type };
}
