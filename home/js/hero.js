// The header band: what the gallery is, where to start, and a live toy model of
// the thing every example extends. The toy is an SEIRS epidemic on a dynamic
// network: ties form between nearby nodes and dissolve at a constant rate, and
// infection passes across whatever ties exist at each step.

import { el, rng, ROLE, onWidth, reducedMotion } from "./common.js";

const P = {
  meanDegree: 2.6,     // target mean degree of the dynamic network
  duration: 25,        // mean tie duration, in steps
  reach: 0.25,         // ties form within this share of the panel width
  beta: 0.25,          // transmission probability per infectious tie per step
  latent: 5,           // mean steps in E
  infectious: 12,      // mean steps in I
  immunity: 60,        // mean steps in R before waning
  stepMs: 140,         // wall-clock time per step
  history: 220         // steps shown in the curve
};
const ORDER = ["s", "e", "i", "r"];
const LABEL = { s: "S", e: "E", i: "I", r: "R" };

function makeSim(n, w, h, seed) {
  const R = rng(seed);
  const nodes = [];
  for (let i = 0; i < n; i++) nodes.push({ hx: R() * w, hy: R() * h, ph: R() * 6.283, f: 0.3 + R() * 0.4, st: "s", t: -1e9 });
  // Spread the home positions out a little.
  for (let it = 0; it < 60; it++) {
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const a = nodes[i], b = nodes[j], dx = a.hx - b.hx, dy = a.hy - b.hy, d = Math.hypot(dx, dy) || 1;
      const min = Math.sqrt((w * h) / n) * 0.8;
      if (d < min) { const k = (min - d) / d * 0.25; a.hx += dx * k; a.hy += dy * k; b.hx -= dx * k; b.hy -= dy * k; }
    }
    for (const a of nodes) { a.hx = Math.max(10, Math.min(w - 10, a.hx)); a.hy = Math.max(10, Math.min(h - 10, a.hy)); }
  }
  const edges = new Map(), dying = [];
  const key = (a, b) => (a < b ? a * 4096 + b : b * 4096 + a);
  const near = nodes.map((a, i) => nodes.map((b, j) => j).filter((j) => j !== i && Math.hypot(a.hx - nodes[j].hx, a.hy - nodes[j].hy) < P.reach * w));
  const target = Math.round((n * P.meanDegree) / 2);
  const counts = [];
  let step = 0, quietSince = -1;

  function form(max) {
    let k = 0, guard = 0;
    while (edges.size < target && k < max && guard++ < 200) {
      const a = Math.floor(R() * n), c = near[a];
      if (!c.length) continue;
      const b = c[Math.floor(R() * c.length)], id = key(a, b);
      if (!edges.has(id)) { edges.set(id, { a, b, born: step }); k++; }
    }
  }
  function introduce(k) {
    for (let m = 0; m < k; m++) { const i = Math.floor(R() * n); if (nodes[i].st === "s") { nodes[i].st = "e"; nodes[i].t = step; } }
  }

  form(Infinity);
  for (const e of edges.values()) e.born = -100;
  introduce(2);
  counts.push(census());

  function census() {
    const c = { s: 0, e: 0, i: 0, r: 0 };
    for (const a of nodes) c[a.st]++;
    return c;
  }

  function tick() {
    step++;
    for (const [id, e] of edges) if (R() < 1 / P.duration) { edges.delete(id); dying.push({ ...e, died: step }); }
    form(6);
    const infect = [];
    for (const e of edges.values()) {
      const A = nodes[e.a], B = nodes[e.b];
      if (A.st === "i" && B.st === "s" && R() < P.beta) infect.push(B);
      else if (B.st === "i" && A.st === "s" && R() < P.beta) infect.push(A);
    }
    for (const a of nodes) {
      if (a.st === "e" && R() < 1 / P.latent) { a.st = "i"; a.t = step; }
      else if (a.st === "i" && R() < 1 / P.infectious) { a.st = "r"; a.t = step; }
      else if (a.st === "r" && R() < 1 / P.immunity) { a.st = "s"; a.t = step; }
    }
    for (const a of infect) if (a.st === "s") { a.st = "e"; a.t = step; }
    const c = census();
    // Reintroduce infection a while after it dies out, as an importation would.
    if (c.e + c.i === 0) { if (quietSince < 0) quietSince = step; if (step - quietSince > 25) { introduce(2); quietSince = -1; } }
    counts.push(c);
    if (counts.length > P.history) counts.shift();
    while (dying.length && step - dying[0].died > 3) dying.shift();
  }

  return { nodes, edges, dying, counts, tick, get step() { return step; } };
}

function liveModel() {
  const wrap = el("figure", { class: "live" });
  const head = el("div", { class: "live-head" });
  const legend = el("div", { class: "live-legend" });
  const cells = {};
  for (const s of ORDER) {
    const k = el("span", { class: "k" }); k.style.background = ROLE[s].color;
    const v = el("b", { text: "0" });
    cells[s] = v;
    legend.append(el("span", { class: "lg", title: ROLE[s].name }, k, el("span", { class: "lb", text: LABEL[s] }), v));
  }
  const btnPlay = el("button", { type: "button", class: "live-btn", "aria-label": "Pause the animation" });
  const btnNew = el("button", { type: "button", class: "live-btn", text: "New run", "aria-label": "Start a new run" });
  head.append(legend, el("div", { class: "live-ctl" }, btnPlay, btnNew));

  const net = el("canvas", { class: "live-net", "aria-hidden": "true" });
  const curve = el("canvas", { class: "live-curve", "aria-hidden": "true" });
  const cap = el("figcaption", {},
    el("strong", { text: "A stochastic network model, live. " }),
    "Ties form between nearby people and dissolve at random while an SEIRS infection passes across them. Immunity wanes, so infection returns in waves. The lines below the network track the share of people in each state. Every example in the gallery starts from this core and adds its own modules.");
  wrap.append(head, net, curve, cap);

  let W = 0, Hn = 0, Hc = 64, sim = null, runs = 0, playing = !reducedMotion(), visible = true, last = 0, raf = 0;
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const gN = net.getContext("2d"), gC = curve.getContext("2d");

  function size(w) {
    W = w; Hn = Math.round(Math.max(200, Math.min(300, w * 0.56)));
    for (const [c, h] of [[net, Hn], [curve, Hc]]) {
      c.width = Math.round(W * dpr); c.height = Math.round(h * dpr);
      c.style.width = "100%"; c.style.height = h + "px";
    }
    restart();
  }
  function restart() {
    const n = Math.round(Math.max(60, Math.min(130, (W * Hn) / 1500)));
    sim = makeSim(n, W, Hn, 20261017 + runs++);
    // Start a little way in, so the first frame already shows spread.
    for (let i = 0, k = playing ? 30 : P.history; i < k; i++) sim.tick();
    draw(performance.now());
    label();
  }
  function label() {
    btnPlay.textContent = playing ? "Pause" : "Play";
    btnPlay.setAttribute("aria-label", playing ? "Pause the animation" : "Play the animation");
  }

  function draw(now) {
    if (!sim) return;
    const t = now / 1000, frac = Math.min(1, (now - last) / P.stepMs);
    const pos = (a) => [a.hx + (playing ? 2.2 * Math.sin(t * a.f + a.ph) : 0), a.hy + (playing ? 2.2 * Math.cos(t * a.f * 0.8 + a.ph) : 0)];
    gN.setTransform(dpr, 0, 0, dpr, 0, 0);
    gN.clearRect(0, 0, W, Hn);
    gN.lineWidth = 1;
    const css = getComputedStyle(wrap);
    const edgeCol = css.getPropertyValue("--live-edge").trim() || "rgba(60,80,100,.45)";
    const drawEdge = (e, alpha) => {
      const [x1, y1] = pos(sim.nodes[e.a]), [x2, y2] = pos(sim.nodes[e.b]);
      gN.globalAlpha = alpha; gN.strokeStyle = edgeCol;
      gN.beginPath(); gN.moveTo(x1, y1); gN.lineTo(x2, y2); gN.stroke();
    };
    for (const e of sim.edges.values()) drawEdge(e, Math.min(1, (sim.step - e.born + frac) / 3));
    for (const e of sim.dying) drawEdge(e, Math.max(0, 1 - (sim.step - e.died + frac) / 3));
    gN.globalAlpha = 1;
    for (const a of sim.nodes) {
      const [x, y] = pos(a), col = ROLE[a.st].hex;
      if (a.st === "i") {
        gN.globalAlpha = 0.2; gN.fillStyle = col;
        gN.beginPath(); gN.arc(x, y, 9, 0, 6.283); gN.fill();
        gN.globalAlpha = 1;
      }
      const age = sim.step - a.t + frac;
      if (a.st === "e" && age < 4) {
        gN.globalAlpha = Math.max(0, 1 - age / 4); gN.strokeStyle = col; gN.lineWidth = 1.5;
        gN.beginPath(); gN.arc(x, y, 4 + age * 4, 0, 6.283); gN.stroke();
        gN.globalAlpha = 1; gN.lineWidth = 1;
      }
      gN.fillStyle = col;
      gN.beginPath(); gN.arc(x, y, a.st === "i" ? 4.6 : 3.4, 0, 6.283); gN.fill();
    }
    drawCurve(css);
  }

  function drawCurve(css) {
    const c = sim.counts, n = sim.nodes.length;
    gC.setTransform(dpr, 0, 0, dpr, 0, 0);
    gC.clearRect(0, 0, W, Hc);
    gC.strokeStyle = css.getPropertyValue("--live-grid").trim() || "#e1e0d9";
    gC.lineWidth = 1;
    gC.beginPath(); gC.moveTo(0, Hc - 0.5); gC.lineTo(W, Hc - 0.5); gC.stroke();
    const now = c[c.length - 1];
    for (const s of ORDER) cells[s].textContent = String(now[s]);
    if (c.length < 2) return;
    const x = (i) => (i / (P.history - 1)) * (W - 4) + 2, y = (v) => Hc - 3 - (v / n) * (Hc - 8);
    gC.lineWidth = 2; gC.lineJoin = "round";
    for (const s of ORDER) {
      gC.strokeStyle = ROLE[s].hex;
      gC.beginPath();
      c.forEach((d, i) => (i ? gC.lineTo(x(i), y(d[s])) : gC.moveTo(x(i), y(d[s]))));
      gC.stroke();
    }
  }

  function loop(now) {
    raf = 0;
    if (!playing || !visible) return;
    if (now - last >= P.stepMs) { sim.tick(); last = now; }
    draw(now);
    raf = requestAnimationFrame(loop);
  }
  function kick() { if (playing && visible && !raf) { last = performance.now(); raf = requestAnimationFrame(loop); } }

  btnPlay.addEventListener("click", () => { playing = !playing; label(); kick(); if (!playing) draw(performance.now()); });
  btnNew.addEventListener("click", () => { restart(); kick(); });
  document.addEventListener("visibilitychange", () => { visible = !document.hidden; kick(); });
  new IntersectionObserver((es) => { visible = es[0].isIntersecting && !document.hidden; kick(); }).observe(wrap);
  onWidth(wrap, (w) => { size(w); kick(); });
  label();
  return wrap;
}

export function hero(G, { onPath } = {}) {
  const exs = G.examples;
  const nFns = Object.keys(G.topics).length;
  const first = exs.find((e) => e.level === "Beginner") || exs[0];
  const pathBtn = el("a", { class: "hbtn ghost", href: "#explore", text: "See the learning path" });
  pathBtn.addEventListener("click", () => onPath && onPath());
  const text = el("div", { class: "hero-text" },
    el("div", { class: "kicker", text: "Template examples for EpiModel" }),
    el("h1", { text: "EpiModel Gallery" }),
    el("p", { class: "lead", text: "Annotated, runnable templates for extending EpiModel to model infectious disease over dynamic networks. Each example walks through the custom modules it adds, runs them, and ships the standalone scripts so you can adapt them to your own question." }),
    el("div", { class: "hero-cta" },
      el("a", { class: "hbtn", href: first.href, text: `Start here: ${first.short}` }),
      pathBtn),
    el("dl", { class: "hero-stats" },
      el("div", {}, el("dt", { text: "Examples" }), el("dd", { text: String(exs.length) })),
      el("div", {}, el("dt", { text: "EpiModel functions in use" }), el("dd", { text: String(nFns) })),
      el("div", {}, el("dt", { text: "Tested in CI against" }), el("dd", { text: G.epimodel ? `EpiModel ${G.epimodel}` : "EpiModel" }))));
  return el("section", { class: "hero", "aria-label": "About the gallery" },
    el("div", { class: "hero-inner" }, text, liveModel()));
}
