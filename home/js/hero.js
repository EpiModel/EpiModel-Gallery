// The header band: what the gallery is, where to start, and the live network model
// shared with www.epimodel.org (brand/js/brand.js).

import { el, liveModel } from "../../brand/js/brand.js";

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
  const live = liveModel({ caption: [el("strong", { text: "A stochastic network model, live. " }),
    "Ties form between nearby people and dissolve at random while an SEIRS infection passes across them. Immunity wanes, so infection returns in waves. The lines below the network track the share of people in each state. Every example in the gallery starts from this core and adds its own modules."] });
  return el("section", { class: "hero", "aria-label": "About the gallery" },
    el("div", { class: "hero-inner" }, text, live));
}
