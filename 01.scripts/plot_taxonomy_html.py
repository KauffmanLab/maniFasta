#!/usr/bin/env python3
"""
plot_taxonomy_html.py
─────────────────────
Interactive, self-contained HTML sunburst of a maniFasta lineage manifest.

Companion to plot_taxonomy.py — it imports that module's tree builder and
colour logic, so the HTML and the static SVG always show the same numbers.

    python3 plot_taxonomy_html.py \
        --metadata AllOralsDB.v2026.214_S.manifest.lineage.tsv \
        --out      AllOralsDB.v2026.214_S.taxonomy_sunburst.html \
        --title    "AllOralsDB v2026.214_S"

Output is ONE .html file with no external dependencies.

Interactions
    hover        tooltip: name, rank, sequence count, % of database, lineage
    click wedge  zoom into that clade (angles rescale to fill the circle)
    click centre zoom out one level
    breadcrumb   click any ancestor to jump back
    search box   dims everything whose lineage doesn't match

stdlib only.
"""

import sys
import json
import argparse
import importlib.util
from pathlib import Path


def info(m): print(f"[INFO] {m}", file=sys.stderr)
def die(m):
    print(f"ERROR: {m}", file=sys.stderr)
    raise SystemExit(1)


def load_plot_taxonomy(path: Path):
    """Import plot_taxonomy.py from an explicit path (it isn't a package)."""
    if not path.is_file():
        die(f"plot_taxonomy.py not found at {path} — pass --plot-taxonomy")
    spec = importlib.util.spec_from_file_location("plot_taxonomy", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for need in ("read_tsv", "build_tree", "make_color_map",
                 "top_group", "shade_hex", "measured_depth", "UNCLASS_COLOR"):
        if not hasattr(mod, need):
            die(f"{path} has no '{need}' — version mismatch?")
    return mod


def tree_to_json(pt, root, *, color_by, cmap, ring_count):
    """Serialise the Node tree to plain dicts, precomputing each wedge's fill."""
    def walk(node, path, depth):
        d = {
            "n": node.name,
            "r": node.rank or "",
            "v": node.value,
        }
        if depth > 0:
            grp = pt.top_group(path, color_by)
            d["c"] = pt.shade_hex(cmap.get(grp, pt.UNCLASS_COLOR),
                                  node.depth, ring_count)
        if depth < ring_count and node.children:
            kids = sorted(node.children.values(), key=lambda c: -c.value)
            d["ch"] = [walk(k, path + [k], depth + 1) for k in kids]
        return d
    return walk(root, [], 0)


HTML = r"""<!DOCTYPE html>
<meta charset="utf-8">
<title>__TITLE__ — taxonomic composition</title>
<style>
  :root { --ink:#1c2426; --muted:#68787c; --rule:#e2e8ea; }
  * { box-sizing:border-box; }
  body { margin:0; padding:24px 28px; background:#fff; color:var(--ink);
         font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif; }
  h1 { font-size:19px; font-weight:600; margin:0 0 2px; letter-spacing:-.01em; }
  .sub { color:var(--muted); font-size:13px; margin-bottom:14px; }
  .bar { display:flex; align-items:center; gap:14px; flex-wrap:wrap;
         padding-bottom:10px; border-bottom:1px solid var(--rule); margin-bottom:12px; }
  #crumb { font-size:13px; color:var(--muted); }
  #crumb span { cursor:pointer; color:#2a7d8f; }
  #crumb span:hover { text-decoration:underline; }
  #crumb b { color:var(--ink); font-weight:600; }
  #q { padding:5px 9px; border:1px solid var(--rule); border-radius:5px;
       font-size:13px; width:210px; font-family:inherit; }
  #q:focus { outline:2px solid #cfe3e6; outline-offset:-1px; }
  .wrap { display:flex; gap:26px; align-items:flex-start; flex-wrap:wrap; }
  svg { flex:0 0 auto; }
  path.w { stroke:#fff; stroke-width:.7; cursor:pointer; }
  path.w:hover { stroke:#1c2426; stroke-width:1.4; }
  path.dim { opacity:.13; }
  #tip { position:fixed; pointer-events:none; opacity:0; transition:opacity .08s;
         background:rgba(20,28,30,.96); color:#fff; padding:9px 11px; border-radius:6px;
         font-size:12.5px; max-width:430px; z-index:9; box-shadow:0 4px 14px rgba(0,0,0,.22); }
  #tip .nm { font-weight:600; font-size:13.5px; }
  #tip .rk { color:#9fb6ba; font-style:italic; }
  #tip .ln { color:#c7d6d9; font-size:11.5px; margin-top:5px;
             border-top:1px solid rgba(255,255,255,.16); padding-top:5px;
             word-break:break-word; }
  #legend { flex:1 1 210px; min-width:200px; font-size:12.5px; }
  #legend div { display:flex; align-items:center; gap:7px; padding:1.5px 0; }
  #legend i { width:11px; height:11px; border-radius:2px; flex:0 0 auto; }
  #ctr { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }
  .hint { color:var(--muted); font-size:12px; margin-top:10px; }
</style>

<h1>__TITLE__</h1>
<div class="sub">__SUBTITLE__</div>

<div class="bar">
  <div id="crumb"></div>
  <input id="q" placeholder="highlight taxon…" autocomplete="off">
</div>

<div class="wrap">
  <svg id="sb" width="760" height="760" viewBox="0 0 760 760"></svg>
  <div id="legend"></div>
</div>
<div class="hint">Click a wedge to zoom · click the centre to go back · hover for full label</div>
<div id="tip"></div>

<script>
const DATA = __DATA__, LEGEND = __LEGEND__, RINGS = __RINGS__, UNIT = "__UNIT__";
const SVG="http://www.w3.org/2000/svg", CX=380, CY=380, R=352, HOLE=64;
const TOTAL = DATA.v, dr = (R-HOLE)/RINGS;
let focus = DATA, focusPath = [DATA], query = "";

const svg = document.getElementById("sb"), tip = document.getElementById("tip");
const fmt = n => Math.round(n).toLocaleString();
const pol = (r,t) => [CX + r*Math.sin(t), CY - r*Math.cos(t)];

function wedge(ri, ro, t0, t1){
  if (t1-t0 >= 2*Math.PI-1e-9){
    const p=(r)=>`M${CX},${CY-r} A${r},${r} 0 1 1 ${CX-.01},${CY-r} Z`;
    return p(ro)+p(ri);
  }
  const big = (t1-t0) > Math.PI ? 1 : 0;
  const [x0,y0]=pol(ro,t0), [x1,y1]=pol(ro,t1), [x2,y2]=pol(ri,t1), [x3,y3]=pol(ri,t0);
  return `M${x0.toFixed(2)},${y0.toFixed(2)} A${ro},${ro} 0 ${big} 1 ${x1.toFixed(2)},${y1.toFixed(2)}`
       + ` L${x2.toFixed(2)},${y2.toFixed(2)} A${ri},${ri} 0 ${big} 0 ${x3.toFixed(2)},${y3.toFixed(2)} Z`;
}

function matches(node, lineage){
  if (!query) return true;
  const q = query.toLowerCase();
  return lineage.concat([node.n]).some(s => s.toLowerCase().includes(q));
}

function draw(){
  svg.textContent = "";
  const frag = document.createDocumentFragment();

  // recurse — angle logic mirrors plot_taxonomy.py draw_svg()
  (function rec(node, depth, a0, a1, lineage){
    if (!node.ch || depth >= RINGS) return;
    const span = a1 - a0;
    let cur = a0;
    for (const ch of node.ch){
      const frac = node.v ? span * (ch.v / node.v) : 0;
      const c0 = cur, c1 = cur + frac;
      cur = c1;
      if (c1 - c0 < 1e-6) continue;
      const ri = HOLE + depth*dr, ro = ri + dr;
      const p = document.createElementNS(SVG, "path");
      p.setAttribute("d", wedge(ri, ro, c0*2*Math.PI, c1*2*Math.PI));
      p.setAttribute("fill", ch.c || "#b9c2c4");
      p.setAttribute("class", "w" + (matches(ch, lineage) ? "" : " dim"));
      const ln = lineage.concat([ch.n]);
      p.addEventListener("mousemove", e => {
        tip.style.opacity = 1;
        tip.style.left = Math.min(e.clientX+14, innerWidth-450) + "px";
        tip.style.top  = (e.clientY+16) + "px";
        tip.innerHTML = `<div class="nm">${esc(ch.n)}</div>`
          + (ch.r ? `<div class="rk">${esc(ch.r.replace(/_/g," "))}</div>` : "")
          + `<div>${fmt(ch.v)} ${UNIT} · ${(100*ch.v/TOTAL).toFixed(2)}% of database</div>`
          + `<div class="ln">${ln.map(esc).join(" › ")}</div>`;
      });
      p.addEventListener("mouseleave", () => tip.style.opacity = 0);
      p.addEventListener("click", ev => {
        ev.stopPropagation();
        if (ch.ch){ focus = ch; focusPath = focusPath.concat([ch]); tip.style.opacity=0; draw(); }
      });
      frag.appendChild(p);
      rec(ch, depth+1, c0, c1, ln);
    }
  })(focus, 0, 0, 1, focusPath.slice(1).map(n => n.n));

  svg.appendChild(frag);

  // centre disc = zoom out
  const c = document.createElementNS(SVG, "circle");
  c.setAttribute("cx",CX); c.setAttribute("cy",CY); c.setAttribute("r",HOLE-2);
  c.setAttribute("fill","#f4f7f8"); c.setAttribute("stroke","#dde5e7");
  c.style.cursor = focusPath.length>1 ? "pointer" : "default";
  c.addEventListener("click", () => {
    if (focusPath.length>1){ focusPath.pop(); focus = focusPath[focusPath.length-1]; draw(); }
  });
  svg.appendChild(c);

  const t1 = document.createElementNS(SVG,"text");
  t1.setAttribute("x",CX); t1.setAttribute("y",CY-4);
  t1.setAttribute("text-anchor","middle"); t1.setAttribute("id","ctr");
  t1.setAttribute("font-size","15"); t1.setAttribute("font-weight","600");
  t1.textContent = fmt(focus.v);
  svg.appendChild(t1);

  const t2 = document.createElementNS(SVG,"text");
  t2.setAttribute("x",CX); t2.setAttribute("y",CY+13);
  t2.setAttribute("text-anchor","middle");
  t2.setAttribute("font-size","10.5"); t2.setAttribute("fill","#68787c");
  t2.textContent = focusPath.length>1 ? "← back" : UNIT;
  svg.appendChild(t2);

  crumbs();
}

function esc(s){ return String(s).replace(/[&<>]/g, m => ({"&":"&amp;","<":"&lt;",">":"&gt;"}[m])); }

function crumbs(){
  const el = document.getElementById("crumb");
  el.innerHTML = focusPath.map((n,i) =>
    i === focusPath.length-1 ? `<b>${esc(n.n)}</b>` : `<span data-i="${i}">${esc(n.n)}</span>`
  ).join(" › ");
  el.querySelectorAll("span").forEach(s => s.addEventListener("click", () => {
    const i = +s.dataset.i;
    focusPath = focusPath.slice(0, i+1);
    focus = focusPath[i];
    draw();
  }));
}

document.getElementById("q").addEventListener("input", e => {
  query = e.target.value.trim(); draw();
});

const lg = document.getElementById("legend");
lg.innerHTML = LEGEND.map(([name,color,val]) =>
  `<div><i style="background:${color}"></i>${esc(name)} <span style="color:#68787c">`
  + `${(100*val/TOTAL).toFixed(1)}%</span></div>`).join("");

draw();
</script>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True, type=Path,
                    help="*.manifest.lineage.tsv from the build")
    ap.add_argument("--out", required=True, type=Path, help="output .html")
    ap.add_argument("--plot-taxonomy", type=Path, default=None,
                    help="path to plot_taxonomy.py (default: alongside this script)")
    ap.add_argument("--prefix", default="tax_")
    ap.add_argument("--color-by", default="phylum",
                    choices=["phylum", "domain_or_realm"], dest="color_by")
    ap.add_argument("--max-depth", type=int, default=0, dest="max_depth")
    ap.add_argument("--max-rings", type=int, default=8, dest="max_rings")
    ap.add_argument("--unclassified-label", default="Unclassified",
                    dest="unclassified_label")
    ap.add_argument("--title", default="Taxonomic composition")
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--unit", default="sequences")
    args = ap.parse_args()

    pt_path = args.plot_taxonomy or (Path(__file__).resolve().parent / "plot_taxonomy.py")
    pt = load_plot_taxonomy(pt_path)

    if not args.metadata.is_file():
        die(f"not found: {args.metadata}")

    header, rows = pt.read_tsv(args.metadata)
    info(f"rows: {len(rows):,}")

    max_depth = args.max_depth if args.max_depth and args.max_depth > 0 else 0
    root, n_rows = pt.build_tree(rows, header, args.prefix,
                                 args.unclassified_label, max_depth)
    cmap = pt.make_color_map(root, args.color_by, args.unclassified_label)

    ring_count = pt.measured_depth(root)
    if args.max_rings:
        ring_count = min(ring_count, args.max_rings)
    ring_count = max(ring_count, 1)
    info(f"tree total: {root.value:,.0f} | rings: {ring_count}")

    tree = tree_to_json(pt, root, color_by=args.color_by,
                        cmap=cmap, ring_count=ring_count)

    # legend: aggregate node values at the colour-by rank
    legend = {}

    def collect(node, path):
        for ch in node.children.values():
            p = path + [ch]
            if ch.rank == args.color_by:
                grp = pt.top_group(p, args.color_by)
                e = legend.setdefault(grp, [cmap.get(grp, pt.UNCLASS_COLOR), 0.0])
                e[1] += ch.value
            collect(ch, p)
    collect(root, [])
    leg = sorted(([k, v[0], v[1]] for k, v in legend.items()),
                 key=lambda x: -x[2])[:24]

    html = (HTML
            .replace("__TITLE__", args.title)
            .replace("__SUBTITLE__", args.subtitle or f"{int(root.value):,} sequences")
            .replace("__UNIT__", args.unit)
            .replace("__RINGS__", str(ring_count))
            .replace("__LEGEND__", json.dumps(leg))
            .replace("__DATA__", json.dumps(tree, separators=(",", ":"))))

    args.out.write_text(html, encoding="utf-8")
    info(f"wrote {args.out}  ({args.out.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
