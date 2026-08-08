#!/usr/bin/env python3

"""
plot_taxonomy.py
────────────────
maniFasta post-build step — static taxonomic sunburst.

Reads an enriched metadata/manifest TSV (the tax_* columns written by
add_lineage_to_metadata.py) and renders a static, publication-ready sunburst
of the taxonomic composition of the protein database.

Backends, chosen by --out's extension:
  .svg  -> pure-Python SVG writer, STDLIB ONLY (no matplotlib/numpy). This is
           the default and always works on a bare cluster python3. SVG is
           vector, so it scales cleanly into manuscripts and slides.
  .png / .pdf -> rendered with matplotlib if it is importable; if matplotlib is
           absent, the tool warns and writes an .svg next to it instead, so a
           figure is always produced.

Rings, inner -> outer: domain/realm, phylum, class, order, family, genus,
species. Wedge angle is proportional to protein count; color is assigned by
top-level group (phylum by default). If a tax_n_proteins / count / n column is
present it is summed; otherwise each row counts as one protein.

Stdlib only for SVG output (csv, math). matplotlib is optional, for raster.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional

RANKS = ["domain_or_realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_LABEL = {"domain_or_realm": "domain/realm", "kingdom": "kingdom", "phylum": "phylum",
              "class": "class", "order": "order", "family": "family", "genus": "genus",
              "species": "species"}
# Pluralised titles for the color legend, keyed by the --color-by rank.
RANK_LEGEND = {"domain_or_realm": "Domains / realms", "kingdom": "Kingdoms",
               "phylum": "Phyla", "class": "Classes", "order": "Orders",
               "family": "Families", "genus": "Genera", "species": "Species"}

# Luminous-on-light categorical palette (distinct, print-legible), hex.
PALETTE = ["#2a9d8f", "#e76f51", "#e9c46a", "#7b6cd9", "#43aa8b", "#4d96d9",
           "#d96bb0", "#90a955", "#f3722c", "#577590", "#b5179e", "#bc6c25",
           "#0fa3b1", "#9e2a2b", "#3a86ff", "#8338ec"]
UNCLASS_COLOR = "#b9c2c4"


def info(msg: str) -> None:
    print(f"[INFO] {msg}", file=sys.stderr)


def warn(msg: str) -> None:
    print(f"[WARN] {msg}", file=sys.stderr)


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


# ── input ────────────────────────────────────────────────────────────────────

def read_tsv(path: Path):
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), [dict(r) for r in reader]


def detect_cols(header: List[str], prefix: str) -> Dict[str, Optional[str]]:
    low = {h.lower(): h for h in header}

    def find(*cands):
        for c in cands:
            if c.lower() in low:
                return low[c.lower()]
        return None

    cols = {}
    for r in RANKS:
        cols[r] = find(prefix + r, r) or (
            find(prefix + "domain", "domain", "superkingdom") if r == "domain_or_realm" else None)
    cols["weight"] = find(prefix + "n_proteins", "n_proteins", "count", "n", "proteins")
    return cols


class Node:
    __slots__ = ("name", "rank", "depth", "value", "children")

    def __init__(self, name, rank="", depth=0):
        self.name = name
        self.rank = rank
        self.depth = depth
        self.value = 0.0
        self.children: Dict[str, "Node"] = {}


def build_tree(rows, header, prefix, unclassified_label, max_depth):
    cols = detect_cols(header, prefix)
    if not cols["phylum"] and not cols["domain_or_realm"]:
        die(f"no taxonomy columns found (looked for {prefix}phylum, {prefix}domain_or_realm, ...). "
            f"Columns present: {header}")
    use_ranks = RANKS[:max_depth] if max_depth else RANKS

    def val(row, r):
        return (row[cols[r]] or "").strip() if cols[r] else ""

    # Keep only ranks that at least one row populates, in canonical order. This
    # drops globally-empty ranks (e.g. 'kingdom' for an all-bacterial DB, which
    # NCBI leaves blank) so they don't add a ring of carry-forward duplicates —
    # while any rank that IS present stays at its true, fixed slot.
    active = [r for r in use_ranks if any(val(row, r) for row in rows)]
    if not active:
        active = ["domain_or_realm"]

    root = Node("Database")
    n_rows = 0
    for row in rows:
        w = float(row[cols["weight"]] or 0) if cols["weight"] else 1.0
        if w <= 0:
            continue
        n_rows += 1

        # Read each active rank into its FIXED slot so rings align by true rank,
        # not by how many ranks happen to be populated for this row.
        vals = [val(row, r) for r in active]
        last = max((i for i, v in enumerate(vals) if v), default=-1)
        if last < 0:
            path = [(active[0], unclassified_label)]
        else:
            # Build a contiguous chain down to the deepest populated rank,
            # filling internal blanks by carrying the nearest known ancestor
            # forward for tree continuity. Gap nodes get an EMPTY rank so they
            # are never counted as a group / labelled at that rank (a carried
            # domain name must not surface as, e.g., a phylum in the legend).
            path = []
            carried = unclassified_label
            for i in range(last + 1):
                if vals[i]:
                    carried = vals[i]
                    path.append((active[i], carried))   # real node at its true rank
                else:
                    path.append(("", carried))          # gap: carried name, no rank

        node = root
        root.value += w
        for i, (rank, name) in enumerate(path):
            child = node.children.get(name)
            if child is None:
                child = Node(name, rank, i + 1)   # depth = ring index (1-based)
                node.children[name] = child
            child.value += w
            node = child
    return root, n_rows


# ── color ────────────────────────────────────────────────────────────────────

def top_group(node_path: List[Node], color_by: str) -> str:
    """Name of the node at the color_by rank in this lineage, or '' if that rank
    is absent for the lineage (so the wedge colours as unclassified rather than
    borrowing a deeper node's name and mislabelling the legend)."""
    for n in node_path:
        if n.rank == color_by:
            return n.name
    return ""


def make_color_map(root: Node, color_by: str, unclassified_label: str):
    groups = []

    def walk(n, path):
        for c in n.children.values():
            g = top_group(path + [c], color_by)
            if g and g != unclassified_label and g not in groups:
                groups.append(g)
            walk(c, path + [c])

    walk(root, [])
    cmap = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}
    cmap[unclassified_label] = UNCLASS_COLOR
    return cmap


def shade_hex(hex_color: str, depth: int, max_depth: int) -> str:
    """Lighten a hex color with increasing ring depth for a subtle gradient."""
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    f = 0.0 if max_depth <= 1 else min(0.45, 0.45 * (depth - 1) / max(1, max_depth - 1))
    r = int(round(r + (255 - r) * f))
    g = int(round(g + (255 - g) * f))
    b = int(round(b + (255 - b) * f))
    return f"#{r:02x}{g:02x}{b:02x}"


# ── ring-depth helper ────────────────────────────────────────────────────────

def measured_depth(root: Node) -> int:
    best = [0]

    def walk(n, d):
        best[0] = max(best[0], d)
        for c in n.children.values():
            walk(c, d + 1)

    walk(root, 0)
    return max(best[0], 1)


def xml_escape(s: str) -> str:
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;"))


# ── SVG backend (stdlib only) ────────────────────────────────────────────────

def draw_svg(root, *, color_by, cmap, unclassified_label, title, subtitle,
             label_min_pct, out_path, max_rings, label_ranks, label_style="callout"):
    total = root.value or 1.0
    ring_count = min(measured_depth(root), max_rings) if max_rings else measured_depth(root)
    ring_count = max(ring_count, 1)

    # Wider-than-tall canvas: reserve generous left/right margins for callout text.
    W, H = 1500, 1140
    CX, CY = W / 2.0, H * 0.545
    R = H * 0.40
    inner_hole = R * 0.18
    dr = (R - inner_hole) / ring_count

    parts: List[str] = []
    legend_groups: Dict[str, List] = {}
    callouts: List[dict] = []      # outer labels deferred for de-overlap

    def polar(r, theta):           # theta: radians from top (12 o'clock), clockwise
        return (CX + r * math.sin(theta), CY - r * math.cos(theta))

    def wedge_path(r_in, r_out, t0, t1):
        if (t1 - t0) >= 2 * math.pi - 1e-9:    # full ring -> annulus (even-odd)
            return (f"M{CX:.2f},{CY - r_out:.2f} "
                    f"A{r_out:.2f},{r_out:.2f} 0 1 1 {CX:.2f},{CY + r_out:.2f} "
                    f"A{r_out:.2f},{r_out:.2f} 0 1 1 {CX:.2f},{CY - r_out:.2f} Z "
                    f"M{CX:.2f},{CY - r_in:.2f} "
                    f"A{r_in:.2f},{r_in:.2f} 0 1 0 {CX:.2f},{CY + r_in:.2f} "
                    f"A{r_in:.2f},{r_in:.2f} 0 1 0 {CX:.2f},{CY - r_in:.2f} Z")
        large = 1 if (t1 - t0) > math.pi else 0
        p1 = polar(r_out, t0); p2 = polar(r_out, t1)
        p3 = polar(r_in, t1);  p4 = polar(r_in, t0)
        return (f"M{p1[0]:.2f},{p1[1]:.2f} "
                f"A{r_out:.2f},{r_out:.2f} 0 {large} 1 {p2[0]:.2f},{p2[1]:.2f} "
                f"L{p3[0]:.2f},{p3[1]:.2f} "
                f"A{r_in:.2f},{r_in:.2f} 0 {large} 0 {p4[0]:.2f},{p4[1]:.2f} Z")

    def add_inplace_label(name, rr, tm, rotate):
        x, y = polar(rr, tm)
        name = name if len(name) <= 26 else name[:25] + "…"
        if rotate:
            ang = math.degrees(math.atan2(-math.cos(tm), math.sin(tm)))
            if ang > 90 or ang < -90:
                ang += 180
            tf = f' transform="rotate({ang:.2f} {x:.2f} {y:.2f})"'
        else:
            tf = ""
        parts.append(
            f'<text x="{x:.2f}" y="{y:.2f}"{tf} text-anchor="middle" '
            f'dominant-baseline="central" font-size="11" fill="#16242a" '
            f'font-family="Helvetica,Arial,sans-serif">{xml_escape(name)}</text>')

    def recurse(node, depth, a0, a1, path):
        span = a1 - a0
        for child in sorted(node.children.values(), key=lambda c: -c.value):
            frac = span * (child.value / node.value) if node.value else 0
            c0, c1 = a0, a0 + frac
            a0 = c1
            if depth + 1 > ring_count:
                continue
            r_in = inner_hole + depth * dr
            r_out = r_in + dr
            t0, t1 = c0 * 2 * math.pi, c1 * 2 * math.pi
            grp = top_group(path + [child], color_by)
            fill = shade_hex(cmap.get(grp, UNCLASS_COLOR), child.depth, ring_count)
            parts.append(f'<path d="{wedge_path(r_in, r_out, t0, t1)}" '
                         f'fill="{fill}" stroke="white" stroke-width="0.8"/>')
            if child.rank == color_by:
                legend_groups.setdefault(grp, [cmap.get(grp, UNCLASS_COLOR), 0.0])
                legend_groups[grp][1] += child.value

            span_deg = (c1 - c0) * 360.0
            pct = 100 * child.value / total
            rank_ok = (not label_ranks) or (child.rank in label_ranks)
            if rank_ok and pct >= label_min_pct and span_deg >= 3.0:
                tm = (t0 + t1) / 2
                rr = (r_in + r_out) / 2
                is_outer = r_out >= 0.55 * R
                if label_style == "radial":
                    add_inplace_label(child.name, rr, tm, rotate=True)
                elif label_style == "horizontal":
                    add_inplace_label(child.name, rr, tm, rotate=False)
                else:  # callout: outer wedges -> leader lines; inner -> in-place horizontal
                    if is_outer:
                        ax, ay = polar(r_out, tm)
                        callouts.append({"name": child.name, "ax": ax, "ay": ay,
                                         "side": "R" if math.sin(tm) >= 0 else "L",
                                         "y": ay})
                    elif span_deg >= 6.0:
                        add_inplace_label(child.name, rr, tm, rotate=False)
            recurse(child, depth + 1, c0, c1, path + [child])

    recurse(root, 0, 0.0, 1.0, [])

    # ── callout layout: stack horizontal labels in left/right columns ────────
    def place(side_items, x_text, x_elbow, anchor):
        gap, lo, hi = 15.0, 70.0, H - 70.0
        side_items.sort(key=lambda d: d["ay"])
        for d in side_items:                       # initialize at anchor point
            d["y"] = d["ay"]
        for i in range(1, len(side_items)):        # push down
            if side_items[i]["y"] - side_items[i - 1]["y"] < gap:
                side_items[i]["y"] = side_items[i - 1]["y"] + gap
        if side_items and side_items[-1]["y"] > hi: # push up from bottom
            side_items[-1]["y"] = hi
            for i in range(len(side_items) - 2, -1, -1):
                if side_items[i + 1]["y"] - side_items[i]["y"] < gap:
                    side_items[i]["y"] = side_items[i + 1]["y"] - gap
        if side_items:
            side_items[0]["y"] = max(side_items[0]["y"], lo)
        for d in side_items:
            name = d["name"] if len(d["name"]) <= 34 else d["name"][:33] + "…"
            parts.append(
                f'<polyline points="{d["ax"]:.1f},{d["ay"]:.1f} '
                f'{x_elbow:.1f},{d["y"]:.1f} {x_text - (6 if anchor=="start" else -6):.1f},{d["y"]:.1f}" '
                f'fill="none" stroke="#9aa9ab" stroke-width="0.8"/>')
            parts.append(
                f'<text x="{x_text:.1f}" y="{d["y"]:.1f}" text-anchor="{anchor}" '
                f'dominant-baseline="central" font-size="11.5" fill="#16242a" '
                f'font-family="Helvetica,Arial,sans-serif">{xml_escape(name)}</text>')

    right = [d for d in callouts if d["side"] == "R"]
    left = [d for d in callouts if d["side"] == "L"]
    place(right, x_text=CX + R + 36, x_elbow=CX + R + 22, anchor="start")
    place(left, x_text=CX - R - 36, x_elbow=CX - R - 22, anchor="end")

    # center hub
    parts.append(f'<circle cx="{CX:.2f}" cy="{CY:.2f}" r="{inner_hole:.2f}" '
                 f'fill="#0c1418" stroke="#2a9d8f" stroke-width="1.4"/>')
    parts.append(f'<text x="{CX:.2f}" y="{CY - 3:.2f}" text-anchor="middle" '
                 f'font-size="20" font-weight="bold" fill="#39e6c4" '
                 f'font-family="monospace">{int(round(total)):,}</text>')
    parts.append(f'<text x="{CX:.2f}" y="{CY + 16:.2f}" text-anchor="middle" '
                 f'font-size="11" fill="#9fc2c4" font-family="monospace">proteins</text>')

    # title + subtitle (top-left)
    parts.append(f'<text x="22" y="40" font-size="23" font-weight="bold" fill="#16242a" '
                 f'font-family="Helvetica,Arial,sans-serif">{xml_escape(title)}</text>')
    if subtitle:
        parts.append(f'<text x="23" y="63" font-size="13" fill="#5a7177" '
                     f'font-family="Helvetica,Arial,sans-serif">{xml_escape(subtitle)}</text>')

    # legend: compact horizontal row across the top (keeps left/right free for callouts)
    items = sorted(legend_groups.items(), key=lambda kv: -kv[1][1])
    leg_title = RANK_LEGEND.get(color_by, color_by.replace("_", " ").capitalize())
    parts.append(f'<text x="22" y="92" font-size="12" font-weight="bold" fill="#16242a" '
                 f'font-family="Helvetica,Arial,sans-serif">{leg_title}:</text>')
    lx = 95.0
    for g, (col, val) in items:
        label = f"{g} ({100 * val / total:.0f}%)"
        parts.append(f'<rect x="{lx:.0f}" y="83" width="13" height="13" rx="3" fill="{col}"/>')
        parts.append(f'<text x="{lx + 19:.0f}" y="92" font-size="12" fill="#16242a" '
                     f'font-family="Helvetica,Arial,sans-serif">{xml_escape(label)}</text>')
        lx += 19 + 7.1 * len(label) + 22       # advance by approx text width

    parts.append(f'<text x="22" y="{H - 16}" font-size="10" fill="#7f9a9e" '
                 f'font-family="Helvetica,Arial,sans-serif">'
                 f'Rings inner→outer: domain · phylum · class · order · family · genus · species. '
                 f'Wedge angle ∝ protein count.</text>')

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}" fill-rule="evenodd">'
           f'<rect width="{W}" height="{H}" fill="white"/>'
           + "".join(parts) + "</svg>")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(svg, encoding="utf-8")
    info(f"Wrote {out_path}  (SVG, stdlib backend, label-style={label_style})")


# ── matplotlib backend (optional, for raster PNG/PDF) ────────────────────────

def draw_matplotlib(root, *, color_by, cmap, unclassified_label, title, subtitle,
                    label_min_pct, dpi, out_path, max_rings, label_ranks):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Wedge
    from matplotlib.lines import Line2D

    total = root.value or 1.0
    ring_count = min(measured_depth(root), max_rings) if max_rings else measured_depth(root)
    ring_count = max(ring_count, 1)

    fig = plt.figure(figsize=(11, 11), dpi=dpi)
    ax = fig.add_axes([0.02, 0.02, 0.74, 0.92], aspect="equal")
    ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.axis("off")
    inner_hole = 0.18
    dr = (1.0 - inner_hole) / ring_count
    START_DEG = 90.0
    legend_groups = {}

    def frac_to_deg(frac):
        return START_DEG - frac * 360.0

    def recurse(node, depth, a0, a1, path):
        span = a1 - a0
        for child in sorted(node.children.values(), key=lambda c: -c.value):
            frac = span * (child.value / node.value) if node.value else 0
            c0, c1 = a0, a0 + frac
            a0 = c1
            if depth + 1 > ring_count:
                continue
            r_in = inner_hole + depth * dr
            grp = top_group(path + [child], color_by)
            fill = shade_hex(cmap.get(grp, UNCLASS_COLOR), child.depth, ring_count)
            ax.add_patch(Wedge((0, 0), r_in + dr, frac_to_deg(c1), frac_to_deg(c0),
                               width=dr, facecolor=fill, edgecolor="white", linewidth=0.7))
            if child.rank == color_by:
                legend_groups.setdefault(grp, [cmap.get(grp, UNCLASS_COLOR), 0.0])
                legend_groups[grp][1] += child.value
            span_deg = (c1 - c0) * 360.0
            pct = 100 * child.value / total
            rank_ok = (not label_ranks) or (child.rank in label_ranks)
            if rank_ok and pct >= label_min_pct and span_deg >= 3.0:
                mid_deg = frac_to_deg((c0 + c1) / 2)
                rr = (r_in + r_in + dr) / 2
                rad = math.radians(mid_deg)
                flip = 90 < (mid_deg % 360) < 270
                rot = mid_deg + 180 if flip else mid_deg
                name = child.name if len(child.name) <= 26 else child.name[:25] + "…"
                ax.text(rr * math.cos(rad), rr * math.sin(rad), name,
                        ha="right" if flip else "left", va="center",
                        rotation=rot, rotation_mode="anchor", fontsize=7.2, color="#16242a")
            recurse(child, depth + 1, c0, c1, path + [child])

    recurse(root, 0, 0.0, 1.0, [])
    ax.add_patch(plt.Circle((0, 0), inner_hole, facecolor="#0c1418", edgecolor="#2a9d8f", linewidth=1.2))
    ax.text(0, 0.045, f"{int(round(total)):,}", ha="center", va="center",
            fontsize=15, color="#39e6c4", fontweight="bold", family="monospace")
    ax.text(0, -0.05, "proteins", ha="center", va="center", fontsize=8.5,
            color="#9fc2c4", family="monospace")
    fig.text(0.02, 0.965, title, fontsize=17, fontweight="bold", color="#16242a")
    if subtitle:
        fig.text(0.02, 0.94, subtitle, fontsize=10, color="#5a7177")
    items = sorted(legend_groups.items(), key=lambda kv: -kv[1][1])
    handles = [Line2D([0], [0], marker="s", linestyle="none", markersize=10,
                      markerfacecolor=col, markeredgecolor="none",
                      label=f"{g}  ({100 * val / total:.0f}%)")
               for g, (col, val) in items]
    if handles:
        leg_title = RANK_LEGEND.get(color_by, color_by.replace("_", " ").capitalize())
        ax2 = fig.add_axes([0.77, 0.02, 0.22, 0.92]); ax2.axis("off")
        ax2.legend(handles=handles, loc="center left", frameon=False,
                   title=leg_title, fontsize=9, title_fontsize=10,
                   handletextpad=0.6, labelspacing=0.8)
    fig.text(0.02, 0.012,
             "Rings inner→outer: domain · phylum · class · order · family · genus · species. "
             "Wedge angle ∝ protein count.", fontsize=7.5, color="#7f9a9e")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    info(f"Wrote {out_path}  (matplotlib backend)")


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True, type=Path,
                    help="Enriched metadata/manifest TSV (tax_* columns).")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output image; format from extension (.svg [stdlib], .png/.pdf [matplotlib]).")
    ap.add_argument("--prefix", default="tax_", help="Lineage column prefix (default tax_).")
    ap.add_argument("--color-by", default="phylum", choices=["phylum", "domain_or_realm"],
                    dest="color_by", help="Top-level group for color/legend (default phylum).")
    ap.add_argument("--max-depth", type=int, default=0,
                    help="Limit rings to this many ranks (e.g. 6 stops at genus; 0 = all 7).")
    ap.add_argument("--label-min-pct", type=float, default=1.5, dest="label_min_pct",
                    help="Only label wedges at least this %% of total (default 1.5).")
    ap.add_argument("--label-ranks", default="", dest="label_ranks",
                    help="Comma-separated ranks to label, e.g. 'phylum,genus,species'. "
                         "Empty = label every rank that passes the size gate.")
    ap.add_argument("--label-style", default="callout", dest="label_style",
                    choices=["callout", "horizontal", "radial"],
                    help="callout: leader lines to stacked horizontal labels for outer rings "
                         "(default, avoids clipping/overrun); horizontal: upright text in each "
                         "wedge; radial: text rotated along the radius. SVG output only.")
    ap.add_argument("--unclassified-label", default="Unclassified",
                    help="Label for rows with no lineage (default Unclassified).")
    ap.add_argument("--title", default="Taxonomic composition of the protein database")
    ap.add_argument("--subtitle", default="", help="Optional subtitle (run label / date).")
    ap.add_argument("--dpi", type=int, default=200, help="Raster DPI for PNG (matplotlib only).")
    args = ap.parse_args()

    header, rows = read_tsv(args.metadata)
    max_depth = args.max_depth if args.max_depth and args.max_depth > 0 else 0
    root, n_rows = build_tree(rows, header, args.prefix, args.unclassified_label, max_depth)
    info(f"Rows used: {n_rows} | total weight: {int(root.value):,} | "
         f"top-level groups: {len(root.children)}")

    cmap = make_color_map(root, args.color_by, args.unclassified_label)
    max_rings = max_depth if max_depth else len(RANKS)
    label_ranks = {r.strip() for r in args.label_ranks.split(",") if r.strip()}
    subtitle = args.subtitle or f"{int(root.value):,} proteins · colored by {RANK_LABEL[args.color_by]}"
    common = dict(color_by=args.color_by, cmap=cmap, unclassified_label=args.unclassified_label,
                  title=args.title, subtitle=subtitle, label_min_pct=args.label_min_pct,
                  max_rings=max_rings, label_ranks=label_ranks)

    ext = args.out.suffix.lower()
    if ext == ".svg":
        draw_svg(root, out_path=args.out, label_style=args.label_style, **common)
    elif ext in (".png", ".pdf"):
        try:
            draw_matplotlib(root, out_path=args.out, dpi=args.dpi, **common)
        except ImportError:
            alt = args.out.with_suffix(".svg")
            warn(f"matplotlib not available for {ext} output; writing {alt} instead "
                 f"(stdlib SVG). Install matplotlib for raster output, or use --out *.svg.")
            draw_svg(root, out_path=alt, label_style=args.label_style, **common)
    else:
        die(f"unsupported output extension {ext!r}; use .svg (stdlib), .png or .pdf (matplotlib).")


if __name__ == "__main__":
    main()
