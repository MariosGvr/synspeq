"""Standalone HTML report (no external assets, inline SVG) from a report.json dictionary."""
from __future__ import annotations
import datetime as dt
from html import escape as E
from typing import Dict, List, Optional

CSS = """
:root{--bg:#fff;--ink:#1b1f24;--mute:#5b6570;--line:#d9dee3;--card:#f6f8fa;--warn:#fff4d6;--warnl:#e3b341;--c1:#2f6f9f;--c2:#4c9f70;--c3:#c9a227;--c4:#b5523b;--c5:#8a6fb0;--c6:#9aa5b1}
@media (prefers-color-scheme:dark){:root{--bg:#14181c;--ink:#e6e9ec;--mute:#9aa5b1;--line:#2c333a;--card:#1b2026;--warn:#3a3217;--warnl:#a3822a}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
main{max-width:980px;margin:0 auto;padding:24px 16px 64px}h1{font-size:26px;margin:0 0 4px}h2{font-size:19px;margin:34px 0 8px;border-bottom:1px solid var(--line);padding-bottom:4px}
.sub{color:var(--mute)}table{border-collapse:collapse;width:100%;margin:8px 0}th,td{text-align:left;padding:5px 10px;border-bottom:1px solid var(--line);font-variant-numeric:tabular-nums}
th{font-weight:600;color:var(--mute);font-size:13px}td.n{text-align:right}.banner{background:var(--warn);border-left:4px solid var(--warnl);padding:10px 14px;margin:14px 0;border-radius:4px}
.cards{display:flex;flex-wrap:wrap;gap:10px;margin:10px 0}.card{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:10px 14px;min-width:150px}
.card b{display:block;font-size:22px}.card span{color:var(--mute);font-size:12px}.note{color:var(--mute);font-size:13px}.legend span{display:inline-block;margin-right:14px;font-size:13px}
.legend i{display:inline-block;width:10px;height:10px;margin-right:5px;border-radius:2px}svg text{fill:var(--ink);font-size:11px}ul{margin:6px 0 6px 20px;padding:0}
"""


def pct(x: Optional[float], d: int = 1) -> str:
    return "n/a" if x is None else f"{100 * x:.{d}f}%"


def num(x) -> str:
    return "n/a" if x is None else (f"{x:,}" if isinstance(x, int) else f"{x:,.1f}")


def card(value: str, label: str) -> str:
    return f'<div class="card"><b>{E(value)}</b><span>{E(label)}</span></div>'


def table(head: List[str], rows: List[List[str]], numeric_from: int = 1) -> str:
    h = "".join(f"<th>{E(c)}</th>" for c in head)
    body = "".join("<tr>" + "".join(f'<td class="{"n" if i >= numeric_from else ""}">{E(str(v))}</td>' for i, v in enumerate(r)) + "</tr>" for r in rows)
    return f"<table><tr>{h}</tr>{body}</table>"


def stacked_bar(parts: List[tuple], total: int, width: int = 900) -> str:
    """parts: [(label, count, colour_var)] -> one horizontal stacked bar with a legend."""
    if total <= 0:
        return ""
    x, rects, leg = 0.0, [], []
    for label, n, col in parts:
        if n <= 0:
            continue
        w = width * n / total
        rects.append(f'<rect x="{x:.1f}" y="0" width="{max(w, 1):.1f}" height="26" fill="var({col})"><title>{E(label)}: {n:,} ({100 * n / total:.1f}%)</title></rect>')
        leg.append(f'<span><i style="background:var({col})"></i>{E(label)} {n:,} ({100 * n / total:.1f}%)</span>')
        x += w
    return f'<svg viewBox="0 0 {width} 26" width="100%" role="img" aria-label="outcome shares">{"".join(rects)}</svg><div class="legend">{"".join(leg)}</div>'


def bars(items: List[tuple], width: int = 900, maxv: Optional[float] = None, fmt=lambda v: f"{v:.1f}") -> str:
    """items: [(label, value, colour_var)] -> horizontal bars."""
    if not items:
        return ""
    m = maxv or max(v for _, v, _ in items) or 1
    h = 24
    out = []
    for i, (lab, v, col) in enumerate(items):
        w = (width - 170) * (v / m)
        out.append(f'<text x="0" y="{i * h + 16}">{E(lab)}</text><rect x="160" y="{i * h + 4}" width="{max(w, 1):.1f}" height="14" fill="var({col})" rx="2"><title>{E(lab)}: {fmt(v)}</title></rect>'
                   f'<text x="{165 + w:.1f}" y="{i * h + 16}">{fmt(v)}</text>')
    return f'<svg viewBox="0 0 {width} {len(items) * h + 4}" width="100%">{"".join(out)}</svg>'


def render_html(rep: dict, out_path: str, title: str = "SYNSPEQ report") -> str:
    pop = rep.get("population", {})
    parts: List[str] = [f"<h1>{E(title)}</h1>",
                        f'<p class="sub">Population <b>{E(str(pop.get("label")))}</b> · SYNSPEQ {E(str(rep.get("synspeq_version")))} · '
                        f'identity policy {E(str(rep.get("identity_policy")))} · {dt.date.today().isoformat()}</p>']
    if rep.get("banner"):
        parts.append('<div class="banner"><b>Not fully assessed</b><ul>' + "".join(f"<li>{E(b)}</li>" for b in rep["banner"]) + "</ul></div>")
    if rep.get("real_queries_pending"):
        parts.append(f'<div class="banner">{rep["real_queries_pending"]:,} intermediates have not been looked up in REAL yet; Pillar 3 is reported as a bound until '
                     f'<code>synspeq real-scan</code> and <code>real-merge</code> have run.</div>')

    # input audit
    parts.append("<h2>Input</h2>")
    parts.append('<div class="cards">' + card(num(pop.get("n_input_records")), "records read") + card(num(pop.get("n_valid_records")), "valid") +
                 card(num(pop.get("n_invalid_records")), "invalid") + card(num(pop.get("n_distinct_molecules")), "distinct molecules") +
                 card(pct(pop.get("duplicate_rate")), "duplicate rate (raw rows)") + "</div>")
    if pop.get("invalid_reasons"):
        parts.append(table(["invalid reason", "records"], [[k, v] for k, v in pop["invalid_reasons"].items()]))
    note = f"Distinct molecule = {'(pocket, identity)' if pop.get('pocket_field') else 'identity'}; "
    note += f"window: {pop['window']}. " if pop.get("window") else ""
    parts.append(f'<p class="note">{E(note)}Identity is the RDKit canonical isomeric SMILES as read; no cleanup is applied. Multi-fragment molecules: {num(pop.get("n_multi_fragment"))}.</p>')

    p1 = rep.get("P1")
    if p1:
        c = p1["counts"]
        n = p1["n_distinct"]
        parts.append("<h2>Pillar 1: route recovery</h2>")
        parts.append('<div class="cards">' + card(pct(p1["find_ge1"]), "Find ≥1 (found, needs a reaction)") + card(pct(p1["find"]), "Find (found)") +
                     card(pct(p1["direct_purchase_share_of_found"]), "direct purchase, share of found") + card(pct(p1["find_raw"]), "Find over raw rows") +
                     card(pct(p1["non_stock_route_rate"]), "route rate among non-stock molecules") + "</div>")
        parts.append(stacked_bar([("direct purchase", c["direct_purchase"], "--c2"), ("routed", c["routed"], "--c1"),
                                  ("found, no saved route", c["found_no_saved_route"], "--c3"), ("unresolved", c["unresolved"], "--c4"),
                                  ("execution error", c["execution_error"], "--c5"), ("pending", c["pending"], "--c6")], n))
        if not p1["complete"]:
            parts.append('<p class="note"><b>Incomplete population:</b> pending or errored molecules are shown; the partition does not cover every molecule.</p>')
        parts.append('<p class="note">Unresolved = not recovered under these settings (a search that ended without a route), not proven unsynthesisable.</p>')

    p2 = rep.get("P2")
    if p2:
        parts.append("<h2>Pillar 2: chemical support (retrieved precedent)</h2>")
        rows = []
        for lab, k in (("primary route", "primary"), ("any saved route", "any_saved"), ("primary, identical records excluded", "primary_analogue"), ("any saved, identical records excluded", "any_saved_analogue")):
            b = p2[k]
            rows.append([lab, num(b["supported"]), num(b["indeterminate"]), pct(b["conditional"]), pct(b["unconditional"])])
        parts.append(table(["scope", "supported", "indeterminate", f"of routed ({p2['n_routed']:,})", f"of all distinct ({p2['n_distinct']:,})"], rows))
        L = p2["by_primary_route_length"]
        parts.append("<p class=\"note\">Support by primary-route length (longer routes have more steps to fail; compare tools at the same length):</p>")
        parts.append(bars([(f"{k} step{'s' if k != '1' else ''} (n={v['n']:,})", 100 * v["primary"], "--c1") for k, v in L.items() if v["primary"] is not None], maxv=100, fmt=lambda v: f"{v:.1f}%"))
        parts.append("<ul>" + "".join(f"<li>{E(x)}</li>" for x in p2.get("notes", [])) + "</ul>")

    p3 = rep.get("P3")
    if p3 and p3.get("n_routed"):
        parts.append("<h2>Pillar 3: deliverability</h2>")
        t = p3["tiers_available"]
        od_status = p3.get("needs_on_demand_status", "value")
        od_card = ("n/a" if p3.get("needs_on_demand") is None else
                   (("at least " if od_status == "bound" else "") + pct(p3["needs_on_demand"])))
        parts.append('<div class="cards">' + card(num(p3["median_residual_steps"]), "median residual steps" + (" (upper bound)" if p3.get("frontier_status") == "bound" else "")) +
                     card(num(p3["median_route_steps"]), "median route steps") +
                     card(pct(p3["frontier_found"]), "route has a purchasable intermediate" + (" (lower bound)" if p3.get("frontier_status") == "bound" else "")) +
                     card(od_card, "frontier is on-demand (REAL)" + (" (not assessed: needs both tiers)" if od_status == "not_assessed" else "")) + "</div>")
        used = " and ".join(x for x, ok in (("stocked catalogue", t["stocked"]), ("on-demand catalogue", t["on_demand"])) if ok)
        parts.append(f'<p class="note">Frontier defined by: {E(used)}. Strict matching (canonical isomeric SMILES). '
                     "'None' means not in the configured catalogues, not unbuyable anywhere. Describes routed molecules only.</p>")
        if p3.get("bound"):
            parts.append(f'<div class="banner">{E(p3["bound"])}</div>')
        d = p3.get("residual_steps_counts")
        if d:
            parts.append(bars([(f"{k} residual step{'s' if str(k) != '1' else ''}", v, "--c1") for k, v in sorted(d.items(), key=lambda kv: int(kv[0]))[:12]], fmt=lambda v: f"{int(v):,}"))

    p4 = rep.get("P4")
    if p4:
        parts.append("<h2>Pillar 4: coverage</h2>")
        rows = []
        for lab, k in (("routed molecules", "routed"), ("all distinct molecules", "all_distinct")):
            b = p4[k]
            rows.append([lab, num(b["n_molecules"]), num(b["n_acyclic"]), num(b["n_distinct_scaffolds"]),
                         "n/a" if b["unique_scaffold_ratio"] is None else f"{b['unique_scaffold_ratio']:.3f}", num(b["n50"]),
                         "n/a" if b["n50_subsample_mean"] is None else f"{b['n50_subsample_mean']:.0f}",
                         "n/a" if b["csr_auc"] is None else f"{b['csr_auc']:.3f}", "n/a" if b.get("reference_affinity") is None else f"{b['reference_affinity']:.3f}"])
        parts.append(table(["set", "molecules", "acyclic", "scaffolds", "unique ratio", "N50", f"N50 ({p4['routed']['subsample_n']:,}-subsample)", "CSR-AUC", "reference affinity"], rows))
        parts.append('<p class="note">Generic Bemis–Murcko scaffolds; acyclic molecules are excluded from scaffold statistics. A small N50 means collapse onto few cores. '
                     "The subsample N50 is only reported for sets at least as large as the subsample. Coverage on the routed set reflects which molecules were findable.</p>")

    res = rep.get("resources", {})
    if res:
        parts.append("<h2>Resources used</h2>")
        parts.append(table(["resource", "state", "snapshot", "size", "sha256"],
                           [[k, v.get("state"), v.get("snapshot_date") or "", num(v.get("size_bytes")) if v.get("size_bytes") else "", (v.get("sha256") or v.get("listing_sha256") or "")[:12]] for k, v in res.items()],
                           numeric_from=99))
    parts.append('<h2>How to read this report</h2><ul>'
                 "<li>Retrieved precedent evidence is not experimental feasibility; no qualifying precedent means support was not retrieved.</li>"
                 "<li>0.6412 is a similarity cutoff chosen on a calibration set, not a probability. The reaction database is shared with the planner's training data.</li>"
                 "<li>Results are comparable only within one stock and one set of planner settings.</li>"
                 "<li>Catalogue membership is a vendor-listed snapshot, not a live quote.</li></ul>")
    html = f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{E(title)}</title><style>{CSS}</style></head><body><main>{"".join(parts)}</main></body></html>'
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return out_path
