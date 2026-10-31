"""The one-page LIMA vs GROMACS report: report.html of `lima gromacscompare`, and of LIMA_data/salad_gromacs_bench/report/build_report.py,
which also makes the website page and social cards from it. Standard library only. Follows LIMA/docs/VisualIdentity.md.

A row is a dict per system: system, atoms, and whichever of these were measured (missing ones show as a dash):
lima, gromacs (ns/day), forces (relative RMS error), rmsdLima, rmsdGmx1, rmsdGmx2 (nm), rmsfR, rmsfRGmx, plus
structureMark (a footnote mark after the RMSD and RMSF r labels) and error (why the system was not completed)."""
import base64
import html
import math

SYSTEMS = {
    # id: display name, description, force field, box [nm], note
    "villin": ("Villin headpiece", "Small protein in water", "AMBER99SB", (4.30, 3.29, 3.52), None),
    "rnase_cubic": ("RNase A", "Protein in water", "AMBER99SB", (6.28, 6.28, 6.28), None),
    "adh_dodec": ("Alcohol dehydrogenase", "Protein tetramer in water", "AMBER99SB-ILDN", (11.03, 11.03, 15.60), "2 copies"),
    "aqp_ensemble": ("Aquaporin", "Membrane channels in a POPC bilayer", "CHARMM36", (11.62, 20.13, 9.13), "2 copies"),
    "ion_channel": ("GluCl ion channel", "Membrane protein in a DOPC bilayer", "AMBER99SB-ILDN + Berger", (11.25, 19.49, 14.83), "2 copies"),
    "stmv": ("Satellite tobacco mosaic virus", "Capsid and RNA in water", "CHARMM", (21.68, 21.68, 21.68), None),
}
DASH = "–"


def data_uri(path, mime):
    return f"data:{mime};base64,{base64.b64encode(path.read_bytes()).decode()}"


def thousands(value):
    return f"{value:,.0f}".replace(",", " ")


def nsday(value):
    return thousands(value) if value >= 100 else f"{value:.1f}"


def scientific(value):
    exponent = math.floor(math.log10(value))
    return f"{value / 10 ** exponent:.1f}<span class=times>×</span>10<sup>{exponent}</sup>"


def speedups(rows):
    return [r["lima"] / r["gromacs"] for r in rows if r.get("lima") and r.get("gromacs")]


def geometric_mean(values):
    return math.exp(sum(math.log(x) for x in values) / len(values))


def headline(rows):
    s = speedups(rows)
    return f"LIMA runs up to {max(s):.1f}× faster than GROMACS" if s and max(s) > 1 else "LIMA vs GROMACS"


def stats(rows):
    s = speedups(rows)
    forces = [r["forces"] for r in rows if r.get("forces")]
    tiles = []
    if s:
        tiles += [f'<div class="stat"><b>{geometric_mean(s):.1f}×</b><span>Mean speedup</span></div>',
                  f'<div class="stat"><b>{max(s):.1f}×</b><span>Largest speedup</span></div>']
    if forces:
        tiles.append(f'<div class="stat"><b>≤ {scientific(max(forces))}</b><span>Force error vs GROMACS</span></div>')
    return '<div class="stats">\n      ' + '\n      '.join(tiles) + '\n    </div>'


def row(r, src):
    """One system's row; src is its thumbnail's URL"""
    name, description, *_ = SYSTEMS[r["system"]]
    lima, gromacs = r.get("lima"), r.get("gromacs")
    if lima and gromacs:
        longest = max(lima, gromacs)
        # A thin line per engine, named on the left and valued on the right, scaled to the row's longest
        line = lambda kind, label, value: (f'<span class="lab {kind}">{label}</span><div class="track"><div class="fill {kind}" '
                                           f'style="width:{100 * value / longest:.2f}%"></div></div><span class="val {kind}">{nsday(value)}</span>')
        bars = line("lima", "LIMA", lima) + line("gromacs", "GROMACS", gromacs)
        speedup = f"{lima / gromacs:.1f}<span>×</span>"
    else:
        bars = f'<p class="error">Not completed: {html.escape(r.get("error", "no performance result").splitlines()[0])}</p>'
        speedup = DASH
    mark = r.get("structureMark", "")
    forces = scientific(r["forces"]) if r.get("forces") else DASH
    if r.get("rmsdLima") is not None:
        low, high = min(r["rmsdGmx1"], r["rmsdGmx2"]), max(r["rmsdGmx1"], r["rmsdGmx2"])
        rmsd = f"{r['rmsdLima']:.2f} nm<small>({low:.2f}–{high:.2f})</small>"
        rmsf = f"{r['rmsfR']:.2f}<small>({r['rmsfRGmx']:.2f})</small>"
    else:
        rmsd = rmsf = DASH
    return f"""
      <article class="row">
        <div class="thumb"><img src="{src}" alt="{html.escape(name)}"></div>
        <div class="system">
          <h2>{html.escape(name)}</h2>
          <p>{html.escape(description)}</p>
          <p class="meta">{thousands(r['atoms'])} atoms</p>
        </div>
        <div class="bars">{bars}</div>
        <div class="speedup">{speedup}</div>
        <div class="parity">
          <span>Forces</span><b>{forces}</b>
          <span>RMSD{mark}</span><b>{rmsd}</b>
          <span>RMSF r{mark}</span><b>{rmsf}</b>
        </div>
      </article>"""


def table(rows_html):
    return f"""<div class="table">
    <div class="columns"><span></span><span>System</span><span>Throughput, ns/day</span><span class="r">Speedup</span><span>Parity · LIMA (GROMACS)</span></div>
    {''.join(rows_html)}
  </div>"""


def page(theme, layout, extra_css, logo_svg, dateline, headline, lede, stats, table, social_note, method):
    """The self-contained page. theme is dark or light; layout is the root element's class, and extra_css is appended to the
    stylesheet (eg. the social cards' canvases)"""
    # The logo is white line art; on light it takes the ink colour
    if theme == "light":
        logo_svg = logo_svg.replace("#FDF6E3", "#14171e")
    logo = f"data:image/svg+xml;base64,{base64.b64encode(logo_svg.encode()).decode()}"
    return f"""<!DOCTYPE html>
<html lang="en" data-theme="{theme}" class="{layout}">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<meta name="color-scheme" content="{theme}">
<title>LIMA vs GROMACS</title>
<style>
:root {{ color-scheme: dark; --ink: #e4e7ef; --paper: #14171e; --soft: #9ca5b7; --line: #303743; --accent: #a5dfc5; --accent-text: #a5dfc5;
  --tile: #1c212b; --dim: #5d6576; --mono: Consolas, 'Courier New', monospace; font-family: Arial, Helvetica, sans-serif; color: var(--ink); background: var(--paper); }}
:root[data-theme="light"] {{ color-scheme: light; --ink: #14171e; --paper: #f6f7f9; --soft: #5b6475; --line: #dde1e8; --accent: #6cc39c;
  --accent-text: #1f7a55; --tile: #ffffff; --dim: #b3b9c5; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--paper); }}
.page {{ max-width: 1160px; margin: auto; padding: 38px 48px 30px; }}
header {{ display: flex; justify-content: space-between; align-items: center; gap: 24px; flex-wrap: wrap; }}
.brand {{ display: inline-flex; align-items: center; gap: 8px; font-size: 25px; font-weight: 700; letter-spacing: -1px; }}
.brand img {{ width: 40px; height: 40px; margin-right: 4px; }}
.brand span, .dateline {{ font: 12px var(--mono); letter-spacing: 0; color: var(--soft); }}
h1 {{ font-size: clamp(32px, 3.6vw, 44px); font-weight: 500; line-height: 1.1; letter-spacing: -2px; margin: 48px 0 12px; }}
.lede {{ color: var(--soft); font-size: 15px; line-height: 1.75; max-width: 680px; margin: 0; }}
.stats {{ display: flex; gap: 56px; flex-wrap: wrap; margin: 28px 0 30px; }}
.stat b {{ display: block; white-space: nowrap; font: 500 40px/1 Arial, Helvetica, sans-serif; letter-spacing: -1.5px; color: var(--accent-text); }}
.stat b sup {{ font-size: 0.5em; letter-spacing: 0; }}
.stat > span {{ display: block; margin-top: 10px; font: 11px var(--mono); color: var(--soft); text-transform: uppercase; letter-spacing: 1.1px; }}
.columns, .row {{ display: grid; grid-template-columns: 96px minmax(160px, 1fr) minmax(260px, 2.4fr) 64px 200px; gap: 26px; align-items: center; }}
.columns {{ font: 11px var(--mono); text-transform: uppercase; letter-spacing: 1.1px; color: var(--soft); padding-bottom: 12px; border-bottom: 1px solid var(--line); }}
.columns .r, .speedup {{ text-align: right; }}
.row {{ padding: 10px 0; border-bottom: 1px solid var(--line); }}
.thumb {{ height: 88px; }}
.thumb img {{ width: 100%; height: 100%; object-fit: contain; }}
.system h2 {{ font-size: 17px; font-weight: 500; letter-spacing: -0.3px; margin: 0 0 4px; }}
.system p {{ margin: 0; color: var(--soft); font-size: 13px; line-height: 1.5; }}
.system .meta {{ font: 11px/1.6 var(--mono); margin-top: 4px; }}
.bars {{ display: grid; grid-template-columns: 74px 1fr 64px; gap: 8px 14px; align-items: center; font: 12px var(--mono); }}
.lab {{ color: var(--soft); white-space: nowrap; }}
.lab small {{ display: block; font-size: 10px; }}
.val {{ text-align: right; color: var(--soft); }}
.lab.lima, .val.lima {{ color: var(--accent-text); }}
.track {{ height: 6px; }}
.fill {{ height: 100%; border-radius: 3px; }}
.fill.lima {{ background: var(--accent); }}
.fill.gromacs {{ background: var(--dim); }}

.speedup {{ font: 500 26px/1 Arial, Helvetica, sans-serif; letter-spacing: -1px; color: var(--accent-text); }}
.speedup span {{ font-size: 16px; margin-left: 1px; }}
.parity {{ display: grid; grid-template-columns: 52px 1fr; gap: 3px 10px; font: 12px var(--mono); text-align: left; }}
.parity span {{ color: var(--soft); }}
.parity b {{ font-weight: 400; color: var(--ink); white-space: nowrap; }}
.parity small {{ color: var(--soft); margin-left: 8px; }}
.parity sup {{ font-size: 9px; }}
.times {{ margin: 0 2px; color: var(--soft); }}
footer {{ display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 32px; margin-top: 26px; font: 11px/1.75 var(--mono); color: var(--soft); }}
footer h3 {{ font: 11px var(--mono); text-transform: uppercase; letter-spacing: 1.1px; color: var(--ink); margin: 0 0 8px; }}
footer p {{ margin: 0 0 8px; }}
footer b {{ font-weight: 400; color: var(--ink); }}
@media (max-width: 900px) {{
  .page {{ padding: 26px 22px; }}
  .columns {{ display: none; }}
  .row {{ grid-template-columns: 72px 1fr auto; grid-template-areas: "thumb system speedup" "bars bars bars" "parity parity parity"; gap: 12px 16px; }}
  .thumb {{ grid-area: thumb; height: 72px; }} .system {{ grid-area: system; }} .bars {{ grid-area: bars; }} .speedup {{ grid-area: speedup; }}
  .parity {{ grid-area: parity; }}
  footer {{ grid-template-columns: 1fr; }}
}}
@media print {{ @page {{ size: A4 landscape; margin: 10mm; }} body {{ -webkit-print-color-adjust: exact; print-color-adjust: exact; }} .page {{ padding: 0; }} .eyebrow {{ margin-top: 28px; }} }}
.social-note {{ display: none; }}
{extra_css}
</style>
</head>
<body>
<div class="page">
  <header>
    <div class="brand"><img src="{logo}" alt="">LIMA <span>/ MD</span></div>
    <div class="dateline">{dateline}</div>
  </header>

  <div class="lead">
    <h1>{headline}</h1>
    <p class="lede">{lede}</p>

    {stats}
  </div>

  {table}

  <p class="social-note">{social_note}</p>

  <footer>{method}</footer>
</div>
</body>
</html>
"""
