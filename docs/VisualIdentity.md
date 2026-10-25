# LIMA visual identity

Rules for anything user facing: the website, reports, figures, the app. Derived from the website (`LIMA_web/public/site.css`) and
the logo. Change a rule here first, then apply it.

## Character
- Calm, technical, precise. "Serious physics. Modern workflows." No gradients, glow, drop shadows, stock imagery or emoji.
- Dark only. Lots of empty space; hierarchy comes from size, weight and the single accent, not from boxes and colour.
- Numbers are the hero. State them plainly, with units, and let the layout make them stand out.

## Colour
| Token | Hex | Use |
|---|---|---|
| `--paper` | `#14171e` | Page background |
| `--tile` | `#1c212b` | Raised surfaces: panels, windows |
| `--line` | `#303743` | Hairlines, borders, gridlines, axes |
| `--ink` | `#e4e7ef` | Primary text |
| `--soft` | `#9ca5b7` | Secondary text, labels, captions |
| `--accent` | `#a5dfc5` | LIMA itself: its data, links, the one thing to look at |
| `--button-ink` | `#14281f` | Text on a filled accent surface |

- One accent. Mint means LIMA; nothing else is mint.
- Dark is the default (website, app). Light is for social media and print, via `[data-theme="light"]` overriding the tokens:

| Token | Light |
|---|---|
| `--paper` | `#f6f7f9` |
| `--tile` | `#ffffff` |
| `--line` | `#dde1e8` |
| `--ink` | `#14171e` |
| `--soft` | `#5b6475` |
| `--accent` (fills) | `#6cc39c` |
| `--accent-text` (accent text) | `#1f7a55` |
| comparison neutral | `#b3b9c5` |

- On dark, `--accent-text` equals `--accent`. Use `--accent-text` for mint text and `--accent` for mint shapes.
- Tints of the accent for fills: `rgba(165,223,197,.09)` background, `.3` border, `.16`/`.55` on hover.
- Comparison data (other engines, references) uses the neutrals: `--soft` and dimmer greys, never a second hue.

## Type
- Body and headings: Arial, Helvetica, sans-serif. Headings weight 500 with tight negative tracking
  (h1 about -2.8px at 60px, h2 about -0.5px at 22px).
- Labels, metadata, units, eyebrows, buttons, code: Consolas, 'Courier New', monospace, 11-13px.
- Eyebrows: monospace, 11px, uppercase, letter-spacing 1.6px, in `--accent`.
- Body text 15px, line-height 1.75. Secondary prose in `--soft`.
- Numbers in tables and charts: monospace, tabular, so digits align.

## Shape
- Hairline borders, 1px `--line`. Radius 5px for controls, 8-9px for panels.
- Layout: max width 1160px; generous padding (38-48px desktop, 22px mobile).

## Logo
- White line-art mark (`LIMA_web/public/resources/limalogo/iconwhite.svg`) with the wordmark "LIMA" in Arial bold,
  -1px tracking, followed by a monospace `/ MD` in `--soft`. On light backgrounds the mark takes `--ink` (`#14171e`); otherwise
  never recoloured or outlined.

## Charts and figures
- Background `--paper` (or `--tile` inside a panel), gridlines and axes in `--line`, labels in `--soft` monospace.
- LIMA's series in `--accent`; every other series in neutrals. Order: LIMA last or on top, so it reads as the result.
- Direct labels on bars and lines instead of legends where possible. No 3D, no chart junk.
- Every number has its unit, and every comparison states the hardware and settings it was measured with.
- Comparisons: one thin line (6px, rounded) per engine, stacked, LIMA first. Each line is named on its left and valued on its
  right in 12px mono, so it reads without a legend: LIMA in accent, others in `--soft` with a neutral line (`#5d6576` dark,
  `#b3b9c5` light). Published or reference values: a dashed neutral outline, with its source under the name.
  Scale each item's lines to its own longest value when items differ by orders of magnitude; the numbers carry the absolutes.
  (Tried and rejected: overlapping fat bars with values inside, which needed a legend and read poorly.)
- Large numbers with thin-space thousands separators (`2 403`); small ones as `6.7×10⁻⁴`, not `6.7e-4`.
- Ratios as the hero number: `2.5×`, accent, weight 500, the `×` smaller.

## Molecule images
- Rendered by LIMA itself (`lima render --whole --hide-water`), cropped to the molecules on a transparent background, so they sit
  on `--paper`. Keep the renderer's element colours; don't recolour to the palette.
- Each image fills its frame; never scale by system size. The apparent atom size already tells the viewer the scale.

## Reports
- One page: header (logo, dateline in mono), accent eyebrow, a one-line claim as h1, a one-sentence lede, up to three hero
  stats, the data rows, and a three-column method footer in 11px mono. Self-contained HTML (images inline), dark, printable on
  A4 landscape.
- Say only what the reader needs; no obvious statements. Per row: name, a few words of description, size (atom count). Details
  such as force fields or setup tricks go in the method footer.
- External reports show user-facing units only (ns/day). Internal metrics such as µs/step stay in internal reports.
