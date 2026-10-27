# Lab Figure Studio

A Shiny app that turns raw assay spreadsheets into publication figures, and
hands the user back the data and the code.

Three assays, one interface:

| Tab | Input | Output |
|---|---|---|
| **Limiting dilution** | cells seeded + organoids counted, *or* wells scored +/- per dose | frequency bar chart with CIs and significance brackets, dose–response / ELDA curve, frequency and comparison tables |
| **RT-qPCR** | long-format Ct table | ΔΔCt bar chart faceted by gene, Welch tests on ΔCt |
| **Volcano** | any DE results table | labelled volcano, hits table |

## Run it

```r
shiny::runApp("labviz")
```

Or open `app.R` in RStudio and click **Run App**. Missing packages install on
first launch.

## Batch mode

The button that actually gets you out of the loop. On the limiting-dilution
tab, set **Split figures by** to a column such as cell line, screen, or
timepoint. The app then writes one frequency chart and one curve *per level*
into a single ZIP:

```
Plots/    HT29_frequency.pdf/.png, HT29_curve.pdf/.png, LS174T_..., SW480_...
Tables/   per-level frequencies + comparisons, plus pooled ALL_*.csv
Inputs/   the data exactly as analysed
R/        the fitting and plotting functions
analysis.R, README.txt
```

Levels whose model will not fit are skipped and named in `README.txt` rather
than failing the whole export. Borrowed from Metabolite AutoPlotter
(Pietzke & Vazquez 2020), which does the same thing one-plot-per-metabolite.

## Protocol sheets

Real lab sheets are not tidy tables. The importer under **Protocol sheet** on
the limiting-dilution tab reads one as written: blocks separated by blank
rows, the cell line and media floating above each block as a title, no column
headers, and no consistent ordering.

It works out which column is which from the numbers themselves:

- the **organoid** column is the one that never exceeds any other count —
  you cannot form more organoids than you seeded cells
- the **protocol** column is round and repeating (100, 200, 300, 400)
- the **counted** column is noisy and near-unique (87, 147, 245, 337)

That last distinction is the point of the feature. A sheet that records both
what you meant to seed and what you actually counted will give frequencies
roughly 20% too low if the model is fitted to the protocol number. The
importer defaults to the counted column, says so, and puts a table on the
figure tab showing exactly what the other choice would do to every group.

Verified against `LDA 0514.xlsx`: two blocks, seven groups, all seven
frequencies reproduced to seven significant figures from the untouched
workbook.

## Layout

```
labviz/
├── app.R                  # navbar + wiring, nothing else
├── R/                     # auto-sourced by Shiny
│   ├── 00_theme.R         # ggplot theme, palettes, p-value stars
│   ├── 01_io.R            # upload / paste / column-mapping module
│   ├── 02_lda_stats.R     # the cloglog GLM
│   ├── 03_lda_plots.R     # bar + curve + ELDA figures
│   ├── 10_mod_lda.R       # ┐
│   ├── 20_mod_qpcr.R      # ├ one Shiny module per assay
│   ├── 30_mod_volcano.R   # ┘
│   └── 40_export.R        # PNG/PDF/SVG + the reproducible-code ZIP
└── www/custom.css
```

Adding a fourth assay means adding one `R/` file with a `*SidebarUI`, `*UI`
and `*Server`, then three lines in `app.R`. Nothing else changes.

The statistics and plotting files (`00`, `02`, `03`) are deliberately free of
any Shiny dependency — they are plain functions. That is what lets the export
button ship them to the user as a working standalone script.

## Deploying to GitHub Pages

The app is built to WebAssembly and served as static files, so there is no
server to run, nothing to pay for, and no R process anywhere that could
receive lab data.

Push this directory to a repo and let CI do the build — `.github/workflows/deploy.yml`
is already written. In the repo, **Settings → Pages → Source → GitHub Actions**,
then push to `main`. The action installs the packages, runs
`shinylive::export("labviz", "_site")`, and publishes. The build takes a few
minutes; the resulting site is about 106 MB, well inside the 1 GB Pages limit.

`docs/` and `_site/` are gitignored on purpose. Committing a 100 MB build
output bloats the repo and lets the deployed app drift from the source — let
the action rebuild it every push.

To build and check it locally before pushing:

```r
shinylive::export("labviz", "docs")
httpuv::runStaticServer("docs", port = 8008)
```

Note that a service worker is required, so `file://` will not work — you need
a real local server, and some embedded browsers block service workers
entirely. Chrome, Firefox and Safari are all fine.

### Verified WebAssembly compatibility

Every package the app uses is available for webR, and three things were
changed so the app works in the browser:

| Problem in webR | Fix |
|---|---|
| `utils::zip()` shells out to a `zip` binary that does not exist | `zip::zip()` — pure C, no shell |
| `png()` wants cairo/X11 | `ragg::agg_png()` |
| `install.packages()` at startup; `font_google()` fetches from Google | skipped under Emscripten; system font stack |

The exported bundle is genuinely self-contained. `library.data.gz` (15 MB)
carries shiny, bslib, htmltools, sass and 32 more; `shinylive/webr/packages/`
carries the other 39, including ggplot2, DT, svglite, colourpicker and ragg.
Nothing is fetched from `repo.r-wasm.org` or any other third party at runtime,
and after the swap above there are no requests to Google Fonts either — the
page talks only to its own origin.

## Security

The privacy property is architectural, not a promise: uploaded data is read
into the browser tab's sandboxed memory and analysed by an R interpreter
running in that same tab. There is no server-side session to leak, because
there is no server-side session.

Anyone can verify it in under a minute — DevTools → Network shows requests
only while the page loads, none afterwards; or load the page, disconnect from
the network, and keep working.

Two real limits, stated plainly:

- **GitHub Pages is public.** On Free, Pro and Team accounts the published
  site is world-readable even when the source repo is private; private Pages
  needs GitHub Enterprise Cloud (Settings → Pages → Visibility → Private).
  This is fine here — the *app* being public costs nothing, because the *data*
  never goes there. But do not treat an unguessable URL as a control.
- **Never commit real data.** The `.gitignore` blocks `*.csv`, `*.tsv`,
  `*.xlsx`, `*.rds` and `.RData` for that reason. The example datasets are
  synthetic and defined inline in `R/`, so there is no data file to whitelist.

## Notes on the statistics

Limiting dilution uses a binomial GLM with a complementary log-log link.
Both supported layouts are the same model; the well-scoring layout is the ELDA
method of Hu & Smyth (2009, *J Immunol Methods* 347:70) expressed as an
`offset(log(dose))`. Frequencies are `exp(β)`, intervals are Wald on the log
scale, group comparisons are Wald tests on `β_A − β_B`.

The residual-deviance goodness-of-fit p-value is shown under the comparison
table. Treat a small value as a warning that one frequency per group does not
describe the data — usually non-linear dose scaling or disagreeing replicates.

qPCR statistics run on ΔCt values, never on fold changes.
