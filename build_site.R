# Build the website: compile the Shiny app to WebAssembly, then add a proper
# loading screen. Run from this folder:   Rscript build_site.R
#
# Output goes to docs/ — that folder is what gets published to GitHub Pages.

shinylive::export("labviz", "docs")

# The stock page is a blank spinner for the ~minute it takes a first-time
# visitor to download R. That reads as "broken". Say what is happening.
loading <- '
<div id="lab-loading" style="position:fixed;inset:0;z-index:9999;background:#fff;
  display:flex;flex-direction:column;align-items:center;justify-content:center;gap:14px;
  font-family:-apple-system,BlinkMacSystemFont,\'Segoe UI\',Roboto,Helvetica,Arial,sans-serif;
  color:#1f2a37;text-align:center;padding:24px;">
  <div style="font-size:26px;font-weight:700;">Lab Figure Studio</div>
  <div style="font-size:17px;">Loading R into your browser&hellip;
    <span id="lab-loading-sec" style="color:#6b7280;font-variant-numeric:tabular-nums;"></span></div>
  <div style="width:280px;height:6px;border-radius:3px;background:#e5e7eb;overflow:hidden;">
    <div style="width:40%;height:100%;background:#0072B2;border-radius:3px;
      animation:labslide 1.4s ease-in-out infinite;"></div></div>
  <div style="font-size:14px;color:#6b7280;max-width:440px;line-height:1.5;">
    The first visit downloads R itself, so it takes up to a minute.
    After that it opens in seconds.<br>
    Your data never leaves this computer.</div>
  <div id="lab-loading-help" style="display:none;font-size:14px;color:#92400e;
    background:#fef3c7;border-radius:8px;padding:10px 14px;max-width:440px;line-height:1.5;">
    Still loading? If you opened this link inside another app (Slack, Teams,
    email, Claude), open it in Chrome, Safari or Firefox instead &mdash; those
    in-app browsers cannot run it.</div>
</div>
<style>@keyframes labslide{0%{transform:translateX(-110%)}100%{transform:translateX(260%)}}</style>
<script>
(function () {
  var t0 = Date.now(), sec = document.getElementById("lab-loading-sec");
  var iv = setInterval(function () {
    var s = Math.round((Date.now() - t0) / 1000);
    sec.textContent = s + "s";
    if (s >= 90) document.getElementById("lab-loading-help").style.display = "block";
    try {
      var f = document.querySelector("#root iframe"), w = f && f.contentWindow;
      if (w && w.Shiny && w.Shiny.shinyapp && w.Shiny.shinyapp.isConnected()) {
        document.getElementById("lab-loading").remove(); clearInterval(iv);
      }
    } catch (e) {}
  }, 500);
})();
</script>'

idx  <- "docs/index.html"
html <- readLines(idx, warn = FALSE)
html <- sub("<title>Shiny App</title>", "<title>Lab Figure Studio</title>", html, fixed = TRUE)
body <- grep('id="root"', html, fixed = TRUE)
html <- append(html, strsplit(loading, "\n")[[1]], after = body)
writeLines(html, idx)
file.create("docs/.nojekyll")   # GitHub Pages must not run Jekyll over this
cat("Built docs/ with loading screen.\n")
