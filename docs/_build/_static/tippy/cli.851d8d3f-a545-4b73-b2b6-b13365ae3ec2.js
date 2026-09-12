selector_to_html = {"a[href=\"#column-auto-detection-options\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Column auto-detection options<a class=\"headerlink\" href=\"#column-auto-detection-options\" title=\"Link to this heading\">\uf0c1</a></h2><p>pycmplot infers column names automatically. Use these flags only if your\ncolumn names fall outside the recognised defaults.</p>", "a[href=\"#manhattan-plot-with-companion-qq-plots\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Manhattan plot with companion QQ plots<a class=\"headerlink\" href=\"#manhattan-plot-with-companion-qq-plots\" title=\"Link to this heading\">\uf0c1</a></h3>", "a[href=\"#qq-plot-options\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">QQ plot options<a class=\"headerlink\" href=\"#qq-plot-options\" title=\"Link to this heading\">\uf0c1</a></h2>", "a[href=\"#supplying-per-file-genome-builds\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Supplying per-file genome builds<a class=\"headerlink\" href=\"#supplying-per-file-genome-builds\" title=\"Link to this heading\">\uf0c1</a></h3><p>When summary statistics files do not carry a <code class=\"docutils literal notranslate\"><span class=\"pre\">BUILD</span></code> column, supply the\nbuilds in the same order as <code class=\"docutils literal notranslate\"><span class=\"pre\">--sum_stats</span></code>:</p>", "a[href=\"#plotting-behaviour-options\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Plotting behaviour options<a class=\"headerlink\" href=\"#plotting-behaviour-options\" title=\"Link to this heading\">\uf0c1</a></h2>", "a[href=\"#three-trait-stacked-linear-manhattan-plot\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Three-trait stacked linear Manhattan plot<a class=\"headerlink\" href=\"#three-trait-stacked-linear-manhattan-plot\" title=\"Link to this heading\">\uf0c1</a></h3>", "a[href=\"#circular-manhattan-plot\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Circular Manhattan plot<a class=\"headerlink\" href=\"#circular-manhattan-plot\" title=\"Link to this heading\">\uf0c1</a></h3>", "a[href=\"#input-output-options\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Input / output options<a class=\"headerlink\" href=\"#input-output-options\" title=\"Link to this heading\">\uf0c1</a></h2>", "a[href=\"#single-trait-linear-manhattan-plot\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Single-trait linear Manhattan plot<a class=\"headerlink\" href=\"#single-trait-linear-manhattan-plot\" title=\"Link to this heading\">\uf0c1</a></h3>", "a[href=\"#command-line-interface\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Command-Line Interface<a class=\"headerlink\" href=\"#command-line-interface\" title=\"Link to this heading\">\uf0c1</a></h1><p>pycmplot exposes a full command-line interface that mirrors the Python API.\nAfter installation, the <code class=\"docutils literal notranslate\"><span class=\"pre\">pycmplot</span></code> command is available in your PATH.</p>", "a[href=\"#example-commands\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Example commands<a class=\"headerlink\" href=\"#example-commands\" title=\"Link to this heading\">\uf0c1</a></h2><h3>Single-trait linear Manhattan plot<a class=\"headerlink\" href=\"#single-trait-linear-manhattan-plot\" title=\"Link to this heading\">\uf0c1</a></h3>"}
skip_classes = ["headerlink", "sd-stretched-link"]

window.onload = function () {
    for (const [select, tip_html] of Object.entries(selector_to_html)) {
        const links = document.querySelectorAll(` ${select}`);
        for (const link of links) {
            if (skip_classes.some(c => link.classList.contains(c))) {
                continue;
            }
            link.classList.add('has-tippy');
            tippy(link, {
                content: tip_html,
                allowHTML: true,
                arrow: true,
                placement: 'auto-start', maxWidth: 500, interactive: true, theme: 'light-border',

            });
        };
    };
    console.log("tippy tips loaded!");
};
