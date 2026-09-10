selector_to_html = {"a[href=\"plotting.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pycmplot.plotting<a class=\"headerlink\" href=\"#pycmplot-plotting\" title=\"Link to this heading\">\uf0c1</a></h1><p>The plotting subpackage contains three modules: one for linear (stacked)\nManhattan plots, one for circular (Circos-style) Manhattan plots, and one\nfor QQ plots.</p>", "a[href=\"plotting.html#pycmplot-plotting-linear\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">pycmplot.plotting.linear<a class=\"headerlink\" href=\"#pycmplot-plotting-linear\" title=\"Link to this heading\">\uf0c1</a></h2><p>Generates single- and multi-track stacked linear Manhattan plots with\noptional significance lines, locus highlighting, cluster-aware label\nspreading, and intelligent arrow-angle calculation for gene annotations.</p>", "a[href=\"plotting.html#pycmplot-plotting-circular\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">pycmplot.plotting.circular<a class=\"headerlink\" href=\"#pycmplot-plotting-circular\" title=\"Link to this heading\">\uf0c1</a></h2><p>Generates multi-track Circos-style circular Manhattan plots. Track radii\nare computed automatically to give each track proportional visual weight\nrelative to its data range.</p>", "a[href=\"#api-reference\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">API Reference<a class=\"headerlink\" href=\"#api-reference\" title=\"Link to this heading\">\uf0c1</a></h1><p>This section documents every public function and class exposed by pycmplot.\nAll symbols listed here are importable directly from the top-level package\n(<code class=\"docutils literal notranslate\"><span class=\"pre\">from</span> <span class=\"pre\">pycmplot</span> <span class=\"pre\">import</span> <span class=\"pre\">...</span></code>) unless noted otherwise.</p>", "a[href=\"plotting.html#pycmplot-plotting-qq\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">pycmplot.plotting.qq<a class=\"headerlink\" href=\"#pycmplot-plotting-qq\" title=\"Link to this heading\">\uf0c1</a></h2><p>Produces QQ plots with 95 % beta-distribution confidence bands, optional\ngenome-wide significance lines, and genomic inflation (\u03bb) annotation.\nSupports log-uniform point thinning for fast plotting of large datasets.\nThree high-level layouts are provided: combined (grid of per-trait\npanels), separate (one file per trait), and overlay (all traits on one\nshared axes).</p>"}
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
