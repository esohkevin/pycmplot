selector_to_html = {"a[href=\"#id4\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.3.1 - 2026-07-31<a class=\"headerlink\" href=\"#id4\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Fixed</strong></p>", "a[href=\"#id8\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.2.5 - 2026-04-20<a class=\"headerlink\" href=\"#id8\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Fixed</strong></p>", "a[href=\"#id6\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.2.8 - 2026-05-30<a class=\"headerlink\" href=\"#id6\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Added</strong></p>", "a[href=\"#changelog\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Changelog<a class=\"headerlink\" href=\"#changelog\" title=\"Link to this heading\">\uf0c1</a></h1><p>All notable changes to <strong>pycmplot</strong> are documented here.</p><p>The format is based on <a class=\"reference external\" href=\"https://keepachangelog.com/en/1.0.0/\">Keep a Changelog</a>\nand this project adheres to <a class=\"reference external\" href=\"https://semver.org/\">Semantic Versioning</a>.</p>", "a[href=\"#id2\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.4.1 - 2026-08-23<a class=\"headerlink\" href=\"#id2\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Added</strong></p>", "a[href=\"#id7\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.2.7 - 2026-04-27<a class=\"headerlink\" href=\"#id7\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Added</strong></p>", "a[href=\"#id10\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.2.1 \u2014 2026-04-16<a class=\"headerlink\" href=\"#id10\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Added</strong></p>", "a[href=\"#id13\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.1.0 \u2014 2026-04-18<a class=\"headerlink\" href=\"#id13\" title=\"Link to this heading\">\uf0c1</a></h2><p>Initial release.</p><p><strong>Added</strong></p>", "a[href=\"#id11\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.1.9 \u2014 2026-04-14<a class=\"headerlink\" href=\"#id11\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Fixed</strong></p>", "a[href=\"#id9\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.2.2 - 2026-04-18<a class=\"headerlink\" href=\"#id9\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Added</strong></p><p>QQ plots (<a class=\"reference internal\" href=\"api/plotting.html#module-pycmplot.plotting.qq\" title=\"pycmplot.plotting.qq\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">pycmplot.plotting.qq</span></code></a>):</p>", "a[href=\"#id5\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.3.0 - 2026-06-01<a class=\"headerlink\" href=\"#id5\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Fixed</strong></p>", "a[href=\"#id12\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.1.8 \u2014 2026-04-14<a class=\"headerlink\" href=\"#id12\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Added</strong></p>", "a[href=\"#id3\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.4.0 - 2026-08-13<a class=\"headerlink\" href=\"#id3\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Changed</strong></p>", "a[href=\"#id1\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">0.4.2 - 2026-09-12<a class=\"headerlink\" href=\"#id1\" title=\"Link to this heading\">\uf0c1</a></h2><p><strong>Added</strong></p>"}
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
