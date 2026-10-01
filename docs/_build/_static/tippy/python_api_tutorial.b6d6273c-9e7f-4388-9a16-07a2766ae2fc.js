selector_to_html = {"a[href=\"#python-api-tutorial\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Python API Tutorial<a class=\"headerlink\" href=\"#python-api-tutorial\" title=\"Link to this heading\">\uf0c1</a></h1><p>The interactive tutorial below demonstrates the full Python API using real\nsummary statistics data. It covers data loading, single- and multi-track\nlinear Manhattan plots, circular Manhattan plots, liftover, and the hits\nsummary table.</p>"}
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
