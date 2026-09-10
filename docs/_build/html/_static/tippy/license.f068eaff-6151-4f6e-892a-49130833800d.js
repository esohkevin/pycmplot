selector_to_html = {"a[href=\"#license\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">License<a class=\"headerlink\" href=\"#license\" title=\"Link to this heading\">\uf0c1</a></h1><p>pycmplot is distributed under the terms of its open-source license.\nSee the <a class=\"reference external\" href=\"https://github.com/esohkevin/pycmplot/blob/main/LICENSE\">LICENSE</a>\nfile in the GitHub repository for the full text.</p>"}
skip_classes = ["headerlink", "sd-stretched-link"]

window.onload = function () {
    for (const [select, tip_html] of Object.entries(selector_to_html)) {
        const links = document.querySelectorAll(` ${select}`);
        for (const link of links) {
            if (skip_classes.some(c => link.classList.contains(c))) {
                continue;
            }

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
