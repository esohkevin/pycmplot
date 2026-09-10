selector_to_html = {"a[href=\"#module-pycmplot.constants\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">pycmplot.constants<a class=\"headerlink\" href=\"#module-pycmplot.constants\" title=\"Link to this heading\">\uf0c1</a></h2><p>Genome-level constants shared across pycmplot modules.</p>", "a[href=\"#pycmplot-constants\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">pycmplot.constants<a class=\"headerlink\" href=\"#pycmplot-constants\" title=\"Link to this heading\">\uf0c1</a></h1><p>Genome-level constants shared across pycmplot modules: GRCh38 chromosome\nlengths, gene biotype priority weights, and standard chromosome ordering.</p>", "a[href=\"#contents\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Contents<a class=\"headerlink\" href=\"#contents\" title=\"Link to this heading\">\uf0c1</a></h3><p class=\"rubric\">Notes</p><p><code class=\"docutils literal notranslate\"><span class=\"pre\">hg38_chr_lengths</span></code> reflects the GRCh38 primary assembly (GCA_000001405).\nValues may differ slightly from builds that include alternate contigs or\npatches.</p>"}
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
