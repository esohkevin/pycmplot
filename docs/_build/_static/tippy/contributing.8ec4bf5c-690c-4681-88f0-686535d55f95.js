selector_to_html = {"a[href=\"#code-style\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Code style<a class=\"headerlink\" href=\"#code-style\" title=\"Link to this heading\">\uf0c1</a></h2><p>pycmplot follows <a class=\"reference external\" href=\"https://peps.python.org/pep-0008/\">PEP 8</a>. Please run\n<code class=\"docutils literal notranslate\"><span class=\"pre\">flake8</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">black</span></code> before submitting a pull request:</p>", "a[href=\"#building-the-documentation-locally\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Building the documentation locally<a class=\"headerlink\" href=\"#building-the-documentation-locally\" title=\"Link to this heading\">\uf0c1</a></h2>", "a[href=\"#contributing\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Contributing<a class=\"headerlink\" href=\"#contributing\" title=\"Link to this heading\">\uf0c1</a></h1><p>Bug reports, feature requests, and pull requests are welcome on the\n<a class=\"reference external\" href=\"https://github.com/esohkevin/pycmplot\">GitHub repository</a>.</p>", "a[href=\"#reporting-issues\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Reporting issues<a class=\"headerlink\" href=\"#reporting-issues\" title=\"Link to this heading\">\uf0c1</a></h2><p>Please open a <a class=\"reference external\" href=\"https://github.com/esohkevin/pycmplot/issues\">GitHub issue</a>\nand include:</p>", "a[href=\"#docstrings\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Docstrings<a class=\"headerlink\" href=\"#docstrings\" title=\"Link to this heading\">\uf0c1</a></h2><p>All public functions should use\n<a class=\"reference external\" href=\"https://numpydoc.readthedocs.io/en/latest/format.html\">NumPy-style docstrings</a>.\nSee the existing modules for examples.</p>", "a[href=\"#development-setup\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Development setup<a class=\"headerlink\" href=\"#development-setup\" title=\"Link to this heading\">\uf0c1</a></h2>"}
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
