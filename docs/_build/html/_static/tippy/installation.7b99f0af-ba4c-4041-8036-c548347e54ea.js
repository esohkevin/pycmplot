selector_to_html = {"a[href=\"#verify-the-installation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Verify the installation<a class=\"headerlink\" href=\"#verify-the-installation\" title=\"Link to this heading\">\uf0c1</a></h2><p>A successful installation prints the full help message with all command-line\noptions.</p>", "a[href=\"#using-a-virtual-environment\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Using a virtual environment<a class=\"headerlink\" href=\"#using-a-virtual-environment\" title=\"Link to this heading\">\uf0c1</a></h2><p>If you cannot or prefer not to install into the system Python, use a virtual\nenvironment:</p>", "a[href=\"#from-pypi-recommended\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">From PyPI (recommended)<a class=\"headerlink\" href=\"#from-pypi-recommended\" title=\"Link to this heading\">\uf0c1</a></h2>", "a[href=\"#requirements\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Requirements<a class=\"headerlink\" href=\"#requirements\" title=\"Link to this heading\">\uf0c1</a></h2><p>pycmplot requires <strong>Python 3.9 or later</strong>. The following packages are installed\nautomatically as dependencies:</p>", "a[href=\"#from-github-development-version\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">From GitHub (development version)<a class=\"headerlink\" href=\"#from-github-development-version\" title=\"Link to this heading\">\uf0c1</a></h2><p>If your system Python is externally managed (common on Ubuntu 22.04+):</p>", "a[href=\"#installation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Installation<a class=\"headerlink\" href=\"#installation\" title=\"Link to this heading\">\uf0c1</a></h1><h2>Requirements<a class=\"headerlink\" href=\"#requirements\" title=\"Link to this heading\">\uf0c1</a></h2><p>pycmplot requires <strong>Python 3.9 or later</strong>. The following packages are installed\nautomatically as dependencies:</p>"}
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
