# Minimal Makefile for Sphinx docs
# Usage:
#   make html       – build HTML docs into _build/html/
#   make clean      – remove _build/
#   make livehtml   – auto-rebuild on file changes (requires sphinx-autobuild)

SPHINXOPTS  ?=
SPHINXBUILD ?= sphinx-build
SOURCEDIR    = .
BUILDDIR     = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help clean html livehtml

html:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@echo "\nBuild finished. Open _build/html/index.html in your browser."

clean:
	rm -rf $(BUILDDIR)

livehtml:
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
