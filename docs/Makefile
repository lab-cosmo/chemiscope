# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?= -W
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = src
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

TYPESCRIPT_SOURCES=$(wildcard ../src/*.ts ../src/*/*.ts ../src/*/*/*.ts)
$(BUILDDIR)/api-docs: $(TYPESCRIPT_SOURCES)
	npm run api-docs
	touch $(BUILDDIR)/api-docs

.PHONY: help Makefile

# Declare a rule for .ts file, so that the rule below is not used
%.ts:
	exit 1

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: $(BUILDDIR)/api-docs Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
