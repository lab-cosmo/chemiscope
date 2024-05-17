from sphinx_gallery.scrapers import matplotlib_scraper

class SphinxGalleryScrapper:
    def __repr__(self):
        return 'SphinxGalleryScrapper'

    def __call__(self, block, block_vars, gallery_conf):
        variables = block_vars.get("example_globals")
        widget = variables.get("___")
        if widget:
            widget._repr_html_ = lambda: self.generate_html_content(widget.value)
        return matplotlib_scraper(block, block_vars, gallery_conf)

    def generate_html_content(self, dataset):
        div_id = "sphinx-gallery-" + str(id(self))
        html_content = f"""
            <!-- Load all dependencies -->
            <!-- jquery -->
            <script
                src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js"
                integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"
            ></script>

            <!-- Chemiscope code and default viewer code -->
            <script src="../../../_static/js/chemiscope.min.js"></script>

            <!-- HTML Viewer code -->
            <script async="async" src="../../../_static/js/chemischope-sphinx-gallery.js"></script>

            <!-- font-awesome -->
            <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.11.2/css/all.min.css" />
            <link rel="icon" type="image/png" href="../../../_static/chemiscope-icon.png" sizes="32x32" />

            <!-- CSS styles -->
            <link rel="stylesheet" href="../../../_static/chemischope-sphinx-gallery.css" type="text/css" />
            
            <div id="{div_id}-warning-display" style="display: none; flex-wrap: wrap; border: 1px solid #f0ad4e; background-color: #fcf8e3; color: #8a6d3b; padding: 15px; border-radius: 4px; margin-bottom: 20px;">
                <p style="flex: 1"></p>
                <button
                    type="button"
                    style="position: relative; float: right; font-size: 3em; line-height: 1; color: inherit; text-shadow: none; background-color: transparent; border: 0; -webkit-appearance: none;"
                    aria-label="Close"
                    onclick=hideElement("{div_id}-warning-display")
                >&times;</button>
            </div>
            <div style="padding: 0.5rem 1.25rem">&nbsp;</div>

            <div id="{div_id}-loading">
                <i class="fa fa-spinner fa-spin"></i>
            </div>
            <div id="{div_id}"></div>
            <script>
                loadChemiscopeSphinxGallery("{div_id}", {dataset});
            </script>
        """
        return html_content
