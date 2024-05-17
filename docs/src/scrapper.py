from sphinx_gallery.scrapers import matplotlib_scraper

class SphinxGalleryScrapper:
    def __repr__(self):
        return 'SphinxGalleryScrapper'

    def __call__(self, block, block_vars, gallery_conf):
        variables = block_vars.get("example_globals")
        widget = variables.get("___")
        if widget:
            div_id = "sphinx-gallery-" + str(id(self))            
            widget._repr_html_ = lambda: self.generate_html_content(div_id, widget.value)
        return matplotlib_scraper(block, block_vars, gallery_conf)

    def generate_html_content(self, div_id, dataset):
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
            
            <div id="{div_id}-loading">
                <i class="fa fa-spinner fa-spin"></i>
            </div>
            <div id="{div_id}"></div>
            <script>
                loadChemiscopeSphinxGallery("{div_id}", {dataset});
            </script>
        """
        return html_content
