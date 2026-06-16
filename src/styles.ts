// These CSS imports return CSSStyleSheet objects for them to be manually
// "adopted" by shadow roots, unlike the default css-loader configuration
// which injects the CSS to the page's <head> element. The inline import syntax
// of Webpack is used here to avoid conflicts with the Jupyter Notebook build
// configuration.

// `url=false` keeps `url(data:...)` icons as literal data URIs instead of
// turning them into `new URL(..., import.meta.url)` asset references. The latter
// breaks when the bundle is loaded as an ES module from a blob/data URL (e.g. by
// anywidget), where `import.meta.url` is not a valid base URL.
import bootstrap from '!css-loader?{"exportType":"css-style-sheet","url":false}!bootstrap/dist/css/bootstrap.min.css';
import chemiscope from '!css-loader?{"exportType":"css-style-sheet","url":false}!./static/chemiscope.css';

export { bootstrap, chemiscope };
