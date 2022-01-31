// These CSS imports return CSSStyleSheet objects for them to be manually
// "adopted" by shadow roots, unlike the default css-loader configuration
// which injects the CSS to the page's <head> element. The inline import syntax
// of Webpack is used here to avoid conflicts with the Jupyter Notebook build
// configuration.

import bootstrap from '!css-loader?exportType=css-style-sheet!bootstrap/dist/css/bootstrap.min.css';
import chemiscope from '!css-loader?exportType=css-style-sheet!./static/chemiscope.css';

export { bootstrap, chemiscope };
