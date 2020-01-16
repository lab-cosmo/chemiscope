// A custom plotly bundle, based on
// - plotly.js/lib/index-gl2d
// - plotly.js/lib/index-gl3d
//
// Using this instead of the standard plotly build allow to shave 1.5 Mib from
// the resulting minified JS.

'use strict';

var Plotly = require('plotly.js/lib/core');

Plotly.register([
    require('plotly.js/lib/scattergl'),
    require('plotly.js/lib/scatter3d'),
]);

module.exports = Plotly;
