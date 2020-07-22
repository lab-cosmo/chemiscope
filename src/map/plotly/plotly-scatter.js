// A custom plotly bundle, based on
// - plotly.js/lib/index-gl2d
// - plotly.js/lib/index-gl3d
//
// Using this instead of the standard plotly build allow to shave ~2 Mib from
// the resulting minified JS. Without this, chemiscope.min.js is ~3.5 Mib, with
// this it is only 1.6 Mib.

/* eslint-disable */

'use strict';

const Plotly = require('plotly.js/lib/core');

Plotly.register([require('plotly.js/lib/scattergl'), require('plotly.js/lib/scatter3d')]);

module.exports = Plotly;
