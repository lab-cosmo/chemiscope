// A custom plotly bundle, based on
// - plotly.js/lib/index-gl2d
// - plotly.js/lib/index-gl3d
//
// Using this instead of the standard plotly build allow to shave ~2 Mib from
// the resulting minified JS. Without this, chemiscope.min.js is ~3.5 Mib, with
// this it is only 1.6 Mib.

/* eslint-disable */

'use strict';

const markers3d = require('./markers3d.js');

// Require the module
const scatter3d = require('plotly.js/lib/scatter3d');
// monkey patch scatter3d to include more (and better!) symbols for 3d plots
// see https://github.com/plotly/plotly.js/issues/4205 in case this ever gets
// patched upstream and becomes unnecessary
for (const [k, v] of Object.entries(markers3d.default)) {
    scatter3d.markerSymbols[k] = v;
}

const Plotly = require('plotly.js/lib/core');

Plotly.register([require('plotly.js/lib/scattergl'), scatter3d]);

module.exports = Plotly;
