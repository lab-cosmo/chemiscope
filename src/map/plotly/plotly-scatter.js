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
// plotly only accepts the symbols it knew when the module was loaded
scatter3d.attributes.marker.symbol.values = Object.keys(scatter3d.markerSymbols);

// below 1e5 points, plotly pads the axes range using the marker size of every
// point. This is quadratic when the size correlates with the position (e.g.
// when using the same property for size and x axis), taking tens of seconds for
// a few 10k points. Use the largest marker size as padding for all points
// instead (plotly itself uses the average size above 1e5 points).
const Axes = require('plotly.js/src/plots/cartesian/axes');
const findExtremes = Axes.findExtremes;
Axes.findExtremes = function (ax, data, opts) {
    if (opts && Array.isArray(opts.ppad) && opts.ppad.length > 1000) {
        opts = Object.assign({}, opts, { ppad: opts.ppad.reduce((a, b) => Math.max(a, b), 0) });
    }
    return findExtremes(ax, data, opts);
};

const Plotly = require('plotly.js/lib/core');

Plotly.register([require('plotly.js/lib/scattergl'), scatter3d]);

module.exports = Plotly;
