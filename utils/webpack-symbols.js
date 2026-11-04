/* eslint-disable */

/** Insert our symbol drawing helper into regl-scatter2d when webpack builds the app */
const path = require('path');

// match the function in the dependency compiled bundle
const DRAW_OPTIONS_HEADER =
    'Scatter.prototype.getMarkerDrawOptions = function (markerId, group, elements) {';

// quote the path so it can be inserted into a javascript require call
const SYMBOL_HELPER = JSON.stringify(path.resolve(__dirname, '../src/map/plotly/regl-symbols.ts'));

module.exports = function (rendererSource) {
    // stop if the insertion point is missing or ambiguous
    if (rendererSource.split(DRAW_OPTIONS_HEADER).length !== 2) {
        throw Error('regl-scatter2d changed, update utils/webpack-symbols.js');
    }

    // try our helper first, null lets the original function continue
    return rendererSource.replace(
        DRAW_OPTIONS_HEADER,
        () => `${DRAW_OPTIONS_HEADER}
  {
    const symbolDrawOptions = require(${SYMBOL_HELPER}).getSymbolDrawOptions(this, group, markerId, elements);
    if (symbolDrawOptions !== null) {
      return [symbolDrawOptions];
    }
  }`
    );
};
