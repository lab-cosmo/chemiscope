import path from 'path';

import { WEBPACK_CONFIG } from '../webpack.config.mjs';

const mode = process.env.CHEMISCOPE_WEBPACK_MODE || 'production';

// Build a single ES module loaded by anywidget through the `_esm` trait. This
// replaces the previous separate JupyterLab labextension and classic-notebook
// nbextension bundles.
const config = {
    ...WEBPACK_CONFIG,
    mode: mode,
    target: 'web',
    entry: {
        'chemiscope-widget': './python/widget/src/anywidget.ts',
    },
    experiments: {
        outputModule: true,
    },
    output: {
        filename: '[name].mjs',
        library: {
            type: 'module',
        },
        path: path.resolve('./python/chemiscope/static'),
        publicPath: '',
    },
    // we know we bundle large dependencies (Plotly.js, 3Dmol.js), we don't need
    // webpack to warn us about the bundle size every time.
    performance: {
        hints: false,
    },
};

export default config;
