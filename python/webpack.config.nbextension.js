import path from 'path';

import { WEBPACK_CONFIG } from '../webpack.config.js';

const config = {
    ...WEBPACK_CONFIG,
    target: 'web',
    entry: {
        chemiscope: './python/jupyter/src/nbextension.ts',
    },
    output: {
        filename: '[name].min.js',
        // Use AMD modules for the jupyter extension
        libraryTarget: 'amd',
        path: path.resolve('./python/jupyter/nbextension'),
        publicPath: '',
    },
    externals: ['@jupyter-widgets/base'],
};

export default config;
