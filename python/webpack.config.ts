import path from 'path';
import webpack from 'webpack';

import { WEBPACK_CONFIG } from '../webpack.config';

const config: webpack.Configuration = {
    ...WEBPACK_CONFIG,
    target: 'web',
    entry: {
        'chemiscope-widget': path.join(__dirname, 'nbextension', 'src', 'index.ts'),
    },
    output: {
        filename: '[name].min.js',
        // Use AMD modules for the jupyter extension
        libraryTarget: 'amd',
        path: path.join(__dirname, 'nbextension', 'build'),
        publicPath: '',
    },
    externals: ['@jupyter-widgets/base'],
};

export default config;
