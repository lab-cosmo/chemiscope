import path from 'path';
import webpack from 'webpack';

import { WEBPACK_CONFIG } from '../webpack.config';

const config: webpack.Configuration = {
    ...WEBPACK_CONFIG,
    target: 'web',
    entry: {
        chemiscope: path.join(__dirname, 'jupyter', 'src', 'nbextension.ts'),
    },
    output: {
        filename: '[name].min.js',
        // Use AMD modules for the jupyter extension
        libraryTarget: 'amd',
        path: path.join(__dirname, 'jupyter', 'nbextension'),
        publicPath: '',
    },
    externals: ['@jupyter-widgets/base'],
};

export default config;
