/* eslint-disable */

import path from 'path';

import { WEBPACK_CONFIG } from '../webpack.config.mjs';

const config = {
    ...WEBPACK_CONFIG,
    entry: './python/streamlit/src/index.ts',
    output: {
        path: path.resolve('./python/streamlit/build'),
        filename: 'chemiscope-streamlit.min.js',
        publicPath: '',
    },
};

export default config;
