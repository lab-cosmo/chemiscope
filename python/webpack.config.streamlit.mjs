/* eslint-disable */

import path from 'path';

import { WEBPACK_CONFIG } from '../webpack.config.mjs';

const mode = process.env.CHEMISCOPE_WEBPACK_MODE || 'production';

const config = {
    ...WEBPACK_CONFIG,
    mode: mode,
    entry: './python/streamlit/src/index.ts',
    output: {
        path: path.resolve('./python/streamlit/build'),
        filename: 'chemiscope-streamlit.min.js',
        publicPath: '',
    },
};

export default config;
