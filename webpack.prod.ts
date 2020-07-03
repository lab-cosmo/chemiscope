import webpack from 'webpack';
import merge from 'webpack-merge';

import base from './webpack.base';

const config: webpack.Configuration = merge(base, {
    mode: 'production',
});

export default config;
