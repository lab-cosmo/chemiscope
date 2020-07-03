import webpack from 'webpack';
import merge from 'webpack-merge';

import base from './webpack.base';

const config: webpack.Configuration = merge(base, {
    devtool: 'inline-source-map',
    mode: 'development',
});

export default config;
