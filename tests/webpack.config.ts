import path from 'path';
import glob from 'glob';
import * as webpack from 'webpack';

import { BASE_CONFIG } from '../webpack.config';

// create one webpack entry for each test
const entry: Record<string, string> = {};
for (const file of glob.sync(path.join(__dirname, '**', '*.test.ts'))) {
    const name = path.basename(file).split('.')[0] + '.test';
    entry[name] = file;
}

BASE_CONFIG.module?.rules?.push({
    test: /\.node?$/,
    loader: 'node-loader',
    options: {
        name: '[path][name].[ext]',
    },
});

const config: webpack.Configuration = {
    ...BASE_CONFIG,
    mode: 'development',
    target: 'node',
    entry: entry,
    output: {
        filename: '[name].js',
        path: path.resolve(__dirname, '..', 'dist', 'tests'),
    },
    externals: {
        canvas: 'commonjs canvas',
    },
};

export default config;
