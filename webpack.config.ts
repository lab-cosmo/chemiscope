import path from "path";
import webpack from "webpack";

const config: webpack.Configuration = {
    mode: "development",
    entry: {
        "chemiscope": "./src/index.ts",
        "jsmol-widget": "./src/structure/widget.ts",
    },
    output: {
        path: path.resolve(__dirname, "dist"),
        filename: "[name].min.js",
        libraryTarget: "var",
        library: "Chemiscope",
    },
    resolve: {
        modules: ['./node_modules'],
        extensions: ['.js', '.ts'],
    },
    module: {
        rules: [
            { test: /\.ts?$/, loader: 'ts-loader' },
            { test: /\.css?$/, use: ['style-loader', 'css-loader'] },
            { test: /\.html?$/, loader: 'html-loader', options: { minimize: true } },
            { test: /\.js?$/, use: ['ify-loader'] },
        ],
    },
};

export default config;
