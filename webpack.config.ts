import path from "path";
import webpack from "webpack";

const config: webpack.Configuration = {
  mode: "development",
  entry: "./src/index.ts",
  output: {
      path: path.resolve(__dirname, "dist"),
      filename: "sketchviz.min.js",
      libraryTarget: "var",
      library: "Sketchviz",
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
        ],
    },
};

export default config;
