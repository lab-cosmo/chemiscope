import js from '@eslint/js';
import prettierEslint from 'eslint-config-prettier';
import tseslint from 'typescript-eslint';

export default tseslint.config({
    extends: [js.configs.recommended, prettierEslint, ...tseslint.configs.strictTypeChecked],
    ignores: ['**/app/', '**/dist/'],
    rules: {
        'arrow-parens': 'error',
        eqeqeq: 'error',
        'grouped-accessor-pairs': 'error',
        'no-console': 'error',
        'no-extra-parens': 'off',
        'no-loss-of-precision': 'error',
        'no-sequences': 'error',
        'no-template-curly-in-string': 'error',
        'no-throw-literal': 'error',
        'no-unused-expressions': 'error',
        'no-useless-concat': 'error',
        'no-var': 'error',
        radix: 'error',
        semi: 'error',
        'sort-imports': [
            'error',
            {
                // only sort inside groups, not between separate import lines
                ignoreDeclarationSort: true,
            },
        ],
        // disabled rules
        'sort-keys': 'off',
        '@typescript-eslint/no-inferrable-types': 'off',
        '@typescript-eslint/no-empty-function': 'off',
        '@typescript-eslint/no-deprecated': 'off',
        '@typescript-eslint/restrict-template-expressions': 'off',
        '@typescript-eslint/no-confusing-void-expression': 'off',
        '@typescript-eslint/no-unnecessary-condition': 'off',
    },
    languageOptions: {
        parserOptions: {
            projectService: {
                allowDefaultProject: [
                    'eslint.config.mjs',
                    'src/map/plotly/plotly-scatter.js',
                    'src/map/plotly/markers3d.js',
                    'python/webpack.config.labextension.js',
                    'python/webpack.config.nbextension.js',
                ],
            },
            tsconfigRootDir: import.meta.dirname,
        },
    },
});
