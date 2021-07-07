import { ViewersGrid } from '../../src/structure/index';

import { default as setupJSDOM } from '../jsdom';
import { assert } from 'chai';

describe('MapOptions', () => {
    before(() => {
        setupJSDOM();
    });

    it('can remove itself from DOM', () => {
        const root = document.createElement('div');
        const grid = new ViewersGrid(
            config.structure,
            this._indexer,
            dataset.structures,
            dataset.environments
        );

        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);
        assert(root.innerHTML !== '');
        assert(document.body.innerHTML !== '');

        options.remove();
        assert(document.body.innerHTML === '');
        assert(root.innerHTML === '');
    });
});
