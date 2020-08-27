import { MapOptions } from '../../src/map/options';

import { default as setupJSDOM } from '../jsdom';
import { assert } from 'chai';

const DUMMY_PROPERTIES = {
    first: {
        values: [],
    },
    second: {
        values: [],
    },
};
const DUMMY_CALLBACK = () => {
    return { top: 0, left: 0 };
};

describe('MapOptions', () => {
    before(() => {
        setupJSDOM();
    });

    it('validate new values for properties');

    it('creates HTML elements with unique id');

    it('can save settings');

    it('can apply saved settings');

    it('can remove itself from DOM', () => {
        const root = document.createElement('div');
        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);
        assert(root.innerHTML !== '');
        assert(document.body.innerHTML !== '');

        options.remove();
        assert(document.body.innerHTML === '');
        assert(root.innerHTML === '');
    });
});
