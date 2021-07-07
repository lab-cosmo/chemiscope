import { PropertiesMap } from '../../src/map';
import { EnvironmentIndexer } from '../../src/indexer';
import { Target } from '../../src/dataset';
import { getByID } from '../../src/utils';

import { assert } from 'chai';

let KARMA_INSERTED_HTML: string;

const DUMMY_PROPERTIES = {
    first: {
        target: 'structure' as Target,
        values: [1.1, 1.2],
    },
    second: {
        target: 'structure' as Target,
        values: [2.1, 2.2],
    },
    third: {
        target: 'structure' as Target,
        values: [3.1, 3.2],
    },
};

const DUMMY_MAP_SETTINGS = {
    x: {
        max: 0,
        min: 10,
        property: 'first',
        scale: 'linear',
    },
    y: {
        max: 0,
        min: 10,
        property: 'first',
        scale: 'linear',
    },
    z: {
        max: 0,
        min: 10,
        property: 'first',
        scale: 'linear',
    },
    color: {
        max: 0,
        min: 10,
        property: 'second',
        scale: 'linear',
    },
    symbol: '',
    palette: 'hsv (periodic)',
    size: {
        factor: 50,
        property: '',
    },
};

const DUMMY_STRUCTURES = [
    {
        size: 2,
        names: ['X', 'Y'],
        x: [0, 1],
        y: [0, 1],
        z: [0, 1],
    },
];

describe('Map', () => {
    before(() => {
        // store karma's default HTML
        KARMA_INSERTED_HTML = document.body.innerHTML;
    });

    it('can remove itself from DOM', () => {
        const root = document.createElement('div');
        const indexer = new EnvironmentIndexer('structure', DUMMY_STRUCTURES);
        const map = new PropertiesMap(
            { element: root, settings: DUMMY_MAP_SETTINGS },
            indexer,
            DUMMY_PROPERTIES
        );

        assert(root.innerHTML !== '');
        assert(document.body.innerHTML !== KARMA_INSERTED_HTML);

        map.remove();
        assert(root.innerHTML === '');

        assert(document.body.innerHTML === KARMA_INSERTED_HTML);
    });

    it('color range resets when clicked', () => {
        const root = document.createElement('div');
        const indexer = new EnvironmentIndexer('structure', DUMMY_STRUCTURES);
        const map = new PropertiesMap(
            { element: root, settings: DUMMY_MAP_SETTINGS },
            indexer,
            DUMMY_PROPERTIES
        );
        const options = map['_options'];

        const minSelectElement = getByID<HTMLSelectElement>(options.getId(`map-color-min`));
        const maxSelectElement = getByID<HTMLSelectElement>(options.getId(`map-color-max`));
        const originalMin = options.color.min.value;
        const originalMax = options.color.max.value;

        minSelectElement.value = '1.0';
        minSelectElement.dispatchEvent(new window.Event('change'));
        assert(options.color.min.value === 1.0);

        maxSelectElement.value = '3.0';
        maxSelectElement.dispatchEvent(new window.Event('change'));
        assert(options.color.max.value === 3.0);

        const resetButton = map['_colorReset'];
        resetButton.click();
        resetButton.dispatchEvent(new window.Event('click'));

        assert(options.color.min.value === originalMin);
        assert(options.color.max.value === originalMax);

        map.remove();
    });
});
