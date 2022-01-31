import { PropertiesMap } from '../../src/map';
import { MarkerData } from '../../src/map/marker';
import { EnvironmentIndexer } from '../../src/indexer';
import { Property } from '../../src/dataset';
import { GUID } from '../../src/utils';

import { assert } from 'chai';

let KARMA_INSERTED_HTML: string;

const DUMMY_PROPERTIES = {
    first: {
        target: 'structure',
        values: [1.1, 1.2, 0.6],
    } as Property,
    second: {
        target: 'structure',
        values: [2.1, 2.2, 1.4],
    } as Property,
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
        const map = new PropertiesMap({ element: root, settings: {} }, indexer, DUMMY_PROPERTIES);

        assert(root.innerHTML !== '');
        assert(document.body.innerHTML !== KARMA_INSERTED_HTML);

        map.remove();
        assert(root.innerHTML === '');

        assert(document.body.innerHTML === KARMA_INSERTED_HTML);
    });

    it('the color range resets when clicked', () => {
        const root = document.createElement('div');
        const indexer = new EnvironmentIndexer('structure', DUMMY_STRUCTURES);
        const map = new PropertiesMap({ element: root, settings: {} }, indexer, DUMMY_PROPERTIES);
        const options = map['_options'];

        // make sure we have colors in the dataset (with 2 properties the
        // default is to have a constant color)
        options.color.property.value = 'first';

        const minSelectElement = options.getModalElement<HTMLSelectElement>(`map-color-min`);
        const maxSelectElement = options.getModalElement<HTMLSelectElement>(`map-color-max`);
        const originalMin = options.color.min.value;
        const originalMax = options.color.max.value;

        maxSelectElement.value = '3.0';
        maxSelectElement.dispatchEvent(new window.Event('change'));
        assert(options.color.max.value === 3.0);

        minSelectElement.value = '1.0';
        minSelectElement.dispatchEvent(new window.Event('change'));
        assert(options.color.min.value === 1.0);

        const resetButton = map['_colorReset'];
        resetButton.dispatchEvent(new window.Event('click'));

        assert(options.color.min.value === originalMin);
        assert(options.color.max.value === originalMax);

        map.remove();
    });
});

describe('map markers', () => {
    let MAP: PropertiesMap;
    const firstGUID = 'foo' as GUID;
    const secondGUID = 'bar' as GUID;
    const getMarker = (id: GUID): MarkerData => {
        const marker = MAP['_selected'].get(id);
        assert(marker !== undefined);
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
        return marker!;
    };

    beforeEach(() => {
        const root = document.createElement('div');
        root.style.height = '200px';
        root.style.width = '200px';
        const indexer = new EnvironmentIndexer('structure', DUMMY_STRUCTURES);
        MAP = new PropertiesMap({ element: root, settings: {} }, indexer, DUMMY_PROPERTIES);

        document.body.append(root);
    });

    afterEach(() => {
        MAP.remove();
    });

    it('can add and remove markers', () => {
        assert(MAP['_selected'].size === 0);
        assert(MAP['_active'] === undefined);

        MAP.addMarker(firstGUID, 'red', { structure: 0, environment: 0 });
        assert(MAP['_selected'].size === 1);
        assert(MAP['_active'] === firstGUID);
        assert(getMarker(firstGUID).color === 'red');

        MAP.addMarker(secondGUID, 'blue', { structure: 1, environment: 1 });
        assert(MAP['_selected'].size === 2);
        assert(MAP['_active'] === secondGUID);
        assert(getMarker(secondGUID).color === 'blue');

        MAP.setActive(firstGUID);
        assert(MAP['_selected'].size === 2);
        assert(MAP['_active'] === firstGUID);

        MAP.removeMarker(firstGUID);
        assert(MAP['_selected'].size === 1);
        // the remaining marker is made active
        assert(MAP['_active'] === secondGUID);

        MAP.removeMarker(secondGUID);
        assert(MAP['_selected'].size === 0);
        assert(MAP['_active'] === undefined);
    });
});
