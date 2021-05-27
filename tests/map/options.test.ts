import { AxisOptions, MapOptions } from '../../src/map/options';
import { GUID, getByID } from '../../src/utils';

import { assert } from 'chai';

let KARMA_INSERTED_HTML: string;

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

function traverseDOM(element: Element, callback: (element: Element) => void) {
    if (element.children) {
        const elements = element.children;
        for (let i = 0; i < elements.length; i++) {
            traverseDOM(elements[i], callback);
        }
    }
    callback(element);
}

describe('MapOptions', () => {
    before(() => {
        // store karma's default HTML
        KARMA_INSERTED_HTML = document.body.innerHTML;
    });

    it('can remove itself from DOM', () => {
        const root = document.createElement('div');
        const options = new MapOptions(
            root,
            'this-is-my-id' as GUID,
            DUMMY_PROPERTIES,
            DUMMY_CALLBACK
        );
        assert(root.innerHTML !== '');

        options.remove();
        assert(root.innerHTML === '');
        assert(document.body.innerHTML === KARMA_INSERTED_HTML);
    });

    it('scale label for min/max switches between linear and log', () => {
        const root = document.createElement('div');
        const guid = 'guid' as GUID;

        const options = new MapOptions(root, guid, DUMMY_PROPERTIES, DUMMY_CALLBACK);

        function checkScaleLabel(axisOptions: AxisOptions, axisName: string): void {
            const min = getByID(`guid-map-${axisName}-min-label`);
            const max = getByID(`guid-map-${axisName}-max-label`);
            const selectElement = getByID<HTMLSelectElement>(`guid-map-${axisName}-scale`);

            // change from linear (default) to log scale
            selectElement.value = 'log';
            selectElement.dispatchEvent(new window.Event('change'));
            options.setLogLabel(axisOptions, axisName);
            assert(min.innerHTML === 'min: 10^');
            assert(max.innerHTML === 'max: 10^');

            // change back from log to linear
            selectElement.value = 'linear';
            selectElement.dispatchEvent(new window.Event('change'));
            options.setLogLabel(axisOptions, axisName);
            assert(min.innerHTML === 'min:');
            assert(max.innerHTML === 'max:');
        }

        checkScaleLabel(options.x, 'x');
        checkScaleLabel(options.y, 'y');
        checkScaleLabel(options.z, 'z');

        options.remove();
    });

    it('has a unique id in the page', () => {
        const root = document.createElement('div');

        const guid = 'this-is-my-id' as GUID;
        const options = new MapOptions(root, guid, DUMMY_PROPERTIES, DUMMY_CALLBACK);
        traverseDOM(document.body, (element) => {
            if (element.id) {
                assert(element.id.includes(guid));
            }
        });

        options.remove();
    });
});
