import { AxisOptions, MapOptions } from '../../src/map/options';
import { getByID } from '../../src/utils';

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

describe('MapOptions', () => {
    before(() => {
        // store karma's default HTML
        KARMA_INSERTED_HTML = document.body.innerHTML;
    });

    it('can remove itself from DOM', () => {
        const root = document.createElement('div');
        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);
        assert(root.innerHTML !== '');

        options.remove();
        assert(root.innerHTML === '');
        assert(document.body.innerHTML === KARMA_INSERTED_HTML);
    });

    it('scale label for min/max switches between linear and log', () => {
        const root = document.createElement('div');
        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);

        checkScaleLabel(options.x, 'x');
        checkScaleLabel(options.y, 'y');
        checkScaleLabel(options.z, 'z');

        function checkScaleLabel(axisOptions: AxisOptions, axisName: string): void {
            const min = getByID(`chsp-${axisName}-min-label`);
            const max = getByID(`chsp-${axisName}-max-label`);
            const selectElement = getByID<HTMLSelectElement>(`chsp-${axisName}-scale`);

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
    });
});
