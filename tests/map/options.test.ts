import { MapOptions, AxisOptions } from '../../src/map/options';
import { getByID } from '../../src/utils';

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

    it('can remove itself from DOM', () => {
        const root = document.createElement('div');
        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);
        assert(root.innerHTML !== '');
        assert(document.body.innerHTML !== '');

        options.remove();
        assert(document.body.innerHTML === '');
        assert(root.innerHTML === '');
    });

    it('scale label switches between linear and log', () => {
        const root = document.createElement('div');
        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);

        assertScaleLabel(options.x, 'x');
        assertScaleLabel(options.y, 'y');
        assertScaleLabel(options.z, 'z');

        function assertScaleLabel(axisOptions: AxisOptions, axisName: string): void {
            const min = getByID(`chsp-${axisName}-min-label`);
            const max = getByID(`chsp-${axisName}-max-label`);

            // change from linear (default) to log scale
            options.setLogLabel(axisOptions, axisName);
            assert(min.innerHTML === 'min: 10^');
            assert(max.innerHTML === 'max: 10^');

            // change back from log to linear
            options.setLogLabel(axisOptions, axisName);
            assert(min.innerHTML === 'min:');
            assert(max.innerHTML === 'max:');
        }
    });
});
