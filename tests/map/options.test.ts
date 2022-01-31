import { AxisOptions, MapOptions } from '../../src/map/options';

import { assert } from 'chai';

let KARMA_INSERTED_HTML: string;

const DUMMY_PROPERTIES = {
    first: {
        values: [],
    },
    second: {
        values: [],
    },
    third: {
        values: [],
    },
};

const DUMMY_CALLBACK = () => {
    return { top: 0, left: 0 };
};

function traverseDOM(element: Element | ShadowRoot, callback: (element: Element) => void) {
    if (element.children) {
        const elements = element.children;
        for (let i = 0; i < elements.length; i++) {
            traverseDOM(elements[i], callback);
        }
    }

    if (element instanceof Element) {
        callback(element);
    }
}

function createShadow(): ShadowRoot {
    const shadowElement = document.createElement('div');
    const shadow = shadowElement.attachShadow({ mode: 'open' });

    return shadow;
}


describe('MapOptions', () => {
    before(() => {
        // store karma's default HTML
        KARMA_INSERTED_HTML = document.body.innerHTML;
    });

    it('can remove itself from DOM', () => {
        const root = createShadow();
        const options = new MapOptions(
            root,
            DUMMY_PROPERTIES,
            DUMMY_CALLBACK
        );
        assert(root.innerHTML !== '');

        options.remove();
        assert(root.innerHTML === '');
        assert(document.body.innerHTML === KARMA_INSERTED_HTML);
    });

    it('scale label for min/max switches between linear and log', () => {
        const root = createShadow();

        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);

        function checkScaleLabel(axisOptions: AxisOptions, axisName: string): void {
            const min = options.getModalElement(`map-${axisName}-min-label`);
            const max = options.getModalElement(`map-${axisName}-max-label`);
            const selectElement = options.getModalElement<HTMLSelectElement>(`map-${axisName}-scale`);

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

    it('axis property changes when modified by the user', () => {
        const root = createShadow();
        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);

        function checkPropertySelect(axisOptions: AxisOptions, axisName: string) {
            const selectElement = options.getModalElement<HTMLSelectElement>(`map-${axisName}-property`);
            selectElement.value = 'first';
            selectElement.dispatchEvent(new window.Event('change'));
            assert(selectElement.value !== 'second');
            assert(axisOptions.property.value !== 'second');
            selectElement.value = 'second';
            selectElement.dispatchEvent(new window.Event('change'));
            assert(axisOptions.property.value === 'second');
        }

        checkPropertySelect(options.x, 'x');
        checkPropertySelect(options.y, 'y');
        checkPropertySelect(options.z, 'z');
        checkPropertySelect(options.color, 'color');

        options.remove();
    });

    it('axis range changes when min/max are modified', () => {
        const root = createShadow();
        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);

        function checkRangeSelect(axisOptions: AxisOptions, axisName: string) {
            const minSelectElement = options.getModalElement<HTMLSelectElement>(`map-${axisName}-min`);
            const maxSelectElement = options.getModalElement<HTMLSelectElement>(`map-${axisName}-max`);

            assert(minSelectElement.value !== '0.01');
            assert(axisOptions.min.value !== 0.01);
            minSelectElement.value = '0.01';
            minSelectElement.dispatchEvent(new window.Event('change'));
            assert(axisOptions.min.value === 0.01);

            assert(maxSelectElement.value !== '1.01');
            assert(axisOptions.max.value !== 1.01);
            maxSelectElement.value = '1.01';
            maxSelectElement.dispatchEvent(new window.Event('change'));
            assert(axisOptions.max.value === 1.01);
        }

        checkRangeSelect(options.x, 'x');
        checkRangeSelect(options.y, 'y');
        checkRangeSelect(options.z, 'z');
        checkRangeSelect(options.color, 'color');

        options.remove();
    });

    it('property used for size changes when modified by the user', () => {
        const root = createShadow();

        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);
        const selectElement = options.getModalElement<HTMLSelectElement>('map-size-property');

        assert(selectElement.value !== 'second');
        assert(options.size.property.value !== 'second');
        selectElement.value = 'second';
        selectElement.dispatchEvent(new window.Event('change'));
        assert(options.size.property.value === 'second');

        options.remove();
    });

    it('size factor changes when the slider is modified by the user', () => {
        const root = createShadow();
        const options = new MapOptions(root, DUMMY_PROPERTIES, DUMMY_CALLBACK);
        const sliderElement = options.getModalElement<HTMLInputElement>('map-size-factor');

        assert(sliderElement.value !== '100');
        assert(options.size.factor.value !== 100);
        sliderElement.value = '100';
        sliderElement.dispatchEvent(new window.Event('change'));
        assert(options.size.factor.value === 100);

        options.remove();
    });
});
