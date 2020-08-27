import { enumerate, getFirstKey, getNextColor } from '../src/utils';

import { assert } from 'chai';

describe('enumerate', () => {
    it('enumerates values in an iterator', () => {
        const array = ['foo', 'bar', 'baz', 5673, 0];
        let i = 0;
        for (const [index, value] of enumerate(array)) {
            assert(i === index);
            assert(array[i] === value);
            i++;
        }
    });
});

describe('getNextColor', () => {
    const MAX_COLORS = 9;

    it('gives different colors each time', () => {
        const colors: string[] = [];
        for (let i = 0; i < MAX_COLORS; i++) {
            const newColor = getNextColor(colors);
            assert(colors.indexOf(newColor) === -1);
            colors.push(newColor);
        }
    });

    it('can only give up to MAX_COLORS values', () => {
        const colors: string[] = [];
        for (let i = 0; i < MAX_COLORS; i++) {
            colors.push(getNextColor(colors));
        }
        assert.throws(() => getNextColor(colors), 'required more colors than available');
    });
});

describe('getFirstKey', () => {
    it('gives the first key of a Map', () => {
        const map = new Map<string, number>();
        map.set('first', 33);
        map.set('second', 44);
        map.set('third', 55);

        assert(getFirstKey(map) === 'first');

        map.delete('first');
        assert(getFirstKey(map) === 'second');
    });

    it('can exclude some values', () => {
        const map = new Map<string, number>();
        map.set('first', 33);
        map.set('second', 44);
        map.set('third', 55);

        assert(getFirstKey(map, 'first') === 'second');
    });
});
