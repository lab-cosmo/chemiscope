import { MapData } from '../../src/map/data';

import { assert } from 'chai';

describe('MapData', () => {
    it('successfully convert from string to numeric properties', () => {
        const data = new MapData({
            first: { values: ['1', '5', '1', '3', '5'], target: 'atom' },
        });

        assert.deepEqual(data.atom.first.values, [0, 1, 0, 2, 1]);
        assert.deepEqual(data.atom.first.string?.string(0), '1');
        assert.deepEqual(data.atom.first.string?.string(1), '5');
        assert.deepEqual(data.atom.first.string?.string(2), '3');
    });

    it('separates atom and structure properties', () => {
        const data = new MapData({
            first: { values: [1], target: 'atom' },
            second: { values: [2], target: 'atom' },
            third: { values: [1], target: 'structure' },
            fourth: { values: [3], target: 'atom' },
            fifth: { values: [5], target: 'structure' },
        });

        assert(Object.keys(data.atom).length === 3);
        assert(Object.keys(data.structure).length === 2);

        for (const name of ['first', 'second', 'fourth']) {
            assert(name in data.atom);
        }

        for (const name of ['third', 'fifth']) {
            assert(name in data.structure);
        }
    });

    it('gives the maximal number of symbols for string properties', () => {
        const data = new MapData({
            first: { values: ['1', '2', '3', '4', '5'], target: 'atom' },
            second: { values: [1, 2, 3, 4, 5], target: 'atom' },
            third: { values: ['1', '2', '3'], target: 'structure' },
            fourth: { values: [1, 2, 3, 4, 5], target: 'atom' },
            fifth: { values: [1, 2, 3], target: 'structure' },
        });

        assert(data.maxSymbols === 5);
    });

    it('fails if some properties have different lengths', () => {
        assert.throws(() => {
            new MapData({
                first: { values: [1, 2, 3], target: 'atom' },
                second: { values: [1, 2, 3, 4], target: 'atom' },
            });
        }, "atom property 'second' do not have the same size as the first property 'first': expected 3, got 4");

        assert.throws(() => {
            new MapData({
                foo: { values: [1, 2, 3], target: 'structure' },
                bar: { values: [1, 2], target: 'structure' },
            });
        }, "structure property 'bar' do not have the same size as the first property 'foo': expected 3, got 2");
    });
});
