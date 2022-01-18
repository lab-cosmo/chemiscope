/* eslint-disable @typescript-eslint/no-non-null-assertion */

import { EnvironmentIndexer } from '../src/indexer';

import { assert } from 'chai';

const DUMMY_STRUCTURES = [
    {
        size: 2,
        data: '',
    },
    {
        size: 3,
        data: '',
    },
];

describe('EnvironmentIndexer', () => {
    it('can give data about the structure sizes', () => {
        const indexer = new EnvironmentIndexer('structure', DUMMY_STRUCTURES);
        assert(indexer.structuresCount() === 2);
        assert(indexer.atomsCount(0) === 2);
        assert(indexer.atomsCount(1) === 3);

        assert.deepEqual(indexer.activeStructures(), [0, 1]);
        assert.deepEqual(indexer.activeAtoms(0), [0, 1]);
        assert.deepEqual(indexer.activeAtoms(1), [0, 1, 2]);
    });

    it('can work in structure mode', () => {
        const indexer = new EnvironmentIndexer('structure', DUMMY_STRUCTURES);

        let indexes = indexer.from_structure_atom(0, undefined);
        assert(indexes !== undefined);
        assert(indexes!.environment === 0);
        assert(indexes!.structure === 0);
        assert(indexes!.atom === undefined);

        indexes = indexer.from_structure_atom(1, undefined);
        assert(indexes !== undefined);
        assert(indexes!.environment === 1);
        assert(indexes!.structure === 1);
        assert(indexes!.atom === undefined);

        indexes = indexer.from_environment(0);
        assert(indexes.environment === 0);
        assert(indexes.structure === 0);
        assert(indexes.atom === undefined);

        indexes = indexer.from_environment(1);
        assert(indexes.environment === 1);
        assert(indexes.structure === 1);
        assert(indexes.atom === undefined);

        assert.throws(() => {
            // unable to get the environment for an atom in 'structure' mode
            indexer.from_structure_atom(0, 1);
        });
    });

    it('can work in atom mode', () => {
        const environments = [
            { structure: 0, center: 0, cutoff: 3.5 },
            { structure: 0, center: 1, cutoff: 3.5 },
            { structure: 1, center: 0, cutoff: 3.5 },
            { structure: 1, center: 1, cutoff: 3.5 },
            { structure: 1, center: 2, cutoff: 3.5 },
        ];
        const indexer = new EnvironmentIndexer('atom', DUMMY_STRUCTURES, environments);

        let indexes = indexer.from_structure_atom(0, 1);
        assert(indexes !== undefined);
        assert(indexes!.environment === 1);
        assert(indexes!.structure === 0);
        assert(indexes!.atom === 1);

        indexes = indexer.from_structure_atom(1, 2);
        assert(indexes !== undefined);
        assert(indexes!.environment === 4);
        assert(indexes!.structure === 1);
        assert(indexes!.atom === 2);

        indexes = indexer.from_environment(2);
        assert(indexes.environment === 2);
        assert(indexes.structure === 1);
        assert(indexes.atom === 0);

        indexes = indexer.from_environment(3);
        assert(indexes.environment === 3);
        assert(indexes.structure === 1);
        assert(indexes.atom === 1);
    });

    it('can work with partial list of environments', () => {
        const environments = [
            { structure: 1, center: 0, cutoff: 3.5 },
            { structure: 1, center: 2, cutoff: 3.5 },
        ];
        const indexer = new EnvironmentIndexer('atom', DUMMY_STRUCTURES, environments);

        assert(indexer.structuresCount() === 2);
        assert(indexer.atomsCount(0) === 2);
        assert(indexer.atomsCount(1) === 3);

        assert.deepEqual(indexer.activeStructures(), [1]);
        assert.deepEqual(indexer.activeAtoms(0), []);
        assert.deepEqual(indexer.activeAtoms(1), [0, 2]);

        let indexes = indexer.from_structure_atom(0, 1);
        assert(indexes === undefined);

        indexes = indexer.from_structure_atom(1, 2);
        assert(indexes !== undefined);
        assert(indexes!.environment === 1);
        assert(indexes!.structure === 1);
        assert(indexes!.atom === 2);

        indexes = indexer.from_environment(0);
        assert(indexes.environment === 0);
        assert(indexes.structure === 1);
        assert(indexes.atom === 0);
    });
});
