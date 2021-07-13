import { assert } from 'chai';

import { ViewersGrid } from '../../src/structure/index';
import { EnvironmentIndexer } from '../../src/indexer';

let KARMA_INSERTED_HTML: string;

const DUMMY_STRUCTURES = [
    {
        size: 2,
        names: ['X', 'Y'],
        x: [0, 1],
        y: [0, 1],
        z: [0, 1],
    },
];

describe('ViewersGrid', () => {
    before(() => {
        // store karma's default HTML
        KARMA_INSERTED_HTML = document.body.innerHTML;
    });

    it('can remove itself from DOM', () => {
        const root = document.createElement('div');
        document.body.appendChild(root);
        root.id = 'gridID';
        const indexer = new EnvironmentIndexer('structure', DUMMY_STRUCTURES);
        const grid = new ViewersGrid('gridID', indexer, DUMMY_STRUCTURES);

        assert(root.innerHTML !== '');
        assert(document.body.innerHTML !== '');

        grid.remove();
        assert(root.innerHTML === '');
        root.remove();
        assert(document.body.innerHTML === KARMA_INSERTED_HTML);
    });
});
