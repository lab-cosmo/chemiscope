import { default as setupJSDOM } from '../jsdom';
import { assert } from 'chai';
setupJSDOM();

import { default as $ } from 'jquery';
(global as any).$ = $;

import { ViewersGrid } from '../../src/structure/index';
import { EnvironmentIndexer } from '../../src/indexer';
import { MoleculeViewer } from '../../src/structure/widget';

const DUMMY_STRUCTURES = [
    {
        size: 2,
        names: ['X', 'Y'],
        x: [0, 1],
        y: [0, 1],
        z: [0, 1],
    },
];

const DUMMY_ENVIRONMENTS = {};

describe('ViewersGrid', () => {
    before(() => {});

    it('can remove itself from DOM', () => {
        const root = document.createElement('div');
        document.body.appendChild(root);
        root.id = 'gridID';
        const indexer = new EnvironmentIndexer('structure', DUMMY_STRUCTURES);
        const grid = new ViewersGrid('gridID', indexer, DUMMY_STRUCTURES);

        assert(root.innerHTML !== '');
        assert(document.body.innerHTML !== '');

        //grid.remove();
        assert(document.body.innerHTML === '');
        document.removeChild(root);
        assert(root.innerHTML === '');
    });
});
