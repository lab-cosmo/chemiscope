import assert from 'assert';
import {Property} from './dataset';
import {EnvironmentIndexer, Indexes} from './indexer';

interface TableProperty {
    values: number[] | string[];
    cell: HTMLTableDataCellElement;
}

/// A table to display the properties of the current selected structure/environement
export class PropertiesTable {
    private _root: HTMLElement;
    private _indexer: EnvironmentIndexer;

    private _atomsTable: HTMLTableElement;
    private _atomHeader: HTMLTableHeaderCellElement;
    private _atoms: TableProperty[];

    private _structureHeader: HTMLTableHeaderCellElement;
    private _structures: TableProperty[];

    constructor(id: string, properties: {[name: string]: Property}, indexer: EnvironmentIndexer) {
        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`)
        }
        this._root = root;
        this._indexer = indexer;

        this._root.innerHTML = `
        <div class="skv-properties">
            <table class="table table-striped table-sm">
                <thead><th colspan=2 style="text-align: center;"></th></thead>
                <tbody></tbody>
            </table>
            <table class="table table-striped table-sm">
                <thead><th colspan=2 style="text-align: center;"></th></thead>
                <tbody></tbody>
            </table>
        </div>`;

        this._atomsTable = this._root.querySelectorAll('table')[0];
        this._atomsTable.style.display = (this._indexer.mode === 'atom') ? 'auto' : 'none';

        const allHeaders = this._root.querySelectorAll('th');
        this._atomHeader = allHeaders[0]
        this._structureHeader = allHeaders[1];

        const allBody = this._root.querySelectorAll('tbody');
        this._atoms = []
        this._structures = []
        for (const name in properties) {
            const tr = document.createElement('tr');
            const td = document.createElement('td');
            td.innerText = name;
            tr.appendChild(td);
            const cell = document.createElement('td');
            tr.appendChild(cell);

            if (properties[name].target === 'atom') {
                allBody[0].appendChild(tr);
                this._atoms.push({
                    values: properties[name].values,
                    cell: cell,
                })
            } else {
                assert(properties[name].target === 'structure');
                allBody[1].appendChild(tr);
                this._structures.push({
                    values: properties[name].values,
                    cell: cell,
                })
            }
        }

        this.show({structure: 0, atom: 0});
    }

    public show({structure, atom}: Indexes) {
        this._structureHeader.innerText = `Properties for structure ${structure}`;
        for (const s of this._structures) {
            s.cell.innerText = s.values[structure].toString()
        }

        if (this._indexer.mode === 'atom') {
            this._atomHeader.innerText = `Properties for atom ${atom}`;
            const env = this._indexer.environment({structure, atom});
            for (const a of this._atoms) {
                a.cell.innerText = a.values[env].toString()
            }
        } else {
            assert(atom === undefined || atom === 0);
        }
    }
}
