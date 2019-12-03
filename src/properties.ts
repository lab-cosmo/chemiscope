import {Property} from './dataset'

/// A table to display the properties of the current selected structure/environement
export class PropertiesTable {
    private _root: HTMLElement;
    private _header: HTMLTableHeaderCellElement;
    private _cells: {
        [name: string]: HTMLTableDataCellElement
    }

    private _structureProperties: {
        [name: string]: number[] | string[]
    }
    private _atomProperties: {
        [name: string]: number[] | string[]
    }

    constructor(id: string, properties: {[name: string]: Property}) {
        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`)
        }
        this._root = root;

        this._structureProperties = {}
        this._atomProperties = {}
        for (const name in properties) {
            if (properties[name].target === 'structure') {
                this._structureProperties[name] = properties[name].values
            } else if (properties[name].target === 'atom') {
                this._atomProperties[name] = properties[name].values
            }
        }

        this._root.innerHTML = `
        <div class="skv-properties">
            <table class="table table-striped table-sm">
                <thead><th colspan=2 style="text-align: center;"></th></thead>
                <tbody></tbody>
            </table>
        </div>`;
        this._header = this._root.querySelector('th')!;

        this._cells = {};
        const tbody = this._root.querySelector('tbody')!;
        for (const name in this._structureProperties) {
            const tr = document.createElement('tr');

            const tc = document.createElement('td');
            tc.innerText = name;
            tr.appendChild(tc);

            this._cells[name] = document.createElement('td');
            tr.appendChild(this._cells[name]);

            tbody.appendChild(tr);
        }
        // TODO: handle atom properties

        this.display(0);
    }

    public display(index: number) {
        const centerType = 'structure';
        this._header.innerText = `Properties for ${centerType} ${index}`;

        for (const name in this._structureProperties) {
            this._cells[name].innerText = this._structureProperties[name][index].toString()
        }
    }
}
