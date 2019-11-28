import {MapInput} from './map_data'

/// A table to display the properties of the current selected structure/environement
export class PropertiesTable {
    private _root: HTMLElement;
    private _data: MapInput;
    private _header: HTMLTableHeaderCellElement;
    private _cells: {
        [name: string]: HTMLTableDataCellElement
    }

    constructor(id: string, data: MapInput) {
        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`)
        }
        this._root = root;
        this._data = data;

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
        for (const name in this._data) {
            const tr = document.createElement('tr');

            const tc = document.createElement('td');
            tc.innerText = name;
            tr.appendChild(tc);

            this._cells[name] = document.createElement('td');
            tr.appendChild(this._cells[name]);

            tbody.appendChild(tr);
        }

        this.display(0);
    }

    public display(index: number) {
        const centerType = 'structure';
        this._header.innerText = `Properties for ${centerType} ${index}`;

        for (const name in this._data) {
            this._cells[name].innerText = this._data[name][index].toString()
        }
    }
}
