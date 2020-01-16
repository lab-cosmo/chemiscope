import assert from 'assert';

import {Property, Target} from '../dataset';
import {Indexes, EnvironmentIndexer} from '../indexer';

interface TableProperty {
    values: number[] | string[];
    cell: HTMLTableDataCellElement;
}

/// A table to display the properties of the current selected
/// structure/environement
export class Table {
    private _target: Target;
    private _indexer: EnvironmentIndexer;

    private _header: HTMLTableHeaderCellElement;
    private _properties: TableProperty[];

    constructor(root: HTMLElement, target: Target, collapseID: string, properties: {[name: string]: Property}, indexer: EnvironmentIndexer) {
        const template = document.createElement('template');
        template.innerHTML = `<div class="collapse" id=${collapseID}>
        <div class="skv-properties">
            <table class="table table-striped table-sm">
                <thead><th colspan=2 style="text-align: center;"></th></thead>
                <tbody></tbody>
            </table>
        </div></div>`;
        const group = template.content.firstChild! as HTMLElement;
        root.appendChild(group);

        this._header = group.querySelector('th')!;
        this._target = target;
        this._indexer = indexer;
        this._properties = [];

        const tbody = group.querySelector('tbody')!;
        for (const name in properties) {
            const tr = document.createElement('tr');
            const td = document.createElement('td');
            td.innerText = name;
            tr.appendChild(td);
            const cell = document.createElement('td');
            tr.appendChild(cell);

            tbody.appendChild(tr);
            this._properties.push({
                values: properties[name].values,
                cell: cell,
            })
        }

        this.show({structure: 0, atom: 0});
    }

    public show(indexes: Indexes) {
        let id, index;
        if (this._target === 'structure') {
            id = indexes.structure + 1;
            index = indexes.structure;
        } else {
            assert(this._target === 'atom');
            assert(indexes.atom !== undefined);
            id = indexes.atom! + 1;
            index = this._indexer.environment(indexes);
        }

        this._header.innerText = `Properties for ${this._target} ${id}`;
        for (const s of this._properties) {
            s.cell.innerText = s.values[index].toString()
        }
    }
}
