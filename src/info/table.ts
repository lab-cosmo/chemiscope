/**
 * @packageDocumentation
 * @module info
 */

import assert from 'assert';

import { Property, Target } from '../dataset';
import { Indexes } from '../indexer';

interface TableProperty {
    values: number[] | string[];
    cell: HTMLTableDataCellElement;
}

/** @hidden
 * A table to display the properties of the current selected structure/environement
 */
export class Table {
    private _target: Target;

    private _header: HTMLTableHeaderCellElement;
    private _properties: TableProperty[];

    /**
     * Create and append a new table inside the given `HTMLElement`.
     *
     * @param root       where to append the new table
     * @param target     is this table related to atom or structure
     * @param collapseID HTML id to use for the root div with class=collapse
     * @param properties properties to display in this table.
     */
    constructor(
        root: HTMLElement,
        target: Target,
        collapseID: string,
        properties: { [name: string]: Property }
    ) {
        const template = document.createElement('template');
        template.innerHTML = `<div class="collapse chsp-info-table" id=${collapseID} data-parent='#info-tables'>
        <div class="chsp-properties-table">
            <table class="table table-striped table-sm">
                <thead><th colspan=2 style="text-align: center;"></th></thead>
                <tbody></tbody>
            </table>
        </div></div>`;
        const group = template.content.firstChild as HTMLElement;
        root.appendChild(group);

        this._header = group.querySelector('th') as HTMLTableHeaderCellElement;
        this._target = target;
        this._properties = [];

        const tbody = group.querySelector('tbody') as HTMLTableSectionElement;
        for (const name in properties) {
            const tr = document.createElement('tr');
            const td = document.createElement('td');

            //  add the units to the property if it exists, this is identical to the _title in ../map/map.ts
            const units = properties[name].units;
            let title = name;
            if (units !== undefined) {
                title += ` [${units}]`;
            }
            td.innerText = title;

            // add a tooltip containing the description of the property and underline if it exists
            const description = properties[name].description;
            if (description !== undefined) {
                td.style.borderBottom = '1px dotted #00f';
                td.style.display = 'inline';
                td.title = description;
            }

            tr.appendChild(td);
            const cell = document.createElement('td');
            tr.appendChild(cell);

            tbody.appendChild(tr);
            this._properties.push({
                cell: cell,
                values: properties[name].values,
            });
        }

        this.show({ environment: 0, structure: 0, atom: 0 });
    }

    /**
     * Show the properties for the given `environment`, corresponding to the
     * given structure/atom `indexes`
     */
    public show(indexes: Indexes): void {
        let displayId;
        let index;
        if (this._target === 'structure') {
            displayId = indexes.structure + 1;
            index = indexes.structure;
        } else {
            assert(this._target === 'atom');
            assert(indexes.atom !== undefined);
            displayId = indexes.atom + 1;
            index = indexes.environment;
        }

        this._header.innerText = `Properties for ${this._target} ${displayId}`;
        for (const s of this._properties) {
            s.cell.innerText = s.values[index].toString();
        }
    }
}
