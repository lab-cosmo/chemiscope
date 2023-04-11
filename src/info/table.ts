/**
 * @packageDocumentation
 * @module info
 */

import assert from 'assert';

import { Parameter, Property, Target } from '../dataset';
import { Indexes } from '../indexer';
import { plotMultiDimensionalProperties } from './plotting';
import { getByID } from '../utils';

/**
 * TableProperty holds the objects to show the properties in the info bar
 * values and cell are common between all types of properties
 * parameter: the parameter associated to the multidimensional property,
 *            acts as the xaxis in the plot
 * xlabel: xlabel of the plot
 * ylabel: ylabel of the plot
 */
interface TableProperty {
    values: number[] | string[] | number[][];
    cell: HTMLTableCellElement;
    parameter?: number[]; // for multidimensional properties
    xlabel?: string; // for multidimensional properties
    ylabel?: string; // for multidimensional properties
}

/** @hidden
 * A table to display the properties of the current selected structure/environement
 */
export class Table {
    private _target: Target;

    private _root: HTMLElement;
    private _header: HTMLTableCellElement;
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
        properties: { [name: string]: Property },
        parameters?: { [name: string]: Parameter }
    ) {
        const template = document.createElement('template');
        template.innerHTML = `<div class="collapse chsp-info-table" id=${collapseID} data-bs-parent='#info-tables'>
        <div class="chsp-properties-table">
            <table class="table table-striped table-sm">
                <thead><th colspan=2 style="text-align: center;"></th></thead>
                <tbody></tbody>
            </table>
        </div></div>`;
        const group = template.content.firstChild as HTMLElement;
        root.appendChild(group);
        this._root = root;

        this._header = group.querySelector('th') as HTMLTableCellElement;
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
                title += `/${units}`;
            }

            // add a tooltip containing the description of the property and underline if it exists
            const description = properties[name].description;
            if (description !== undefined) {
                td.innerHTML = `<span style="border-bottom: 1px dotted #00f; cursor: help"></span>`;

                const span = td.firstChild as HTMLSpanElement;
                span.innerText = title;
                span.setAttribute('title', description);
            } else {
                td.innerText = title;
            }

            tr.appendChild(td);
            const cell = document.createElement('td');
            tr.appendChild(cell);
            tbody.appendChild(tr);

            const propertyParameter = properties[name].parameters;
            if (typeof propertyParameter === 'undefined') {
                this._properties.push({
                    cell: cell,
                    values: properties[name].values,
                });
            } else if (parameters && typeof propertyParameter[0] === 'string') {
                let xlabel = propertyParameter[0];
                const parameterUnits = parameters[propertyParameter[0]].units as string;
                if (parameterUnits !== undefined) {
                    xlabel += `/${parameterUnits}`;
                }

                const plotHolder = document.createElement('div');
                plotHolder.style.display = 'block';
                plotHolder.style.width = '100%';

                cell.appendChild(plotHolder);

                this._properties.push({
                    cell: cell,
                    values: properties[name].values,
                    parameter: parameters[propertyParameter[0]].values,
                    xlabel: xlabel,
                    ylabel: title,
                });

                // add show/hide button to td
                const button = document.createElement('button');
                button.classList.add('btn', 'btn-secondary', 'btn-sm', 'chsp-toggle-plot-btn');
                button.textContent = 'Show/Hide';
                button.onclick = () => {
                    if (plotHolder.style.display === 'block') {
                        plotHolder.style.display = 'none';
                    } else {
                        plotHolder.style.display = 'block';
                    }
                };
                td.appendChild(button);
            }
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
            if (!Array.isArray(s.values[index])) {
                s.cell.innerText = s.values[index].toString();
            } else {
                // now we plot!!
                const widthPlotCell = this._root.offsetWidth / 1.5;

                assert(s.cell.firstElementChild !== null);

                plotMultiDimensionalProperties(
                    s.parameter as number[],
                    s.values[index] as number[],
                    s.cell.firstElementChild as HTMLElement,
                    widthPlotCell,
                    s.xlabel,
                    s.ylabel
                );
            }
        }
    }
}
