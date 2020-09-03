/**
 * @packageDocumentation
 * @module structure
 */

import * as linalg from './linalg';
import { default as $3Dmol } from './3dmol';
import { assignBonds } from './3Dmol/assignBonds';

import { Structure } from '../dataset';

/** @hidden
 * A UnitCell, usable to convert between cartesian and fractional coordinates
 */
class UnitCell {
    private matrix: linalg.Matrix;
    private inverse: linalg.Matrix;

    /**
     * Create a new unit cell from an array containing `[ax ay az bx by bz cx
     * cy cz]`, where **a**, **b**, and **c** are the unit cell vectors. All
     * values should be expressed in Angströms.
     */
    constructor(data: number[]) {
        if (data.length !== 9) {
            throw Error(`invalid length for cell: expected 9, got ${data.length}`);
        }

        const vx = data.slice(0, 3) as linalg.Vector3D;
        const vy = data.slice(3, 6) as linalg.Vector3D;
        const vz = data.slice(6, 9) as linalg.Vector3D;
        this.matrix = [vx, vy, vz];
        if (linalg.determinant(this.matrix) < 1e-9) {
            throw Error('invalid unit cell');
        }
        this.inverse = linalg.invert(this.matrix);
    }

    /** convert from cartesian to fractional coordinates */
    public fractional(cartesian: linalg.Vector3D): linalg.Vector3D {
        return linalg.dot(this.inverse, cartesian);
    }

    /** convert from fractional to cartesian coordinates */
    public cartesian(fractional: linalg.Vector3D): linalg.Vector3D {
        return linalg.dot(this.matrix, fractional);
    }

    /** Get the lenghts of the unitcell vectors */
    public lengths(): linalg.Vector3D {
        return [
            linalg.norm(this.matrix[0]),
            linalg.norm(this.matrix[1]),
            linalg.norm(this.matrix[2]),
        ];
    }

    /** Get the angles between the unitcell vectors, in degrees */
    public angles(): linalg.Vector3D {
        return [
            linalg.angle(this.matrix[1], this.matrix[2]),
            linalg.angle(this.matrix[0], this.matrix[2]),
            linalg.angle(this.matrix[0], this.matrix[1]),
        ];
    }
}

/** @hidden
 * Add data from the `structure` to the `model`
 *
 * @param model 3Dmol GLModel that will contain structure data
 * @param structure the structure to convert
 */
export function $3DmolStructure(model: $3Dmol.GLModel, structure: Structure): void {
    const rotation: linalg.Matrix = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ];
    if (structure.cell !== undefined) {
        const cell = new UnitCell(structure.cell);
        const [a, b, c] = cell.lengths();
        const [alpha, beta, gamma] = cell.angles();
        model.setCrystData(a, b, c, alpha, beta, gamma);
        // TODO: deal with the case where the cell matrix is not upper
        // triangular by computing the right rotation to align along x/y/z
    }

    const atoms = [];
    for (let i = 0; i < structure.size; i++) {
        const [x, y, z] = linalg.dot(rotation, [structure.x[i], structure.y[i], structure.z[i]]);
        atoms.push({
            serial: i,
            elem: structure.names[i],
            x: x,
            y: y,
            z: z,
        });
    }
    assignBonds(atoms as $3Dmol.AtomSpec[]);

    model.addAtoms(atoms);
}
