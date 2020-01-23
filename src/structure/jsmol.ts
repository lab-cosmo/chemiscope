/**
 * @packageDocumentation
 * @module structure
 */

import * as linalg from './linalg'

import {Structure} from '../dataset';

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
            throw Error('invalid length for cell: expected 9, got ' + data.length);
        }

        const vx = data.slice(0, 3) as linalg.Vector3D;
        const vy = data.slice(3, 6) as linalg.Vector3D;
        const vz = data.slice(6, 9) as linalg.Vector3D;
        this.matrix = [vx, vy, vz];
        if (linalg.determinant(this.matrix) < 1e-9) {
            throw Error("invalid unit cell");
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
 * Convert a `structure` to a format that JSmol can read: XYZ if there is no cell
 * data, and BCS if there is cell data
 *
 * @param structure the structure to convert
 * @return a string representing the structure that JSmol is able to read
 */
export function structure2JSmol(structure: Structure): string {
    let buffer = new Array();
    if (structure.cell === undefined) {
        // use XYZ format
        const natoms = structure.names.length;
        buffer.push(`${natoms}\n\n`);
        for (let i=0; i<natoms; i++) {
            buffer.push(`${structure.names[i]} ${structure.x[i]} ${structure.y[i]} ${structure.z[i]}\n`);
        }
    } else {
        // use BCS format
        const cell = new UnitCell(structure.cell);
        buffer.push(`1\n`);
        const [a, b, c] = cell.lengths();
        const [alpha, beta, gamma] = cell.angles();
        buffer.push(`${a} ${b} ${c} ${alpha} ${beta} ${gamma}\n`);
        const natoms = structure.names.length;
        buffer.push(`${natoms}\n`);
        for (let i=0; i<natoms; i++) {
            const [x, y, z] = cell.fractional([structure.x[i], structure.y[i], structure.z[i]])
            buffer.push(`${structure.names[i]} ${i} - ${x} ${y} ${z}\n`);
        }
    }
    return buffer.join('');
}
