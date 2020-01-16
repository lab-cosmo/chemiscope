import * as linalg from './linalg'

import {Structure} from '../dataset';

class UnitCell {
    private matrix: linalg.Matrix;
    private inverse: linalg.Matrix;

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

    /// convert from cartesian to fractional coordinates
    public fractional(cartesian: linalg.Vector3D): linalg.Vector3D {
        return linalg.dot(this.inverse, cartesian);
    }

    /// convert from cartesian to fractional coordinates
    public cartesian(fractional: linalg.Vector3D): linalg.Vector3D {
        return linalg.dot(this.matrix, fractional);
    }

    /// Get the lenghts of the unitcell
    public lengths(): linalg.Vector3D {
        return [
            linalg.norm(this.matrix[0]),
            linalg.norm(this.matrix[1]),
            linalg.norm(this.matrix[2]),
        ];
    }

    /// Get the angles of the unitcell, in degrees
    public angles(): linalg.Vector3D {
        return [
            linalg.angle(this.matrix[1], this.matrix[2]),
            linalg.angle(this.matrix[0], this.matrix[2]),
            linalg.angle(this.matrix[0], this.matrix[1]),
        ];
    }
}

/// Convert a structure to a format that JSMol can read: XYZ if there is no cell
/// data, and BCS if there is cell data
export function structure2JSmol(s: Structure): string {
    let buffer = new Array();
    if (s.cell === undefined) {
        // use XYZ format
        const natoms = s.names.length;
        buffer.push(`${natoms}\n\n`);
        for (let i=0; i<s.names.length; i++) {
            buffer.push(`${s.names[i]} ${s.x[i]} ${s.y[i]} ${s.z[i]}\n`);
        }
    } else {
        // use BCS format
        const cell = new UnitCell(s.cell);
        buffer.push(`1\n`);
        const [a, b, c] = cell.lengths();
        const [alpha, beta, gamma] = cell.angles();
        buffer.push(`${a} ${b} ${c} ${alpha} ${beta} ${gamma}\n`);
        const natoms = s.names.length;
        buffer.push(`${natoms}\n`);
        for (let i=0; i<s.names.length; i++) {
            const [x, y, z] = cell.fractional([s.x[i], s.y[i], s.z[i]])
            buffer.push(`${s.names[i]} ${i} - ${x} ${y} ${z}\n`);
        }
    }
    return buffer.join('');
}
