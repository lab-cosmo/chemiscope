/**
 * @packageDocumentation
 * @module structure
 */

import * as $3Dmol from '3dmol';
import { Environment, Structure, UserStructure } from '../dataset';
import { computeSecondaryStructure } from './topology';
import { MoleculeViewer } from './viewer';
import { Indexes } from '../indexer';

/** We need to use the same opacity everywhere since 3Dmol uses 1 opacity per model */
export function defaultOpacity(): number {
    // The actual visual appearance seems to depend on browser and platform.
    // This is a value that seems to work well across the board, but it might
    // need to become platform-dependent in case of changes in the browsers.
    return 0.85;
}

/**
 * Add data from the `structure` to the `model`
 *
 * @param model 3Dmol GLModel that will contain structure data
 * @param structure the structure to convert
 */
export function setup3DmolStructure(model: $3Dmol.GLModel, structure: Structure): void {
    if (structure.cell !== undefined) {
        const cell = structure.cell;
        // prettier-ignore
        const matrix = new $3Dmol.Matrix3(
            cell[0], cell[3], cell[6],
            cell[1], cell[4], cell[7],
            cell[2], cell[5], cell[8]
        );
        model.setCrystMatrix(matrix);
    }

    const atoms = [];
    for (let i = 0; i < structure.size; i++) {
        const x = structure.x[i];
        const y = structure.y[i];
        const z = structure.z[i];
        const atom = {
            serial: i,
            elem: structure.elements ? structure.elements[i] : structure.names[i],
            atom: structure.names[i],
            x: x,
            y: y,
            z: z,
            hetflag: true,
            ss: 'c', // initialize to coil by default, the same as in 3Dmol
        } as unknown as $3Dmol.AtomSpec;

        if (structure.resnames !== undefined) {
            atom.resn = structure.resnames[i];
        }

        if (structure.chains !== undefined) {
            atom.chain = structure.chains[i];
        }

        if (structure.resids !== undefined) {
            atom.resi = structure.resids[i];
        }

        if (structure.hetatom !== undefined) {
            atom.hetflag = structure.hetatom[i];
        }

        atoms.push(atom);
    }

    computeSecondaryStructure(atoms, /*hbondCutoff=*/ 3.2);
    model.addAtoms(atoms);
}

export interface LabeledArrow {
    label: $3Dmol.Label;
    arrow: $3Dmol.GLShape;
}

export interface ColorBar {
    /** The numbers scaling the color bar for the minimum, the middle
     * and the maximum values respectively displayed from left to right
     * under the color bar */
    min: $3Dmol.Label;
    mid: $3Dmol.Label;
    max: $3Dmol.Label;
    /** The property used in the color bar,
     * displayed under the number labels */
    property: $3Dmol.Label;
    // The color gradient (i.e. the main element) of the color bar
    gradient: $3Dmol.Label;
}

/** Possible options passed to `MoleculeViewer.load` */
export interface LoadOptions {
    /** Supercell to display (default: [1, 1, 1]) */
    supercell: [number, number, number];
    /** Should preserve we the current camera orientation (default: false) */
    keepOrientation: boolean;
    /** Are we loading a file part of a trajectory (default: false) */
    trajectory: boolean;
    /** List of atom-centered environments in the current structure, potentially
     * undefined if the environnement is not part of the dataset.
     * The `structure` field of `Environment` is ignored */
    environments: (Environment | undefined)[];
    /**
     * Index of the environment to highlight, this is only considered if
     * `environments` is set.
     */
    highlight: number;
}

/// Extension of `Environment` adding the global index of the environment
export interface NumberedEnvironment extends Environment {
    index: number;
}

/**
 * Create a list of environments grouped together by structure.
 *
 *
 * This function returns `undefined` if `environment` is undefined, else it
 * returns a list of list of environments, such as `list[0]` contains all
 * environments in structure 0; `list[33]` all environments in structure 33, etc.
 *
 * @param  structures Expected number of structures
 * @param  environments Full list of environments
 *
 * @return              The list of environments grouped by structure
 */
export function groupByStructure(
    structures: (Structure | UserStructure)[],
    environments?: Environment[]
): NumberedEnvironment[][] | undefined {
    if (environments === undefined) {
        return undefined;
    }

    const result = Array.from({ length: structures.length }).map((_, i) =>
        Array.from({ length: structures[i].size })
    );

    for (let i = 0; i < environments.length; i++) {
        const env = environments[i];
        result[env.structure][env.center] = {
            index: i,
            ...env,
        };
    }

    return result as NumberedEnvironment[][];
}

export interface ViewerGridData {
    /// the viewer itself
    viewer: MoleculeViewer;
    /// color associated with this viewer
    color: string;
    /// set of indexes currently displayed in this viewer
    current: Indexes;
}

/**
 *  Creates a download request from a URI
 * @param uri       URI of the image
 * @param name      Name of the downloaded image
 */
export function downloadURI(uri: string, name: string) {
    const link = document.createElement('a');
    link.download = name;
    link.href = uri;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    link.remove();
}
