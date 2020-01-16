/**
 * @packageDocumentation
 * @module main
 */

interface JSmolEnvironment {
    cutoff: number;
}

/// All the data needed to create a Vizualizer
export interface Dataset {
    /// Dataset metadata
    meta: Metadata;
    /// List of structures in the dataset
    structures: Structure[];
    /// List of properties for the structure
    properties: {
        [name: string]: Property;
    };
    /// List of atom-centered environments
    environments?: Environment[];
}

/// Metadata that can be associated with a dataset
interface Metadata {
    /// The dataset name
    name: string;
    // TODO: add more here: author, journal, ...
}

export interface Structure {
    names: number[],
    x: number[],
    y: number[],
    z: number[],
    cell?: number[],
}

export type Target = 'structure' | 'atom';

/// A single property of a dataset
export interface Property {
    /// Property target: are we considering atomic properties or global
    /// properties
    target: Target;
    /// Property value: string values should represent classification results
    /// numeric values can be use for everything else.
    values: string[] | number[];
}

/// An atom-centered environments
export interface Environment extends JSmolEnvironment {
    /// Index of the related structure in Dataset.structure
    structure: number;
    /// Index of the central atom in the structure, 0-based
    center: number;
}

/// Check that the given object, potentially comming from javascript,
/// has all reaquired properties to be a dataset
export function checkDataset(o: any) {
    if (!('meta' in o && typeof o['meta'] === 'object')) {
        throw Error("missing 'meta' in dataset");
    }
    checkMetadata(o['meta']);

    if (!('structures' in o && o['structures'].length !== undefined)) {
        throw Error("missing 'structures' in dataset");
    }
    const [n_structures, n_atoms] = checkStructure(o['structures']);

    if (!('properties' in o && typeof o['properties'] === 'object')) {
        throw Error("missing 'properties' in dataset");
    }
    checkProperties(o['properties'], n_structures, n_atoms);

    if ('environments' in o) {
        if (typeof o['environments'].length === undefined) {
            throw Error("'environments' must be an array in dataset");
        }

        if (o['environments'].length !== n_atoms) {
            throw Error(`epxpected ${n_atoms} environments, got ${o['environments'].length} instead`);
        }
        checkEnvironments(o['environments'], o['structures']);
    }
}

function checkMetadata(o: any) {
    if (!('name' in o && typeof o['name'] === 'string')) {
        throw Error("missing 'meta.name' in dataset");
    }
}

function checkStructure(o: any): [number, number] {
    let n_atoms = 0;
    for (let i=0; i<o.length; i++) {
        const structure = o[i];
        if (!('names' in structure && structure['names'].length !== undefined)) {
            throw Error(`missing 'names' for structure ${i}`);
        }
        const size = structure['names'].length;
        n_atoms += size;

        for (const name of ['x', 'y', 'z']) {
            if (!(name in structure && structure[name].length !== undefined)) {
                throw Error(`missing '${name}' for structure ${i}`);
            }
            if (structure['x'].length !== size) {
                throw Error(`wrong size for '${name}' in structure ${i}, expected ${size}, got ${structure[name].length}`);
            }
        }

        if ('cell' in structure) {
            if (structure['cell'].length !== 9) {
                throw Error(`'cell' must be an array of size 9 for structure ${i}`);
            }
        }
    }
    return [o.length, n_atoms];
}

function checkProperties(o: any, n_structures: number, n_atoms: number) {
    for (const key in o['properties']) {
        const property = o['properties'][key];

        if (!('target' in property && typeof property['target'] === 'string')) {
            Error(`'properties['${key}'].target' should be a string`);
        }

        if (property.target !== 'atom' && property.target !== 'structure') {
            throw Error(`'properties['${key}'].target', should be 'atom' | 'structure'`);
        }

        if (!('values' in property && property.values.length !== undefined)) {
            Error(`'properties['${key}'].values' should be an array`);
        }

        if (property.target == 'atom' && property.values.length != n_atoms) {
            throw Error(`wrong size for 'properties['${key}'].values': expected ${n_atoms}, got ${property.values.length}`);
        } else if (property.target == 'structure' && property.values.length != n_structures) {
            throw Error(`wrong size for 'properties['${key}'].values': expected ${n_structures}, got ${property.values.length}`);
        }

        const initial = typeof property.values[0]
        if (initial != 'string' && initial != 'number') {
            throw Error(`'properties['${key}'].values' should contain string or number`);
        }

        for (let i=0; i<property.values.length; i++) {
            if (typeof property.values[i] !== initial) {
                throw Error(`'properties['${key}'].values' should be of a single type`);
            }
        }
    }
}

function checkEnvironments(o: any, structures: Structure[]) {
    for (let i=0; i<o.length; i++) {
        const env = o[i];

        if (!('structure' in env && typeof env['structure'] === 'number')) {
            throw Error(`missing 'structure' for environment ${i}`);
        }

        if (!isPositiveInteger(env['structure']) || env['structure'] >= structures.length) {
            throw Error(`out of bounds 'structure' for environment ${i}: index is ${env['structure']}, we have ${structures.length} structures`);
        }

        const n_atoms = structures[env['structure']].names.length;
        if (!('center' in env && typeof env['center'] === 'number')) {
            throw Error(`missing 'center' for environment ${i}`);
        }

        if (!isPositiveInteger(env['center']) || env['center'] >= n_atoms) {
            throw Error(`out of bounds 'center' for environment ${i}: index is ${env['center']}, we have ${n_atoms} atoms`);
        }

        if (!('cutoff' in env && typeof env['cutoff'] === "number")) {
            throw Error(`missing 'cutoff' for environment ${i}`);
        }
    }
}

function isPositiveInteger(number: number): boolean {
    return number !== Infinity && Math.floor(number) === number && number >= 0;
}
