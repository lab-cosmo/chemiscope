/**
 * @packageDocumentation
 * @module main
 */

/** A dataset containing all the data to be displayed. */
export interface Dataset {
    /** metadata for this dataset */
    meta: Metadata;
    /**
     * List of structures in the dataset.
     *
     * The structures can either follow the `Structure` interface, in which
     * case they will be loaded as-defined; or contain any kind of data as a
     * [[UserStructure]], in which case the [[StructureViewer.loadStructure]]
     * callback should be set to translate from whatever is stored to a
     * [[Structure]].
     */
    structures: Structure[] | UserStructure[];
    /**
     * List of properties for the structures (`target == "structure"`), or
     * atom-centered environments in the structures (`target == "atom"`).
     *
     * For structure properties, the `values` array of the property should have
     * the same size as the [[Dataset.structures|structure list]].
     *
     * For atomic properties, the `values` array of the property should have
     * the same size as the [[Dataset.environments|environments list]].
     */
    properties: {
        [name: string]: Property;
    };
    /**
     * List of atom-centered environments in the dataset.
     *
     * Currently, the code assumes that every atom is associated with an
     * environment. This may change in the future.
     */
    environments?: Environment[];
}

/** Various metadata associated with a dataset */
export interface Metadata {
    /** dataset name */
    name: string;
    /** authors of the dataset */
    authors?: string[];
    /** academic references related to this dataset */
    references?: string[];
    /** description of the dataset */
    description?: string;
    /** style settings **/
    presets?: string; // probably better type is called for 
}

/** A single atomic structure */
export interface Structure {
    /** Number of atoms in the structure */
    size: number;
    /** Names of all atoms in the structure */
    names: number[];
    /**
     * x position (cartesian coordinate) of all atoms in the structure
     *
     * This array should have the same size as [[Structure.names]], and contain
     * values expressed in Angströms.
     */
    x: number[];
    /**
     * y position (cartesian coordinate) of all atoms in the structure
     *
     * This array should have the same size as [[Structure.names]], and contain
     * values expressed in Angströms.
     */
    y: number[];
    /**
     * z position (cartesian coordinate) of all atoms in the structure
     *
     * This array should have the same size as [[Structure.names]], and contain
     * values expressed in Angströms.
     */
    z: number[];
    /**
     * Unit cell of the system, given as `[ax ay az bx by bz cx cy cz]`, where
     * **a**, **b**, and **c** are the unit cell vectors. All values should be
     * expressed in Angströms.
     */
    cell?: number[];
}

/**
 * User-defined data to allow dynamic loading of the structures.
 *
 * The main use-case of this is making the initial loading time of chemiscope
 * faster by loading structure on-demand, from files, a database or even a
 * javascript program.
 */
export interface UserStructure {
    /** Number of atoms in the structure */
    size: number;
    /**
     * User-defined data which can be turned into a [[Structure]].
     *
     * [[StructureViewer.loadStructure]] must be set to be able to load this
     * data.
     */
    data: any;
}

/** Possible types of properties: full structure property, or atomic property */
export type Target = 'structure' | 'atom';

/**
 * A single property in a dataset.
 *
 * Properties can be physical properties (energy, number of atoms, density,
 * *etc.*); values associated with the structure (such as SOAP vectors),
 * projected onto a lower dimensionality sub-space (through PCA or any other
 * algorithm); or any other value associated with every structure/environment in
 * the dataset.
 */
export interface Property {
    /** is this property associated with a full structure or a single atom? */
    target: Target;
    /**
     * values of the property
     *
     * string values should represent classification results (category "A", "B"
     * or "C"); and numeric values should be use for everything else.
     */
    values: string[] | number[];
}

/**
 * An atom-centered environment.
 *
 * Currently, only spherical (i.e. cutoff-based) environments are supported
 */
export interface Environment {
    /** Index of the related structure in [[Dataset.structures]] */
    structure: number;
    /** Index of the central atom in the structure, 0-based */
    center: number;
    /** Spherical cutoff radius, expressed in Angströms */
    cutoff: number;
}

/** @hidden
 * Check that the given object, potentially comming from javascript, has all
 * required properties to be a dataset.
 */
export function checkDataset(o: any) {
    if (!('meta' in o && typeof o.meta === 'object')) {
        throw Error('missing "meta" key in dataset');
    }
    checkMetadata(o.meta);

    if (!('structures' in o && o.structures.length !== undefined)) {
        throw Error('missing "structures" key in dataset');
    }
    const [structureCount, envCount] = checkStructures(o.structures);

    if (!('properties' in o && typeof o.properties === 'object')) {
        throw Error('missing "properties" key in dataset');
    }
    checkProperties(o.properties, structureCount, envCount);

    if ('environments' in o) {
        if (typeof o.environments.length === undefined) {
            throw Error('"environments" must be an array in dataset');
        }

        if (o.environments.length !== envCount) {
            throw Error(`expected ${envCount} environments, got ${o.environments.length} instead`);
        }
        checkEnvironments(o.environments, o.structures);
    }
}

function checkMetadata(o: any) {
    if (!('name' in o && typeof o.name === 'string')) {
        throw Error('missing "meta.name" key in dataset');
    }

    if ('description' in o && typeof o.description !== 'string') {
        throw Error('"meta.description" should be a string in dataset');
    }

    if ('authors' in o) {
        if (typeof o.authors.length === undefined) {
            throw Error('"meta.authors" must be an array in dataset');
        }

        for (const a of o.authors) {
            if (typeof a !== 'string') {
                throw Error('"meta.authors" must be an array of strings in dataset');
            }
        }
    }

    if ('references' in o) {
        if (typeof o.references.length === undefined) {
            throw Error('"meta.references" must be an array in dataset');
        }

        for (const a of o.references) {
            if (typeof a !== 'string') {
                throw Error('"meta.references" must be an array of strings in dataset');
            }
        }
    }
}

export function checkStructures(o: any): [number, number] {
    let envCount = 0;
    for (let i = 0; i < o.length; i++) {
        const structure = o[i];
        if (!('size' in structure && typeof structure.size === 'number' && isPositiveInteger(structure.size))) {
            throw Error(`missing 'size' for structure ${i}`);
        }
        envCount += structure.size;

        if ('names' in structure && 'x' in structure && 'y' in structure && 'z' in structure) {
            for (const key of ['names', 'x', 'y', 'z']) {
                if (!(key in structure && structure[key].length !== undefined)) {
                    throw Error(`missing '${name}' for structure ${i}`);
                }
                if (structure[key].length !== structure.size) {
                    throw Error(`wrong size for '${name}' in structure ${i}, expected ${structure.size}, got ${structure[name].length}`);
                }
            }

            if ('cell' in structure) {
                if (structure.cell.length !== 9) {
                    throw Error(`'cell' must be an array of size 9 for structure ${i}`);
                }
            }
        } else if ('data' in structure) {
            // nothing to check
        } else {
            throw Error(`structure ${i} must contains either 'names'/'x'/'y'/'z' or 'data'`);

        }
    }

    return [o.length, envCount];
}

function checkProperties(o: any, structureCount: number, envCount: number) {
    for (const key in o.properties) {
        const property = o.properties[key];

        if (!('target' in property && typeof property.target === 'string')) {
            Error(`'properties['${key}'].target' should be a string`);
        }

        if (property.target !== 'atom' && property.target !== 'structure') {
            throw Error(`'properties['${key}'].target', should be 'atom' | 'structure'`);
        }

        if (!('values' in property && property.values.length !== undefined)) {
            Error(`'properties['${key}'].values' should be an array`);
        }

        // check size is possible
        let expected;
        if (property.target === 'atom') {
            expected = envCount;
        } else if (property.target === 'structure') {
            expected = structureCount;
        } else {
            throw Error(`invalid property target: ${property.target}`);
        }

        if (property.values.length !== expected) {
            throw Error(`wrong size for 'properties['${key}'].values': expected ${expected}, got ${property.values.length}`);
        }

        const initial = typeof property.values[0];
        if (initial !== 'string' && initial !== 'number') {
            throw Error(`'properties['${key}'].values' should contain string or number`);
        }

        for (const value of property.values) {
            if (typeof value !== initial) {
                throw Error(`'properties['${key}'].values' should be of a single type`);
            }
        }
    }
}

function checkEnvironments(o: any, structures: Structure[]) {
    for (let i = 0; i < o.length; i++) {
        const env = o[i];

        if (!('structure' in env && typeof env.structure === 'number')) {
            throw Error(`missing 'structure' for environment ${i}`);
        }

        if (!isPositiveInteger(env.structure) || env.structure >= structures.length) {
            throw Error(
                `out of bounds 'structure' for environment ${i}: index is \
                ${env.structure}, we have ${structures.length} structures`,
            );
        }

        if (!('center' in env && typeof env.center === 'number')) {
            throw Error(`missing 'center' for environment ${i}`);
        }

        const size = structures[env.structure].size;
        if (!isPositiveInteger(env.center) || env.center >= size) {
            throw Error(
                `out of bounds 'center' for environment ${i}: index is \
                ${env.center}, we have ${size} atoms in structure ${env.structure}`,
            );
        }

        if (!('cutoff' in env && typeof env.cutoff === 'number')) {
            throw Error(`missing 'cutoff' for environment ${i}`);
        }
    }
}

/** @hidden */
export function isPositiveInteger(number: number): boolean {
    return Number.isInteger(number) && number >= 0;
}
