/**
 * @packageDocumentation
 * @module main
 */

import { Arrow, CustomShape, Ellipsoid, Sphere } from './structure/shapes';
import { ShapeParameters } from './structure/shapes';

/** A dataset containing all the data to be displayed. */
export interface Dataset {
    /** metadata for this dataset */
    meta: Metadata;
    /**
     * List of structures in the dataset.
     *
     * The structures can either follow the `Structure` interface, in which
     * case they will be loaded as-defined; or contain any kind of data as a
     * {@link UserStructure}, in which case the {@link ViewersGrid.loadStructure}
     * callback should be set to translate from whatever is stored to a
     * {@link Structure}.
     */
    structures: Structure[] | UserStructure[];
    /**
     * Dictionary containing shape information about the dataset.
     * Each shape is associated with a name that will also be used in the GUI,
     * and the definition of the type and parameters of each shape is given by
     * in the {@link ShapeParameters} dictionary.
     */
    shapes?: { [name: string]: ShapeParameters };
    /**
     * List of properties for the structures (`target == "structure"`), or
     * atom-centered environments in the structures (`target == "atom"`).
     *
     * For structure properties, the `values` array of the property should have
     * the same size as the structure list in `Dataset.structures`.
     *
     * For atomic properties, the `values` array of the property should have
     * the same size as the environments list in `Dataset.environments`.
     */
    properties: { [name: string]: Property };
    /**
     * List of atom-centered environments in the dataset.
     *
     * Currently, the code assumes that every atom is associated with an
     * environment. This may change in the future.
     */
    environments?: Environment[];
    /** Settings for visualization of this dataset */
    settings?: Partial<Settings>;
    /** Parameters of multidimensional properties */
    parameters?: { [name: string]: Parameter };
}

/**
 * Type definition for settings that can be saved with a dataset. THey should be
 * a simple object with string keys, scalar values, array values or nested
 * Settings objects.
 */
// eslint-disable-next-line @typescript-eslint/no-empty-interface
export interface Settings
    extends Record<
        string,
        string | string[] | number | number[] | boolean | boolean[] | Settings | Settings[]
    > {}

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
}

/** A single atomic structure */
export interface Structure {
    /** Number of atoms in the structure */
    size: number;
    /** Names of all atoms in the structure */
    names: string[];
    /**
     * x position (cartesian coordinate) of all atoms in the structure
     *
     * This array should have the same size as {@link Structure.names}, and contain
     * values expressed in Angströms.
     */
    x: number[];
    /**
     * y position (cartesian coordinate) of all atoms in the structure
     *
     * This array should have the same size as {@link Structure.names}, and contain
     * values expressed in Angströms.
     */
    y: number[];
    /**
     * z position (cartesian coordinate) of all atoms in the structure
     *
     * This array should have the same size as {@link Structure.names}, and contain
     * values expressed in Angströms.
     */
    z: number[];
    /**
     * Unit cell of the system, given as `[ax ay az bx by bz cx cy cz]`, where
     * **a**, **b**, and **c** are the unit cell vectors. All values should be
     * expressed in Angströms.
     */
    cell?: number[];
    /**
     * possible shapes to display, multiple groups of shapes with different
     * names are allowed
     */
    shapes?: {
        /**
         * dictionary containing shape data
         */
        [name: string]: ShapeParameters;
    };
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
     * User-defined data which can be turned into a {@link Structure}.
     *
     * {@link ViewersGrid.loadStructure} must be set to be able to load this
     * data.
     */
    data: unknown;
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
     *
     * 2D properties are also supported, with the goal of displaying one curve
     * for each atom/structure. In this case, the first dimension of the 2D
     * array corresponds to number of structures/atoms and the second dimension
     * corresponds to the values for each structure/atom.
     */
    values: string[] | number[] | number[][];
    /** user-facing description of the property */
    description?: string;
    /** unit of the property values */
    units?: string;
    /** parameter name associated to the multidimensional property */
    parameters?: string[];
}

/**
 * An atom-centered environment.
 *
 * Currently, only spherical (i.e. cutoff-based) environments are supported
 */
export interface Environment {
    /** Index of the related structure in {@link Dataset.structures} */
    structure: number;
    /** Index of the central atom in the structure, 0-based */
    center: number;
    /** Spherical cutoff radius, expressed in Angströms */
    cutoff: number;
}

/**
 * Parameters are used to describe multidimensional properties.
 *
 * The first dimension of a multidimensional property always corresponds to
 * structure/atoms in the dataset, and additional dimensions should be described
 * with a `Parameter`.
 */
export interface Parameter {
    /** name of the parameter */
    name?: string;
    /** units of the elements in the property array */
    units?: string;
    /** values of the parameter */
    values: number[];
}

/** Arbitrary javascript object, to be validated */
export type JsObject = Record<string, unknown>;

/** @hidden
 * Check that the given object, potentially comming from javascript, has all
 * required properties to be a dataset.
 */
export function validateDataset(o: JsObject): void {
    if (typeof o !== 'object') {
        throw Error('the dataset must be a JavaScript object');
    }

    if (!('meta' in o)) {
        throw Error('missing "meta" key in the dataset');
    } else if (!(typeof o.meta === 'object' && o.meta !== null)) {
        throw Error('"meta" must be an object in the dataset');
    }
    checkMetadata(o.meta as JsObject);

    if (!('structures' in o)) {
        throw Error('missing "structures" key in the dataset');
    } else if (!Array.isArray(o.structures)) {
        throw Error('"structures" must be an array in the dataset');
    }
    const [structureCount, atomsCount] = checkStructures(o.structures as JsObject[]);

    let envCount = atomsCount;
    if ('environments' in o) {
        if (!Array.isArray(o.environments)) {
            throw Error('"environments" must be an array in the dataset');
        }

        envCount = o.environments.length;
        checkEnvironments(
            o.environments as JsObject[],
            o.structures as (Structure | UserStructure)[]
        );
    }

    if ('shapes' in o) {
        const check_shape = checkShapes(
            o.shapes as Record<string, JsObject>,
            structureCount,
            envCount
        );
        if (check_shape !== '') {
            throw Error('Error checking shape definitions: ' + check_shape);
        }
        const check_assign = assignShapes(
            o.shapes as { [name: string]: ShapeParameters },
            o.structures as Structure[]
        );
        if (check_assign !== '') {
            throw Error('Error assigning shapes to structures: ' + check_assign);
        }
    }

    if (!('properties' in o)) {
        throw Error('missing "properties" key in then dataset');
    } else if (!(typeof o.properties === 'object' && o.properties !== null)) {
        throw Error('"properties" must be an object in the dataset');
    }

    if ('parameters' in o) {
        if (!(typeof o.parameters === 'object' && o.parameters !== null)) {
            throw Error('"parameters" must be an object in the dataset');
        }
    }

    checkProperties(
        o.properties as Record<string, JsObject>,
        structureCount,
        envCount,
        o.parameters as Record<string, JsObject> | undefined
    );
}

function checkMetadata(o: JsObject) {
    if (!('name' in o)) {
        throw Error('missing "meta.name" key in the dataset');
    } else if (typeof o.name !== 'string') {
        throw Error('"meta.name" must be a string in the dataset');
    }

    if ('description' in o && typeof o.description !== 'string') {
        throw Error('"meta.description" should be a string in the dataset');
    }

    if ('authors' in o) {
        if (!Array.isArray(o.authors)) {
            throw Error('"meta.authors" must be an array in the dataset');
        }

        for (const a of o.authors) {
            if (typeof a !== 'string') {
                throw Error('"meta.authors" must be an array of strings in the dataset');
            }
        }
    }

    if ('references' in o) {
        if (!Array.isArray(o.references)) {
            throw Error('"meta.references" must be an array in the dataset');
        }

        for (const a of o.references) {
            if (typeof a !== 'string') {
                throw Error('"meta.references" must be an array of strings in the dataset');
            }
        }
    }
}

function checkStructures(o: JsObject[]): [number, number] {
    let atomsCount = 0;
    for (let i = 0; i < o.length; i++) {
        const structure = o[i];
        if (
            !(
                'size' in structure &&
                typeof structure.size === 'number' &&
                isPositiveInteger(structure.size)
            )
        ) {
            throw Error(`missing 'size' for structure ${i}`);
        }
        atomsCount += structure.size;

        if ('data' in structure) {
            // user-specified structure, nothing to do
        } else {
            const message = checkStructure(structure);
            if (message !== '') {
                throw Error(`error in structure ${i}: ${message}`);
            }
        }
    }

    return [o.length, atomsCount];
}

function checkShapes(
    shapes: Record<string, JsObject>,
    structureCount: number,
    envCount: number
): string {
    if (typeof shapes !== 'object' || shapes === null) {
        return "'shapes' must be an object";
    }

    for (const [key, o_shape] of Object.entries(shapes as object)) {
        const shape = o_shape as ShapeParameters;
        if (!('kind' in shape)) {
            return `missing "kind" in shape ${key}`;
        }

        if (typeof shape.kind !== 'string') {
            return `shapes 'kind' must be a string for shape ${key}`;
        }

        if (!('parameters' in shape)) {
            return `missing "parameters" in shape ${key}`;
        }

        const parameters = shape.parameters;

        if ('structure' in parameters) {
            const s_parameters = parameters.structure;
            if (!Array.isArray(s_parameters)) {
                return `'structure' parameters should be an array in shape ${key}`;
            }

            if (s_parameters.length !== structureCount) {
                return `'structure' parameters in shape ${key} contain ${s_parameters.length} entries, but there are ${structureCount} structures.`;
            }
        }

        if ('atom' in parameters) {
            const a_parameters = parameters.atom;
            if (!Array.isArray(a_parameters)) {
                return `'atom' parameters should be an array in shape ${key}`;
            }

            if (a_parameters.length !== envCount) {
                return `'atom' parameters in shape ${key} contain ${a_parameters.length} entries, but there are ${envCount} environments.`;
            }
        }
    }

    return '';
}

function validateShape(kind: string, parameters: Record<string, unknown>): string {
    if (kind === 'sphere') {
        return Sphere.validateParameters(parameters);
    } else if (kind === 'ellipsoid') {
        return Ellipsoid.validateParameters(parameters);
    } else if (kind === 'arrow') {
        return Arrow.validateParameters(parameters);
    } else if (kind === 'custom') {
        return CustomShape.validateParameters(parameters);
    }
    return '';
}

function assignShapes(
    shapes: { [name: string]: ShapeParameters },
    structures: Structure[]
): string {
    // creates shapes associated with actual structures by combining all the information given in the definition
    let atomsCount = 0;
    for (let i_structure = 0; i_structure < structures.length; i_structure++) {
        const structure = structures[i_structure];
        structure.shapes = {};
        for (const [name, shape] of Object.entries(shapes)) {
            const parameters = {
                global: shape.parameters.global,
                structure: shape.parameters.structure,
                atom: shape.parameters.atom,
            };

            let full_parameters = shape.parameters.global;
            if (parameters.structure) {
                parameters.structure = [parameters.structure[i_structure]];
                full_parameters = { ...full_parameters, ...parameters.structure[0] };
            }

            if (parameters.atom) {
                parameters.atom = parameters.atom.slice(atomsCount, atomsCount + structure.size);

                for (const atom of parameters.atom) {
                    const atom_parameters = { ...full_parameters, ...atom };
                    const check = validateShape(shape.kind, atom_parameters);
                    if (check !== '') {
                        return `Validation error for an atom in shape ${name}: ${check}`;
                    }
                }
            } else {
                const check = validateShape(shape.kind, full_parameters);
                if (check !== '') {
                    return `Validation error for a structure in shape ${name}: ${check}`;
                }
            }

            structure.shapes[name] = {
                kind: shape.kind,
                parameters: parameters,
            };
        }

        atomsCount += structure.size;
    }

    return ''; // success!
}

/**
 * Check that the given object is a structure. Return a string describing the
 * issue with `s` if any, or the empty string if `s` looks like a valid
 * structure.
 */
export function checkStructure(s: JsObject): string {
    if (typeof s !== 'object' || s === null) {
        throw Error('the structure must be a JavaScript object');
    }

    if (!('size' in s && typeof s.size === 'number' && isPositiveInteger(s.size))) {
        return 'missing "size"';
    }

    for (const key of ['names', 'x', 'y', 'z']) {
        if (!(key in s)) {
            return `missing "${key}"`;
        }
        const array = s[key];
        if (!Array.isArray(array)) {
            return `"${key}" must be an array`;
        }

        if (s.size > 0 && array.length !== s.size) {
            return `wrong size for "${key}", expected ${s.size}, got ${array.length}`;
        }
    }

    if ('cell' in s) {
        if (!(Array.isArray(s.cell) && s.cell.length === 9)) {
            return '"cell" must be an array of size 9';
        }
    }

    return '';
}

function checkProperties(
    properties: Record<string, JsObject>,
    structureCount: number,
    envCount: number,
    parameters?: Record<string, JsObject> | undefined
) {
    for (const key in properties) {
        const property = properties[key];

        if (!('target' in property && typeof property.target === 'string')) {
            Error(`'properties['${key}'].target' should be a string`);
        }

        if (property.target !== 'atom' && property.target !== 'structure') {
            throw Error(`'properties['${key}'].target', should be 'atom' | 'structure'`);
        }

        if (!('values' in property && Array.isArray(property.values))) {
            throw Error(`'properties["${key}"].values' should be an array`);
        }

        // check size if possible
        let expected = 0;
        if (property.target === 'atom') {
            expected = envCount;
        } else if (property.target === 'structure') {
            expected = structureCount;
        }

        if (expected > 0 && property.values.length !== expected) {
            throw Error(
                `wrong size for 'properties['${key}'].values': expected ${expected}, got ${property.values.length}`
            );
        }

        const initial = typeof property.values[0];
        if (
            initial !== 'string' &&
            initial !== 'number' &&
            !isMultidimensional(property.values as number[][])
        ) {
            throw Error(
                `'properties['${key}'].values' should contain string or number or an array of numbers`
            );
        }

        for (const value of property.values) {
            if (typeof value !== initial) {
                throw Error(`'properties['${key}'].values' should be of a single type`);
            }
        }

        // few checks on multidimensional properties
        if (isMultidimensional(property.values as number[][])) {
            // check if parameters exists
            if (!parameters) {
                throw Error(
                    `'parameters' should be provided for multidimensional properties '${key}'`
                );
            }
            // check if parameter keyword exists and has the right format
            const propertyParameters = property.parameters as string[];
            if (!isArrayString(propertyParameters)) {
                throw Error(`'properties['${key}'].parameters' should be an array of strings`);
            }
            // check if the length of parameters is 1 TODO: remove when support for multiple parameters is ready
            if (propertyParameters.length !== 1) {
                throw Error(`'properties['${key}'].parameters' should contain a single parameter`);
            }
            // check if parameters of the property exists in the parameters
            //for (const value of propertyParameters) {
            for (const value of propertyParameters) {
                if (!(value in parameters)) {
                    throw Error(
                        `parameter '${value}' of 'properties['${key}']' does not appear in the list provided parameters`
                    );
                }
            }
            // check if the length of the first array matches the length of the parameters
            const initialValues = property.values[0] as number[];
            for (const value of propertyParameters) {
                const parameterValues = parameters[value].values as number[];
                if (initialValues.length !== parameterValues.length) {
                    throw Error(
                        `'properties['${key}'].values' and 'parameters['${value}'].values' should have the same length`
                    );
                }
            }

            // check if all the multidimensial array elements have the same length
            if (!isConsistent2DArray(property.values as number[][])) {
                throw Error(
                    `'properties['${key}].values' should contain arrays of the same length`
                );
            }
        }

        // check that units & description are valid
        if ('description' in property && typeof property.description !== 'string') {
            throw Error(`'properties['${key}'].description' should contain a string`);
        }

        if ('units' in property && typeof property.units !== 'string') {
            throw Error(`'properties['${key}'].units' should contain a string`);
        }
    }
}

function checkEnvironments(o: JsObject[], structures: (Structure | UserStructure)[]) {
    for (let i = 0; i < o.length; i++) {
        const env = o[i];

        if (!('structure' in env && typeof env.structure === 'number')) {
            throw Error(`missing 'structure' for environment ${i}`);
        }

        if (!isPositiveInteger(env.structure) || env.structure >= structures.length) {
            throw Error(
                `out of bounds 'structure' for environment ${i}: index is \
                ${env.structure}, we have ${structures.length} structures`
            );
        }

        if (!('center' in env && typeof env.center === 'number')) {
            throw Error(`missing 'center' for environment ${i}`);
        }

        const size = structures[env.structure].size;
        if (!isPositiveInteger(env.center) || env.center >= size) {
            throw Error(
                `out of bounds 'center' for environment ${i}: index is \
                ${env.center}, we have ${size} atoms in structure ${env.structure}`
            );
        }

        if (!('cutoff' in env && typeof env.cutoff === 'number')) {
            throw Error(`missing 'cutoff' for environment ${i}`);
        }
    }
}

function isPositiveInteger(number: number): boolean {
    return Number.isInteger(number) && number >= 0;
}

function isMultidimensional(array: number[][]): boolean {
    // check if an array is 2D
    let result = true;
    for (const value of array) {
        result = Array.isArray(value) && result;
    }
    return result;
}

function isConsistent2DArray(array: number[][]): boolean {
    // check if the elements of 2D array have the same length
    const initial = array[0];
    let result = true;
    for (const value of array) {
        result = value.length === initial.length && result;
    }
    return result;
}

function isArrayString(array: unknown): boolean {
    if (Array.isArray(array)) {
        let result = true;
        for (const value of array) {
            result = typeof value === 'string' && result;
        }
        return result;
    } else {
        return false;
    }
}
