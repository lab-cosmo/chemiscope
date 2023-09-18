"""
A simple demonstration of the construction of a PCA map based on SOAP
features computed by librascal, and exported as a chemiscope json.
"""

import numpy as np

import ase.io
import urllib.request

from rascal.representations import SphericalInvariants
from rascal.neighbourlist.structure_manager import mask_center_atoms_by_species
from sklearn.decomposition import PCA

import chemiscope

#  fetch structures from a librascal examples repo
url = "https://raw.githubusercontent.com/cosmo-epfl/librascal-example-data/833b4336a7daf471e16993158322b3ea807b9d3f/inputs/molecule_conformers_dftb.xyz"
structures_fn, _ = urllib.request.urlretrieve(url)

frames = ase.io.read(structures_fn, "::10")

# Only use Carbon and Oxygen atoms as active centers.
# `mask_center_atoms_by_species` add data inside the frames to indicate which
# atoms should be considered by librascal.
active_atoms_count = []
for frame in frames:
    mask_center_atoms_by_species(frame, [6, 8])

    n_active_atoms = len(np.where(frame.numbers == 6)[0])
    n_active_atoms += len(np.where(frame.numbers == 8)[0])
    active_atoms_count.append(n_active_atoms)

# compute SOAP features using librascal
SOAP_HYPERPARAMETERS = {
    "soap_type": "PowerSpectrum",
    "interaction_cutoff": 3,
    "max_radial": 8,
    "max_angular": 6,
    "gaussian_sigma_constant": 0.3,
    "gaussian_sigma_type": "Constant",
    "cutoff_smooth_width": 0.5,
    "radial_basis": "GTO",
}

soap = SphericalInvariants(**SOAP_HYPERPARAMETERS)
atom_features = soap.transform(frames).get_features(soap)

# structure features are just the mean over the environments in each structure
structure_features = np.zeros((len(frames), atom_features.shape[1]))
atom_idx_start = 0
for i, (frame, n_active_atoms) in enumerate(zip(frames, active_atoms_count)):
    atom_idx_stop = atom_idx_start + n_active_atoms
    structure_features[i] = atom_features[atom_idx_start:atom_idx_stop].mean(axis=0)
    atom_idx_start = atom_idx_stop

# The most simple dimensionality reduction possible, using Principal Components Analysis
atom_pca = PCA(n_components=2).fit_transform(atom_features)
structure_pca = PCA(n_components=2).fit_transform(structure_features)

properties = {
    "structure PCA": structure_pca,
    "atom PCA": atom_pca,
}

# pass the list of included environments/active centers to chemiscope
# `librascal_atomic_environments` uses data created by
# `mask_center_atoms_by_species` to get the list of active centers
environments = chemiscope.librascal_atomic_environments(
    frames, cutoff=SOAP_HYPERPARAMETERS["interaction_cutoff"]
)

chemiscope.write_input(
    path="C3H5OH-chemiscope.json.gz",
    meta={"name": "C3H5OH"},
    frames=frames,
    properties=properties,
    environments=environments,
)
