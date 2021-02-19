from copy import deepcopy
import mdtraj as md
import numpy as np

__all__ = ['wrap']

def wrap(t, whole_molecules=True, center=None, inplace=False):
    assert isinstance(t, md.Trajectory), "Trajectory must be a MDTraj Trajectory object"
    xyz = t.xyz
    box = t.unitcell_lengths
    res_atoms = [[a.index for a in residue.atoms]
                 for residue in t.top.residues]

    if center is not None:
        assert (len(center) == 3), "Center must be a list of length 3"
        try:
            center = np.array(center, dtype=float)
        except ValueError:
            raise ValueError("Center must be a list of floats")
        xyz = xyz + box[:,None,:] * 0.5 - center[None,None,:]

    res_atoms = np.array(res_atoms, dtype=object)
    if whole_molecules:
        res_cogs = _get_cogs(xyz, res_atoms)
        res_images = _get_image(res_cogs, box)
        images = np.zeros_like(xyz)
        for res_idx, res in enumerate(res_atoms):
                images[:,res,:] = res_images[:,res_idx,:][:,None,:]
        wrapped_xyz = xyz - box[:,None,:] * images

    else:
        images = _get_image(xyz, box)
        wrapped_xyz = xyz - box[:,None,:] * images

    if center is not None:
        wrapped_xyz = wrapped_xyz - box[:,None,:] * 0.5 + center[None,None,:]

    if inplace:
        t.xyz = wrapped_xyz
    else:
        new_traj = deepcopy(t)
        new_traj.xyz = wrapped_xyz

    return new_traj


def _get_cogs(xyz, res_atoms):
    res_cogs = np.zeros((xyz.shape[0], len(res_atoms), 3), dtype=float)
    for i, res in enumerate(res_atoms):
        res_cog = np.mean(xyz[:,res,:], axis=1)
        res_cogs[:,i,:] = res_cog
    return res_cogs

def _get_image(xyz, box):
    images = (xyz // box[:,None,:]).astype(int)
    return images
