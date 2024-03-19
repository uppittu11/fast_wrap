from copy import deepcopy
import mdtraj as md
import numpy as np

__all__ = ['wrap']

def wrap(t, whole_molecules=True, center=None, inplace=False):
    """ Wrap the coordinates in the simulation box for atoms in a trajectory

    Parameters:
    -----------
    t : MDTraj.Trajectory
        The trajectory containing coordinates to be wrapped.
    whole_molecules : bool, optional, default=True
        Set to true if the residues in the trajectory should be kept whole. See
        disclaimer about polymer systems at the top of this script.
    center : np.ndarray, shape=(3), optional, default=None
        The coordinates of the box center. If `center` is `None`, the center is
        assumed to be located at half the box lengths
        (i.e. 0.5 * t.unitcell_lengths).
    inplace : bool, optional, default=False
        Set to true if the input trajectory `t` should be modified. Setting this
        to `True` can be more time and memory efficient as the topology will not
        be copied.

    Returns:
    --------
    new_traj : MDTraj.Trajectory
        The wrapped trajectory containing coordinates to be wrapped. If `inplace`
        is `True`, this will point to the original trajectory `t`.
    """

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
            res = np.array(res, dtype=int)
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
    """ Calculates the center of geometry (cog) for each residue.

    For a given frame,
    cog of a residue = mean of the coordinates of the atoms in residue

    Parameters:
    -----------
    xyz : np.ndarray, shape=(f, n, 3)
        The atom coordinates with `f` frames and `n` atoms.
    res_atoms : np.ndarray, shape=(m,)
        The indices corresponding to each of `m` residues. Note the dtype of
        this array is `object` since the number of atoms in each residue is
        variable

    Returns:
    --------
    res_cogs : np.ndarray, shape=(f, m, 3)
        The residue cog coordinates with `f` frames and `m` residues.
    """

    res_cogs = np.zeros((xyz.shape[0], len(res_atoms), 3), dtype=float)
    for i, res in enumerate(res_atoms):
        res = np.array(res, dtype=int)
        res_cog = np.mean(xyz[:,res,:], axis=1)
        res_cogs[:,i,:] = res_cog
    return res_cogs

def _get_image(xyz, box):
    """ Calculates the image of each supplied coordinate.

    The image measures how many box lengths away from the simulation box a
    particle is.

    For a given frame,
    image = number of box lengths you can fit between the origin and the coordinate

    Parameters:
    -----------
    xyz : np.ndarray, shape=(f, n, 3)
        The atom coordinates with `f` frames and `n` particles.
    box : np.ndarray, shape=(f, 3)
        The box lengths for each frame (output of traj.unitcell_lengths)

    Returns:
    --------
    images : np.ndarray, shape=(f, m, 3)
        The image residue cog coordinates with `f` frames and `m` residues.
    """

    images = (xyz // box[:,None,:]).astype(int)
    return images
