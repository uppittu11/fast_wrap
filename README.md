# fast_wrap
A fast script to wrap atoms or molecules into a simulation box using MDTraj as a backend.

# Usage
Use the bash command `wrap_traj`

## Command line arguments for `wrap_traj`
`-f` unwrapped trajectory file (must be a MDTraj supported file format)
`-c` topology file (must be a MDTraj supported file format)
`-o` output file (must be a MDTraj supported file format)
`--nowhole` [don't] keep molecules whole
`--center` center the box at the origin
