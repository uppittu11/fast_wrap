import argparse

import mdtraj as md

import fast_wrap

def wrap_func():
    parser = argparse.ArgumentParser(description="Wrap a trajectory.")
    parser.add_argument("-f", dest="traj",
            help="unwrapped trajectory file; must be a MDTraj supported file format.")
    parser.add_argument("-c", dest="top", default=None,
            help="topology file; must be a MDTraj supported file format.")
    parser.add_argument("-o", dest="output",
            help="output file; must be a MDTraj supported file format.")
    parser.add_argument("--nowhole", dest="whole", action="store_false", default=True,
            help="[Don't] keep molecules whole.")

    args = parser.parse_args()
    if args.top:
        traj = md.load(args.traj, top=args.top)
    else:
        traj = md.load(args.traj)
    wrapped = fast_wrap.wrap(traj, whole_molecules=args.whole)
    wrapped.save(args.output)
