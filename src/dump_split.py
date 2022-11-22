#!/usr/bin/env python
import numpy as np
import argparse, os

def Split_Main(Iargs):
    """Splits a LAMMPS dumpfile by frame

    This takes a LAMMPS custom dump file and splits it into separate frames. This can significantly
    speed up the rerun command and lead to signficant savings (factor of orders of magnitude in wall-clock time).

    Writes files to dump_dir in the current directory.

    Args:
        Iargs (Argparse object): Input arguments read through argparse.

    """

    if not os.path.exists("./dump_dir"):
        os.makedirs("./dump_dir")
    with open(Iargs.trj,'r') as f:
        for i in range(Iargs.fcount):
            lines = []
            for j in range(Iargs.nper):
                lines.append(f.readline())
            fr = int(lines[1].strip())
            outf = open('./dump_dir/dump.lammpsdump-%d'%fr,'w')
            for line in lines:
                outf.write("%s\n" % line.strip())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate Radial Distribution Function')
    parser.add_argument('-trj',     default="dump.lammpsdump",type=str, help = 'Dump file')
    parser.add_argument('-fcount',  default=10000,          type=int,   help = 'Number of fps')
    parser.add_argument('-nper',    default=3465,           type=int,   help = 'Number of lines per frame')

    Iargs = parser.parse_args()

    Split_Main(Iargs)