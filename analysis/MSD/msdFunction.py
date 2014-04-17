#!/usr/bin/env python

import MDAnalysis
import MDAnalysis.KDTree.NeighborSearch as NS
import MDAnalysis.analysis.align as AL
import numpy as np

def msd_win(lenH, startframe, stopframe):
    msd=0
    for i in range(0, lenH):
        msd+=1.0/3.0*((stopframe[i][0]-startframe[i][0])**2+(stopframe[i][1]-startframe[i][1])**2+(stopframe[i][2]-startframe[i][2])**2)
    msd/=float(lenH)
    return msd

def msd_block(selection, u1, u2, begin, timeframe, blocklen):
    msd_bl=0
    tmpbegin=begin
    for i in selection:
        begin=tmpbegin
        segment=u2.selectAtoms("segid %s" %(i))
        ref=u1.selectAtoms("segid %s" %(i))
        hy=segment.selectAtoms("name H*")
        car=segment.selectAtoms("name C*")
        ns_hy=NS.AtomNeighborSearch(hy)
        H_atoms=ns_hy.search_list(car, 1.5)
        lenH=len(H_atoms)
        for ts in u2.trajectory[begin:begin+blocklen-timeframe+1]:         #u2.trajectory[0,3]: means 0, 1, 2 --> last one doesn't count!
            u2.trajectory[begin]
            AL.alignto(segment, ref, select="name CA")
            startframe=H_atoms.atoms.coordinates()
            u2.trajectory[begin+timeframe-1]
            AL.alignto(segment, ref, select="name CA")
            stopframe=H_atoms.atoms.coordinates()
            msd_bl+=msd_win(lenH, startframe, stopframe)
            begin+=1
        msd_bl/=float(blocklen-timeframe+1)
    msd_bl/=len(selection)
    return msd_bl 


def main():
    pdb="unwrapped_disulfide_run.pdb"
    psf="unwrapped_disulfide_run.psf"
    dcd="powtraj.dcd"
    f=open('powder_250.dat', 'w')
    u1 = MDAnalysis.Universe(pdb)
    u2 = MDAnalysis.Universe(psf, dcd)
    timeframe=300
    begin=0
    blocklen=20000
    selection='ABCDEFGH'

    for i in range(0, len(u2.trajectory)/blocklen*blocklen, blocklen):
        begin=i
        print >> f, msd_block(selection, u1, u2, begin, timeframe, blocklen)


if __name__ == '__main__':
    main()

