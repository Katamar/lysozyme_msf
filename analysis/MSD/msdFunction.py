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
    for i in selection:
        segment=u2.selectAtoms("segid %s" %(i))
        ref=u1.selectAtoms("segid %s" %(i))
        hy=segment.selectAtoms("name H*")
        car=segment.selectAtoms("name C*")
        ns_hy=NS.AtomNeighborSearch(hy)
        H_atoms=ns_hy.search_list(car, 1.5)
        lenH=len(H_atoms)
        #ctr=0  #delete later
        for ts in u2.trajectory[begin:begin+blocklen-1]:
            #ctr+=1
            u2.trajectory[begin]
            AL.alignto(segment, ref, select="name CA")
            startframe=H_atoms.coordinates()
            u2.trajectory[begin+timeframe-1]
            AL.alignto(segment, ref, select="name CA")
            stopframe=H_atoms.atoms.coordinates()
            #print 'begin, end', begin, begin+timeframe-1
            msd_bl+=msd_win(lenH, startframe, stopframe)
            #print msd_bl
            begin+=1
        #print 'i went through this many timewindows, theoretical value', ctr, blocklen-timeframe+1 
        msd_bl/=float(blocklen-timeframe+1)
    msd_bl/=len(selection)
    return msd_bl 


def main():
    #pdb="2lym_wbi.pdb"
    #psf="2lym_wbi.psf"
    #dcd="traj.dcd"
    pdb="2lym_dis_run.pdb"
    psf="2lym_dis_run.psf"
    dcd1="250.production.4.dcd"
    dcd2="250.production.5.dcd"
    #dcd2="u.250.production.5.dcd"
    f=open('msd_250_monomer.dat', 'w')
    u1 = MDAnalysis.Universe(pdb)
    u2 = MDAnalysis.Universe(psf, [dcd1, dcd2])
    print len(u2.trajectory)
    timeframe=300
    begin=0
    blocklen=10000
    selection='A'

    #print msd_block(selection, u1, u2, begin, timeframe, blocklen)
    for i in range(0, (len(u2.trajectory)/blocklen)*blocklen, blocklen):
        print 'new block starts here', i
        begin=i
        print >> f, msd_block(selection, u1, u2, begin, timeframe, blocklen)


if __name__ == '__main__':
    main()

