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
        ctr=0  #delete later
        #print begin, begin+blocklen-2
        for ts in u2.trajectory[begin:begin+blocklen-timeframe+1]:         #u2.trajectory[0,3]: means 0, 1, 2 --> last one doesn't count!
            ctr+=1
            u2.trajectory[begin]
            AL.alignto(segment, ref, select="name CA")
            startframe=H_atoms.atoms.coordinates()
            u2.trajectory[begin+timeframe-1]
            AL.alignto(segment, ref, select="name CA")
            stopframe=H_atoms.atoms.coordinates()
            #if begin>35000:
            print 'begin, end', begin, begin+timeframe-1
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
    #pdb="2lym_dis_run.pdb"
    #psf="2lym_dis_run.psf"
    #dcd1="250.production.4.dcd"
    #dcd2="250.production.5.dcd"
    #dcd="250.production.5.dcd"
    pdb="unwrapped_disulfide_run.pdb"
    psf="unwrapped_disulfide_run.psf"
    dcd="powtraj.dcd"
    f=open('powder_250.dat', 'w')
    u1 = MDAnalysis.Universe(pdb)
    u2 = MDAnalysis.Universe(psf, dcd)
    #u2 = MDAnalysis.Universe(psf, dcd)
    print len(u2.trajectory)
    timeframe=3
    begin=0
    blocklen=7
    selection='ABCDEFGH'

    #print msd_block(selection, u1, u2, begin, timeframe, blocklen)
    for i in range(0, len(u2.trajectory)/blocklen*blocklen, blocklen):
        print 'new block starts here', i
        print 'ovo je zadnja petlja ikada'
        begin=i
        print >> f, msd_block(selection, u1, u2, begin, timeframe, blocklen)


if __name__ == '__main__':
    main()

