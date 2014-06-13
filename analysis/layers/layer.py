#!/usr/bin/env python

import MDAnalysis
import MDAnalysis.KDTree.NeighborSearch as NS
import MDAnalysis.analysis.align as AL
import numpy as np
import math

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
        tmpsum=0
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
            tmpsum+=msd_win(lenH, startframe, stopframe)
            begin+=1
        tmpsum/=float(blocklen-timeframe+1)    
        msd_bl+=tmpsum
    msd_bl/=float(len(selection))
    return msd_bl


def main():
    pdb="2lym_wb.pdb"
    ##psf="2lym_dis_run.psf"
    #dcd="250.production.dcd"
    #f=open('msd_monomer_250.dat', 'w')
    u1 = MDAnalysis.Universe(pdb)
    #u2 = MDAnalysis.Universe(psf, dcd)
    #timeframe=300
    #begin=0
    #blocklen=20000
    selection=['A']
    #print u1
    #for i in range(0, len(u2.trajectory)/blocklen*blocklen, blocklen):
    #    begin=i
    #    print >> f, msd_block(selection, u1, u2, begin, timeframe, blocklen)
    segment=u1.selectAtoms("segid U")
    hy=segment.selectAtoms("name H*")
    car=segment.selectAtoms("name C*")
    ns_hy=NS.AtomNeighborSearch(hy)
    H_atoms=ns_hy.search_list(car, 1.5)
    water=u1.selectAtoms("segid WT1" and "name OH2")
    
    Hidx=H_atoms.indices()
    Hpos=H_atoms.coordinates()
    
    Hdic={}
    count=0
    for i in Hidx:
        Hdic[i]=Hpos[count]
        count+=1

    Widx=water.indices()
    Wpos=water.coordinates()
    
    Wdic={}
    count2=0
    for j in Widx:
        Wdic[j]=Wpos[count2]
        count2+=1

    water_list={}
    for i in Hdic:
        list=[]
        for j in Wdic:
            list.append(math.sqrt((Hdic[i][0]-Wdic[j][0])**2+(Hdic[i][1]-Wdic[j][1])**2+(Hdic[i][2]-Wdic[j][2])**2))
        water_list[i]=list
    
    minimum_distance={}
    for i in water_list:
        minimum_distance[i]=min(water_list[i])

    valuemax=max(minimum_distance.values())
    valuemin=min(minimum_distance.values())

    d=(valuemax-valuemin)/3.0

    inner=[]
    middle=[]
    outer=[]
    for i in minimum_distance:
        if minimum_distance[i]<=d:
            outer.append(i)
        elif minimum_distance[i]<=2*d:
            middle.append(i)
        else:
            inner.append(i)

    ###### vmd selector_conversions:

    str_inner=''
    for i in inner:
        str_inner+=str(i)
        str_inner+=' '

    str_middle=''
    for i in middle:
        str_middle+=str(i)
        str_middle+=' '

    str_outer=''
    for i in outer:
        str_outer+=str(i)
        str_outer+=' '

    print 'str_inner', str_inner
    print 'str_middle', str_middle
    print 'str_outer', str_outer

    #water_distance={}
    #for i in Hdic:
    #    for j in Wdic:
    #        water_distance[i]=math.sqrt((Hdic[i][0]-Wdic[j][0])**2+(Hdic[i][1]-Wdic[j][1])**2+(Hdic[i][2]-Wdic[j][2])**2)


    #print Hdic.keys(), Hdic.values()
       
if __name__ == '__main__':
    main()


