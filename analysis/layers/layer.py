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
    msd_bl=[0]*3
    tmpbegin=begin
    for i in selection:
        listi={}
        print i
        dicsel=select_layers(u1, u2, i="A")
        for item in dicsel:
            tmpsum=0
            begin=tmpbegin
            segment=u2.selectAtoms("segid %s" %(i))
            ref=u1.selectAtoms("segid %s" %(i))
            hy=segment.selectAtoms("name H*")
            car=segment.selectAtoms("name C*")
            ns_hy=NS.AtomNeighborSearch(hy)
            allH_atoms=ns_hy.search_list(car, 1.5)
            H_atoms=allH_atoms.selectAtoms(*dicsel[item])
            lenH=len(H_atoms)
            print lenH
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
            listi[item]=tmpsum
        msd_bl[0]+=listi['inner']
        msd_bl[1]+=listi['middle']
        msd_bl[2]+=listi['outer']
    msd_bl[0]/=float(len(selection))
    msd_bl[1]/=float(len(selection))
    msd_bl[2]/=float(len(selection))
    open("250_inner.dat", 'a').write("%s\n" %(msd_bl[0]))
    open("250_middle.dat", 'a').write("%s\n" %(msd_bl[1]))
    open("250_outer.dat", 'a').write("%s\n" %(msd_bl[2]))


def select_layers(u1, u2, i="A"):
    
    segment=u1.selectAtoms("segid %s" %(i))
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
    
    inner_selection=[]
    middle_selection=[]
    outer_selection=[]
    for i in range(len(H_atoms.atoms)):
        if H_atoms.atoms[i].number in inner:
            inner_selection.append("segid %s and resid %d and name %s" %(H_atoms.atoms[i].segid,  H_atoms.atoms[i].resid,  H_atoms.atoms[i].name))
        if H_atoms.atoms[i].number in middle:
            middle_selection.append("segid %s and resid %d and name %s" %(H_atoms.atoms[i].segid,  H_atoms.atoms[i].resid,  H_atoms.atoms[i].name))
        if H_atoms.atoms[i].number in outer:
            outer_selection.append("segid %s and resid %d and name %s" %(H_atoms.atoms[i].segid,  H_atoms.atoms[i].resid,  H_atoms.atoms[i].name))
        


    print len(inner_selection)
    inner_selection=tuple(inner_selection)
    
    print len(middle_selection)
    middle_selection=tuple(middle_selection)
    
    print len(outer_selection)
    outer_selection=tuple(outer_selection)

    dicsel={}
    dicsel['inner']=inner_selection
    dicsel['middle']=middle_selection
    dicsel['outer']=outer_selection
    return dicsel


def main():
    pdb="2lym_dis_run.pdb"
    psf="2lym_dis_run.psf"
    dcd="shortest.dcd"
    u1 = MDAnalysis.Universe(pdb)
    u2 = MDAnalysis.Universe(psf, dcd)
    #timeframe=300
    #begin=0
    #blocklen=20000
    timeframe=3
    begin=0
    blocklen=5
    selection=['A']
    ###### vmd selector_conversions:

    #str_inner=''
    #for i in inner:
    #    str_inner+=str(i)
    #    str_inner+=' '

    #str_middle=''
    #for i in middle:
    #    str_middle+=str(i)
    #    str_middle+=' '

    #str_outer=''
    #for i in outer:
    #    str_outer+=str(i)
    #    str_outer+=' '

    #print 'str_inner', str_inner
    #print 'str_middle', str_middle
    #print 'str_outer', str_outer

    for i in range(0, len(u2.trajectory)/blocklen*blocklen, blocklen):
        begin=i
        msd_block(selection, u1, u2, begin, timeframe, blocklen)

if __name__ == '__main__':
    main()


