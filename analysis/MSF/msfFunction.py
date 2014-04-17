import MDAnalysis
import MDAnalysis.KDTree.NeighborSearch as NS
import MDAnalysis.analysis.align as AL
import numpy as np

def accumulate_positions(lenc, coord, x_avg, y_avg, z_avg, x_a, y_a, z_a):
    for j in range(0,lenc):
        x_avg[j]+=coord[j][0]
        y_avg[j]+=coord[j][1]
        z_avg[j]+=coord[j][2]
        x_a[j]+=coord[j][0]*coord[j][0]
        y_a[j]+=coord[j][1]*coord[j][1]
        z_a[j]+=coord[j][2]*coord[j][2]
    return x_avg, y_avg, z_avg, x_a, y_a, z_a

def average_positions(lenc, fixedWindowWidth, x_avg, y_avg, z_avg, x_a, y_a, z_a):
    finalval=0
    for j in range(0,lenc):
        x_avg[j]=(x_avg[j]/fixedWindowWidth)*(x_avg[j]/fixedWindowWidth)
        y_avg[j]=(y_avg[j]/fixedWindowWidth)*(y_avg[j]/fixedWindowWidth)
        z_avg[j]=(z_avg[j]/fixedWindowWidth)*(z_avg[j]/fixedWindowWidth)
        x_a[j]=x_a[j]/fixedWindowWidth
        y_a[j]=y_a[j]/fixedWindowWidth
        z_a[j]=z_a[j]/fixedWindowWidth
        finalval+=1.0/3.0*(x_a[j]-x_avg[j]+y_a[j]-y_avg[j]+z_a[j]-z_avg[j])
    return finalval

def msf_win(begin, end, fixedWindowWidth, u1, u2, selection):
    segment=u2.selectAtoms("segid %s" %(selection[0:1]))
    hy = segment.selectAtoms("name H*")
    car = segment.selectAtoms("name C*")
    ns_hy=NS.AtomNeighborSearch(hy)
    H_Catoms=ns_hy.search_list(car, 1.5)
    lenH=len(H_Catoms)
    protein=0                     
    for i in selection:
        x_avg=[0]*lenH
        y_avg=[0]*lenH
        z_avg=[0]*lenH
        x_a=[0]*lenH
        y_a=[0]*lenH
        z_a=[0]*lenH
        finalval=0

        segment=u2.selectAtoms("segid %s" %(i))
        hy = segment.selectAtoms("name H*")
        car = segment.selectAtoms("name C*")
        ns_hy=NS.AtomNeighborSearch(hy)
        H_Catoms=ns_hy.search_list(car, 1.5)
        lenH=len(H_Catoms)
        ref = u1.selectAtoms("segid %s" %(i))
        trj = u2.selectAtoms("segid %s" %(i))
                                                                                                

        for ts in u2.trajectory[begin:end]:
            result=AL.alignto(trj, ref, select="name CA")
            coord=H_Catoms.atoms.coordinates()
            accumulate_positions(len(coord), coord, x_avg, y_avg, z_avg, x_a, y_a, z_a)
        win_avg=average_positions(len(coord), fixedWindowWidth, x_avg, y_avg, z_avg, x_a, y_a, z_a)
        win_avg/=float(len(coord))
        protein+=win_avg
    return protein/float(len(selection)) 

def blocksum(blocklen, fixedWindowWidth, begin, u1, u2, selection):
    blocks=0
    count=0
    for i in range(begin, begin+(blocklen/fixedWindowWidth)*fixedWindowWidth, fixedWindowWidth):
        count+=1
        begin=i
        end=i+fixedWindowWidth
        blocks+=msf_win(begin, end, fixedWindowWidth, u1, u2, selection)
    blocks/=float(count)
    return blocks




def main():
    pdb="unwrapped_disulfide_run.pdb"
    psf="unwrapped_disulfide_run.psf"
    dcd1="u.250.production.4.dcd"
    dcd2="u.250.production.5.dcd"
    dcd="powtraj.dcd"
    f=open('pow_250.dat', 'w') 

    u1 = MDAnalysis.Universe(pdb)
    u2 = MDAnalysis.Universe(psf, [dcd1, dcd2])
    
    fixedWindowWidth=300
    begin=0
    blocklen=20000
    selection='ABCDEFGH'

    for i in range(0, (len(u2.trajectory)/blocklen)*blocklen, blocklen):
        begin=i
        print >> f, blocksum(blocklen, fixedWindowWidth, begin, u1, u2, selection)

if __name__ == "__main__":
    main()
