#an attempt to write a code for the powder system
from scipy import *
from numpy import *
from MDAnalysis import *
from MDAnalysis.analysis.align import *
from string import ascii_uppercase
import numpy as np
import math
import MDAnalysis.KDTree.NeighborSearch as NS 

f=open('250_d_m_rmsd.dat', 'w')

#pdb="/home/sterpone_team/katava/results/disulfide/solution/2lym_dis_run.pdb"
#psf="/home/sterpone_team/katava/results/disulfide/solution/2lym_dis_run.psf"
#dcd1="/home/sterpone_team/katava/results/disulfide/solution/250.production.4.dcd"
#dcd2="/home/sterpone_team/katava/results/disulfide/solution/250.production.5.dcd"
pdb="2lym_wbi.pdb"
psf="2lym_wbi.psf"
dcd="traj.dcd"
u1 = MDAnalysis.Universe(pdb)
u2 = MDAnalysis.Universe(psf, dcd)
print(len(u2.trajectory))
inc=2
begin=0
end=begin+inc
fixedWindowWidth=float(end-begin)
window=len(u2.trajectory)/fixedWindowWidth
print('totalnrofwindows', window*fixedWindowWidth)

for i in 'U':
   segmentA=u2.selectAtoms("segid %s" %(i))
   print(segmentA)
   hy = segmentA.selectAtoms("name H*")
   car = segmentA.selectAtoms("name C*")
   ns_hy=NS.AtomNeighborSearch(hy)
   H_Catoms=ns_hy.search_list(car, 1.5)
   coordA=H_Catoms.atoms.coordinates()

   print(H_Catoms)
   finalval=[[0]*int(window) for x in range(len(H_Catoms))]
   finalval=np.array(finalval, dtype='float')

   timestepcount=-1
   begin=0
   end=inc
   ref = u1.selectAtoms("segid %s" %(i))
   trj = u2.selectAtoms("segid %s" %(i))
   refcoord=ref.atoms.coordinates()
   while(end<=fixedWindowWidth*window):
      x_avg=[0]*(len(H_Catoms))
      y_avg=[0]*(len(H_Catoms))
      z_avg=[0]*(len(H_Catoms))
      x_a=[0]*(len(H_Catoms))
      y_a=[0]*(len(H_Catoms))
      z_a=[0]*(len(H_Catoms))

      for ts in u2.trajectory[begin:end]:
         timestepcount+=1
         result=MDAnalysis.analysis.align.alignto(trj, ref, select="name CA") 
         print>>f, result
         coord=H_Catoms.atoms.coordinates()
        
         for j in range(0,len(coordA)):
            x_avg[j]+=coord[j][0]
            y_avg[j]+=coord[j][1]
            z_avg[j]+=coord[j][2]
            x_a[j]+=coord[j][0]*coord[j][0]
            y_a[j]+=coord[j][1]*coord[j][1]
            z_a[j]+=coord[j][2]*coord[j][2]
              
      for j in range(0,len(coordA)):
         x_avg[j]=(x_avg[j]/fixedWindowWidth)*(x_avg[j]/fixedWindowWidth)
         y_avg[j]=(y_avg[j]/fixedWindowWidth)*(y_avg[j]/fixedWindowWidth)
         z_avg[j]=(z_avg[j]/fixedWindowWidth)*(z_avg[j]/fixedWindowWidth)
         x_a[j]=x_a[j]/fixedWindowWidth
         y_a[j]=y_a[j]/fixedWindowWidth
         z_a[j]=z_a[j]/fixedWindowWidth
         finalval[j][begin/fixedWindowWidth]=1.0/3.0*(x_a[j]-x_avg[j]+y_a[j]-y_avg[j]+z_a[j]-z_avg[j])
      begin+=inc
      end+=inc
      np.savetxt('250_d_m_150.dat', finalval)
 
      if end>fixedWindowWidth*window:
         break
f.close()
