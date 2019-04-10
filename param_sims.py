#!/usr/bin/env python3
#
import numpy as np
import os
import sys

sys.stdout = open('out', 'w')

#change to gamma shape=2
nprop = np.array([0.05, 0.10])
ndeff = np.array([10,100])
pprop = np.array([0.001, 0.005])
pdeff = np.array([100])
a1 =np.stack(np.meshgrid(nprop,ndeff, pprop, pdeff),-1).reshape(-1,4)
a2 = np.stack(np.meshgrid(np.array([0]),np.array([0]), pprop, pdeff),-1).reshape(-1,4)
a3 = np.stack(np.meshgrid(nprop,ndeff, np.array([0]), np.array([0])),-1).reshape(-1,4)
a4 = np.array([[0.,0.,0.,0.]])

params = np.concatenate([a4,a3,a2,a1])
print(params)
N=str(10000)
total_mut = 1e-8
rec=str(2e-8)
L = str(21000000)
L0 = str(7000000)
L1 = str(7000000)
m=str(0)
nsim=10

for row in params:
    for sim in range(nsim):
        sim = str(sim)
        #print(row[0],row[2])
        mut = (float(row[0])+float(row[2]))*total_mut
        if float(row[0])+float(row[2]) >  0:
            mut = '{:.2e}'.format(mut)
        row = row.astype(np.str)
        mut=str(mut)
        ans = os.system(" ".join(['sbatch --open-mode=append run_dual.srun ',sim,N,mut,rec,L,L0,L1,m," ".join(row)]))
        print(" ".join(['sbatch --open-mode=append run_dual.srun ',sim,N,mut,rec,L,L0,L1,m," ".join(row)]), "Response: ", ans)

### params local adapt
params2 = np.stack(np.meshgrid(pprop,pdeff),-1).reshape(-1,2)
recipe_name="mimulus_local.slim"
for row in params2:
    for sim in range(nsim):
        sim = str(sim)
        #print(row[0],row[2])
        mut = float(row[0])*total_mut*3#bc there are three mut types and at one point only one of them is not neutral for that particular subpop
        mut = '{:.2e}'.format(mut)
        row = row.astype(np.str)
        mut=str(mut) #
        ans1 = os.system(" ".join(['sbatch --open-mode=append run_one.srun ',sim,N,mut,rec,L,L0,L1,'0.1',str(row[1]),recipe_name]))
        ans2 = os.system(" ".join(['sbatch --open-mode=append run_one.srun ',sim,N,mut,rec,L,L0,L1,'1',str(row[1]),recipe_name]))
        print(" ".join(['sbatch --open-mode=append run_one.srun ',sim,N,mut,rec,L,L0,L1,'0.1',str(row[1])]), "Response: ", ans1)
        print(" ".join(['sbatch --open-mode=append run_one.srun ',sim,N,mut,rec,L,L0,L1,'1',str(row[1])]), "Response: ", ans2)

### params BDM introgression
params3 = np.stack(np.meshgrid(nprop,ndeff),-1).reshape(-1,2)
recipe_name="mimulus_intro.slim"
for row in params3:
    for sim in range(nsim):
        sim = str(sim)
        #print(row[0],row[2])
        mut = float(row[0])*total_mut/0.9#bc there are three mut types and at one point only one of them is not neutral for that particular subpop
        mut = '{:.2e}'.format(mut)
        row = row.astype(np.str)
        mut=str(mut) #
        ans1 = os.system(" ".join(['sbatch --open-mode=append run_one.srun ',sim,N,mut,rec,L,L0,L1,'0.1',str(row[1]), recipe_name]))
        ans2 = os.system(" ".join(['sbatch --open-mode=append run_one.srun ',sim,N,mut,rec,L,L0,L1,'1',str(row[1]), recipe_name]))
        #ans1 = ans2 = "la"
        print(" ".join(['sbatch --open-mode=append run_one.srun ',sim,N,mut,rec,L,L0,L1,'0.1',str(row[1])]), "Response: ", ans1)
        print(" ".join(['sbatch --open-mode=append run_one.srun ',sim,N,mut,rec,L,L0,L1,'1',str(row[1])]), "Response: ", ans2)
