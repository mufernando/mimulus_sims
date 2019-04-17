import allel
import subprocess, msprime, pyslim
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import re
from glob import glob
import pickle
import sys
import itertools
from timeit import default_timer as timer

def ac_from_ts(ts, n_pops):
    '''
    This function takes a tree sequence,  and returns tuple with a list of allele counts for each subpop and the positions'''
    acs=[]
    hap = allel.HaplotypeArray(ts.genotype_matrix())
    geno = hap.to_genotypes(ploidy=2)
    for i in range(n_pops):
        subpop_indexes = list(ts.samples(population=i))
        acs.append(geno.count_alleles(subpop=subpop_indexes))
    pos=np.array([s.position for s in ts.sites()])
    return(acs, pos)

def win_pi_sims(path, neut_mut, n_pops, n_sims, T, win_size, L):
    foname = os.path.basename(path[:-1])
    print(("Base filename:"+foname), flush=True)
    x = np.arange(n_pops)
    combs = list(itertools.combinations(x, 2))
    pis=np.zeros((len(T),n_sims,n_pops,int(L/win_size)))
    div=np.zeros((len(T), n_sims,len(combs),int(L/win_size)))
    fst=np.zeros((len(T), n_sims,len(combs),int(L/win_size)))
    tajd=np.zeros((len(T),n_sims,n_pops,int(L/win_size)))
    for t in range(len(T)):
        for i in range(n_sims):
            files = glob(path+str(T[t])+"N_sim_"+str(i)+"_RAND_*[0-9]_overlaid.trees")
            print(files)
            assert (len(files) == 1), str(len(files))+" file(s) found with glob T: "+str(T[t])+" sim:"+str(i)
            filename= files[0]
            print(filename)
            ts = pyslim.load(filename).simplify()
            #print(("Pi0: ", ts.pairwise_diversity(samples=ts.samples(population=0)),"Pi1: ", ts.pairwise_diversity(samples=ts.samples(population=1))), flush=True)
            s1 = timer()
            acs, pos = ac_from_ts(ts, n_pops)
            for j in range(n_pops):
                pi, windows, n_bases, counts = allel.windowed_diversity(pos, acs[j], size=win_size, start=1, stop=L)
                pis[t,i,j,:] = pi
                D, windows, counts = allel.windowed_tajima_d(pocs, acs[j], size=win_size, start=1, stop=L)
                tajd[t,i,j,:] = D
            s2 = timer()
            print(("Calculating windowed Pi/TajD... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
            s1=timer()
            for k in range(len(combs)):
                dxy, windows, n_bases, counts = allel.windowed_divergence(pos, acs[combs[k][0]], acs[combs[k][0]], size=win_size, start=1, stop=L)
                div[t,i,k,:] = dxy
                fstat, windows, counts = allel.windowed_hudson_fst(pos, acs[combs[k][0]], acs[combs[k][0]], size=win_size, start=1, stop=L)
                fst[t,i,k,:] = fstat
            s2 = timer()
            print(("Calculating windowed Dxy and Fst... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

    s1 = timer()
    print((pis.shape), flush=True)
    print((tajd.shape), flush=True)
    print((div.shape), flush=True)
    output = open(path+foname+'_pis.pkl', 'wb')
    pickle.dump(pis, output)
    output.close()
    output = open(path+foname+'_tajd.pkl', 'wb')
    pickle.dump(tajd, output)
    output.close()
    output = open(path+foname+'_div.pkl', 'wb')
    pickle.dump(div, output)
    output.close()
    output = open(path+foname+'_fst.pkl', 'wb')
    pickle.dump(fst, output)
    output.close()

    plt.subplot(2, 1, 1)
    plt.plot(np.transpose(pis[0,0,:]), "-")
    plt.title('0N after split')
    plt.ylabel('Pi')
    plt.subplot(2, 1, 2)
    plt.plot(np.transpose(pis[9,0,:]), "-")
    plt.title('10N after split')
    plt.xlabel('Window')
    plt.ylabel('Pi')
    plt.tight_layout()
    plt.savefig(path+foname+'_landscape.pdf')
    plt.close()

    s2 = timer()
    print(("Saving stats and plots to file... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

#sys.stdout.flush()

s1 = timer()
path = sys.argv[1]
win_size=int(sys.argv[2])
L=int(sys.argv[3])
n_sims=int(sys.argv[4])

neut_mut = 1e-8
n_pops = 2

#T=np.arange(0.1,1.1,step=0.1)
#T=np.concatenate([T,np.array([2.0])])
#T=np.around(T, 1)
#T=np.arange(0,11,step=1)
T=np.concatenate([np.arange(0,2.2,step=0.4), np.arange(4,11,step=2)])
T = [float("%.1f"%t) for t in T]

s2 = timer()
print(("Initializing... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

win_pi_sims(path, neut_mut, n_pops, n_sims, T, win_size, L)

"""

path = "/Users/murillo/Drive/phd/w19/rotation/trees/"
folders = glob("N_*_mu_*_L0_*")
paths = [path+folders[i]+"/" for i in range(len(folders))]

total_mut = 1e-8
n_pops = 2
n_sims=1
win_size = 500000
T=np.arange(1,11,step=1)

for i in range(len(paths)):
    path = paths[i]
    foname = os.path.basename(path[:-1])
    print((foname), flush=True)
    matches = re.match( r'N_(\d+)_mu_(.*)_r_(.*)_deff_(\d+)_L_(\d+)_L0_(.+)_L1_(\d+)', foname)
    N, mu, r, deff, L, L0, L1 = matches.groups()
    del_mut = mu
    win_pi_sims(path, neut_mut, del_mut, n_pops, n_sims, T, win_size)"""

#python3 mimulus_win_pi.py /Users/murillo/Drive/phd/w19/rotation/trees/N_10000_mu_0_r_2e-08_deff_0_L_20000000_L0_0_L1_0_m_0/ 500000 1
