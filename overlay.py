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

def remove_mutations(ts, start, end, proportion):
    '''
    This function will return a new tree sequence the same as the input,
    but after removing each non-SLiM mutation within regions specified in lists
    start and end with probability `proportion`, independently. So then, if we
    want to add neutral mutations with rate 1.0e-8 within the regions and 0.7e-8
    outside the regions, we could do
      ts = pyslim.load("my.trees")
      first_mut_ts = msprime.mutate(ts, rate=1e-8)
      mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)
    :param float proportion: The proportion of mutations to remove.
    '''
    new_tables = ts.dump_tables()
    new_tables.mutations.clear()
    mutation_map = [-1 for _ in range(ts.num_mutations)]
    for j, mut in enumerate(ts.mutations()):
        keep_mutation = True
        for i in range(len(start)):
            left = start[i]
            right = end[i]
            assert(left < right)
            if i > 0:
                assert(end[i - 1] <= left)
            if mut.position >= left and mut.position < right and len(mut.metadata) == 0:
                keep_mutation = (random.uniform(0, 1) > proportion)
        if keep_mutation:
            mutation_map[j] = new_tables.mutations.num_rows
            if mut.parent < 0:
                new_parent = -1
            else:
                new_parent = mutation_map[mut.parent]
            new_tables.mutations.add_row(site = mut.site, node = mut.node,
                    derived_state = mut.derived_state,
                    parent = new_parent,
                    metadata = mut.metadata)
    return new_tables.tree_sequence()

def extract_meta(fpath, gene_info, type):
    start = []
    end = []
    flist = fpath.split('/')
    foname, fname = flist[-2], flist[-1]
    if (type == "gene"):
        print(foname, file=sys.stderr)
        matches = re.match( r'N_(\d+)_mu_(.*)_r_(.*)_deff_(\d+)_L_(\d+)', foname)
        N, del_mut, r, deff, L = matches.groups()
        del_mut = float(del_mut)
        L = int(L)
        if del_mut > 0 :
            start_end = np.loadtxt(gene_info, dtype="int")
            start = start_end[:,0]
            end = start_end[:,1]
    else:
        print(foname, file=sys.stderr)
        #matches = re.match( r'N_(\d+)_mu_(.*)_r_(.*)_deff_(.+)_L_(\d+)_L0_(.+)_L1_(\d+)', foname)
        #N, del_mut, r, deff, L, L0, L1 = matches.groups()
        #matches = re.match( r'N_(\d+)_mu_(.*)_r_(.*)_L_(\d+)_L0_(.+)_L1_(\d+)_m_(.+)_ndeff_(.+)_nprop_(.+)_pdeff_(.+)_pprop_(.+)', foname)
        #N, del_mut, r, L, L0, L1, m,ndeff,nprop,pdeff,pprop  = matches.groups()
        matches = re.match( r'N_(\d+)_mu_(.*)_r_(.*)_L_(\d+)_L0_(.+)_L1_(\d+)_m_*', foname)
        N, del_mut, r, L, L0, L1  = matches.groups()
        L, L0, L1 = int(L), int(L0),int(L1)
        del_mut = float(del_mut)
        if del_mut > 0:
            start = []
            end = []
            for i in range(0, L, (L0+L1)):
                start.append(i)
                end.append((i+L1)-1)
    return(start, end, del_mut)

def overlay_varmut(fpath, neut_mut, gene_info, type="gene"):
    flist = fpath.split('/')
    foname, fname = flist[-2], flist[-1]
    ts_slim = pyslim.load(fpath).simplify()
    ts_mut = msprime.mutate(ts_slim, neut_mut, keep=True)
    print("Mutated", fpath, "in msprime...", flush=True)
    start, end, del_mut = extract_meta(fpath, gene_info, type)
    if del_mut > 0:
        s1=timer()
        ts = remove_mutations(ts_mut, start, end, del_mut/neut_mut)
        s2 = timer()
        print(("Removed extra mutations in genic regions from", fpath, "... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
    else:
        ts = ts_mut
    s1=timer()
    ts.dump(fpath[:-6]+"_overlaid.trees")
    s2 = timer()
    print(("Dumped overlaid trees to file", fpath, "... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

fpath = sys.argv[1]
type = sys.argv[2]
gene_info = sys.argv[3]

neut_mut=1e-8

overlay_varmut(fpath, neut_mut, gene_info, type)
