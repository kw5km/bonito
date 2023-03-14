# env: Nnaopore
#!/usr/bin/python
import sys
import matplotlib.pyplot as plt
# sys.path.append('/home/kw5km/work/mnt/kw5km/Research/Nanopore_Channel/Code/Viterbi/')
# import viterbi
# sys.path.append('/home/kw5km/work/mnt/kw5km/Research/Nanopore_Channel/Code/utils/')
# import squiggle_helper
import numpy as np
import os


name = 'chiron_ecoli' #'bonito'
# dat_ref_ = np.load('/home/work/mnt/entropy/code/Nanopore/bonito/bonito/data/dna_r9.4.1/references.npy', mmap_mode='r')

dat_ref_ = []
lengths = []
for f in range(1,1000):
    if f%100==0:print('f:', f)
    bs = np.genfromtxt('/home/kw5km/work/mnt/entropy/code/Nanopore/bonito/bonito/data/train/ecoli_{0:04d}.label'.format(f), usecols=[2], dtype='str')
    lengths.append(bs.size)
    pad = bs.size%480
    bs = np.append(bs,['N']*(480-pad))

    for ti, t in enumerate(range(0, bs.size, 480)):
        dat_ref_.append(bs[t:t+480])

dat_ref_ = np.array(dat_ref_)
print('length sum', np.sum(lengths))
print('length means, std', np.mean(lengths), np.std(lengths))
for b1, b2 in zip(range(0, 900, 100), range(100, 1000, 100)):
    k = 6

    dat_bases = dat_ref_[b1:b2, :]
    # dat_ref = dat_ref_[b1:b2, :]
    # labels = ["N", "A", "C", "G", "T"]
    # dat_bases = np.array(labels)[np.array(dat_ref)]

    kmers_dat = np.empty_like(dat_bases, dtype=object)
    for i in range(dat_bases.shape[0]):
        ks = np.array(squiggle_helper.rolling_window(dat_bases[i,:], k))
        kmers_dat[i, :] = ks

    # print(kmers_dat)

    kmers = viterbi.get_all_kmers(k, ['A', 'C', 'G', 'T'])
    edges_, edges_set_ = viterbi.get_debruijn_edges_from_kmers(kmers)
    scores, scores_vals = viterbi.get_scores_m(edges_, kmers, score_file='/home/work/mnt/kw5km/Research/Nanopore_Channel/Regression_Output/6mer_means.txt')

    for delta in range(1, 10):
        

        break_idx = np.zeros((kmers_dat.shape[0], kmers_dat.shape[1]-1))

        for j in range(kmers_dat.shape[0]):
            for i, (kmer1, kmer2) in enumerate(zip(kmers_dat[j, :-1], kmers_dat[j, 1:])):
                # if (kmer1=='N') or (kmer2=='N'): break
                if ('N' in kmer1) or ('N' in kmer2):
                    break_idx[j, i] = -1
                elif (np.abs(scores_vals[kmer1]-scores_vals[kmer2])<delta) and (kmer1!=kmer2):
                    break_idx[j, i] = 1

        # Calculate (height_of_image / width_of_image)
        f2, ax2 = plt.subplots()
        im = ax2.imshow(break_idx, vmin = np.min(break_idx), vmax = np.max(break_idx))  
        im_ratio = break_idx.shape[0]/break_idx.shape[1]

        cbar = f2.colorbar(im, fraction=0.075*im_ratio, ticks=[-1, 0, 1])
        cbar.ax.set_yticklabels(['Contains N', 'No Constraint Break', 'Constraint Break'])
        # cbar.set_label('')

        plt.title(r'{} Training Data $\delta$ = {}'.format(name, delta))
        plt.tight_layout()
        plt.savefig('/home/work/mnt/entropy/data/Nanopore/bonito/data_stats/{}_{}_{}_delta{}.eps'.format(name, b1, b2, delta))
        plt.close()