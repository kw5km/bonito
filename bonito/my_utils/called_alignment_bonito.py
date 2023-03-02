# env: Nanopore
#!/usr/bin/python
import sys
import numpy as np
import os
from pathlib2 import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import SeqIO
from Bio.Seq import Seq
import edlib



def format_cigar(cigar):

    ins_idx = [int(idx)-1 for idx, val in enumerate(cigar) if val=='I']
    dels_idx = [int(idx)-1 for idx, val in enumerate(cigar) if val=='D']
    subs_idx = [int(idx)-1 for idx, val in enumerate(cigar) if val=='X']

    # print('ins_idx:', ins_idx, 'dels_idx:', dels_idx, 'subs_idx:', subs_idx)

    cigar_array = np.array(list(cigar))
    ins_counts, dels_counts, subs_counts = cigar_array[ins_idx].astype(float), cigar_array[dels_idx].astype(float), cigar_array[subs_idx].astype(float)
    # print('ins_counts:', ins_counts, 'dels_counts:', dels_counts, 'subs_counts:', subs_counts)
    return np.sum(ins_counts), np.sum(dels_counts), np.sum(subs_counts)

plt.rcParams.update({'font.size': 14})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

n, bp, ver = 1000, 200, '3-4'

bon_mean = []
bon_std = []
bon_indelsub = []

rates = ["106", "106", "106", "106", "106", "106", "86", "86", "86", "86", "66"]
deltas = ['000', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1000']

for delt, r in zip(deltas, rates):
    
    print('DELTA:', delt)
    cwd = os.getcwd()

    # 0 is entirely unconstrained so no state split
    if delt == '000':
        path_var_ref = '/home/kw5km/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/deBruijn_6mer_mean/n{}_bp{}_delta{}_DeepSimu/fasta/'.format(n, bp, delt)
        basecall_file = 'n{}_bp{}_delta{}@v{}.fastq'.format(n, bp, delt, ver)
        f_out = 'n{}_bp{}_delta{}@v{}/'.format(n, bp, delt, ver)

    else:
        path_var_ref = '/home/kw5km/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/deBruijn_6mer_mean/Encoded/n{}_bp{}_delta{}_statesplit_rate{}_DeepSimu/fasta/'.format(n, bp, delt, r)
        basecall_file = 'n{}_bp{}_delta{}_rate{}@v{}.fastq'.format(n, bp, delt, r, ver)
        f_out = 'n{}_bp{}_delta{}_rate{}@v{}/'.format(n, bp, delt, r, ver)

    score_f, fig_f = '/scores/', '/Figs/'
    path_var_long = '/home/kw5km/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/bonito/basecalls/'
    path_var_long_out = path_var_long+ f_out 

    bon_calls = {}
    filename_bon = path_var_long+basecall_file  
    print('Bonito file:', filename_bon)
    with open(filename_bon) as bon_file:
        for record in SeqIO.parse(bon_file, "fastq"):
            r_desc = record.description
            id_start, id_end = r_desc.find('signal'), r_desc.find(".fast5")
            id = r_desc[id_start:id_end]
            bon_calls[id] = record.seq
        bon_file.close() 
        
    scores_bon = []
    lev_bon = []
    errs_bon = []

    for signal_id in bon_calls.keys():

        print(signal_id)
        ''' Bonito '''
        seq_bon = bon_calls[signal_id]

        ''' True '''
        filename_true = path_var_ref+signal_id+'.fasta'
        with open(filename_true) as true_file:
            for line in true_file:
                if '>' not in line: bases_true = line
            true_file.close()
        seq_true = Seq(bases_true)

        
        ''' Align '''

        bon_align = edlib.align(str(seq_bon), str(seq_true), task = "path")
        bon_alignment = edlib.getNiceAlignment(bon_align, str(seq_bon), str(seq_true))
        lev_bon.append(bon_align["editDistance"])
        bon_ins, bon_dels, bon_subs = format_cigar(bon_align['cigar'])
        errs_bon.append((bon_ins, bon_dels, bon_subs))


        filepath_align = path_var_long_out+score_f
        if not os.path.exists(filepath_align):
            os.makedirs(filepath_align)

        filename_align = filepath_align+'align_{}.txt'.format(signal_id)
        with open(filename_align, "w") as txt_file:
            txt_file.write('(top seq<-true) \n')
            txt_file.write('Bonito Version {}: \n'.format(ver))
            txt_file.write("\n".join(bon_alignment.values())) 
            txt_file.write('\n')
            txt_file.write('Lev Dist: {}'.format(bon_align["editDistance"]))
            txt_file.write('\n')
            txt_file.write('Edits: {}'.format(bon_align["cigar"]))
            txt_file.write('\n')
            txt_file.write('Insertions: {}, Deletions: {}, Substitutions: {}'.format(bon_ins, bon_dels, bon_subs))
            txt_file.write('\n')
            txt_file.close()

#         # sig = np.loadtxt(path_var_long+'signal/signal_{}.txt'.format(j))
#         # sig = (sig-14)/5.7
#         # time = np.arange(sig.size)
#         # for ti in range(0, time.size-200, 200):
            
#         #     t1, t2 = ti, ti+200
            
#         #     fig, axs = plt.subplots(3, 1)
#         #     fig.tight_layout()
            
#         #     axs[0].set_ylim(50,125)
#         #     axs[1].set_ylim(50,125)
#         #     axs[2].set_ylim(50,125)
            
#         #     axs[0].plot(time[t1:t2],sig[t1:t2], label='DeepSimu Current')
#         #     axs[1].plot(time[t1:t2],sig[t1:t2], label='DeepSimu Current')
#         #     axs[2].plot(time[t1:t2],sig[t1:t2], label='DeepSimu Current')
            
#         #     axs[0].title.set_text('True Segments')
#         #     axs[1].title.set_text('Viterbi Segments')
#         #     axs[2].title.set_text('ONT Segments')
            
            
#         #     for l in range(time.size):
#         #         if(time[l]>t2):break
            
#         #         if l < df_b.shape[0]:
#         #             if (df_b[l]<=t2) and (df_b[l]>=t1):
#         #                 axs[0].axvline(df_b[l], color='green', alpha=0.15, label='True Boundary')
#         #                 axs[1].axvline(df_b[l], color='green', alpha=0.15, label='True Boundary')
#         #                 axs[2].axvline(df_b[l], color='green', alpha=0.15, label='True Boundary')
                        
                        
#         #         ## Viterbi Segs
#         #         if l < len(fn_viterbi):
#         #             if (fn_viterbi[l]<=t2) and (fn_viterbi[l]>=t1):
#         #                 axs[1].axvline(fn_viterbi[l], color='black', label='False Negaive')
#         #         if l < db_viterbi.size:
#         #             if (db_viterbi[l]<=t2) and (db_viterbi[l]>=t1):                       
#         #                 axs[1].axvline(db_viterbi[l], color='blue', alpha=0.5, label='False Positive')
                
#         #         # ONT Segs
#         #         # if l < len(fn_ont):
#         #         #     if (fn_ont[l]<=t2) and (fn_ont[l]>=t1):
#         #         #         axs[2].axvline(fn_ont[l], color='black', label='False Negaive')
#         #         # if l < db_ont.size:
#         #         #     if (db_ont[l]<=t2) and (db_ont[l]>=t1):                       
#         #         #         axs[2].axvline(db_ont[l], color='blue', alpha=0.5, label='False Positive')
                                
#         #     handles, labels = plt.gca().get_legend_handles_labels()
#         #     by_label = OrderedDict(zip(labels, handles))
#         #     plt.legend(by_label.values(), by_label.keys(), loc='upper right')
            
#         #     if not os.path.exists(filename_figs_base+'/'+str(t1)+'_'+str(t2)):
#         #         os.makedirs(filename_figs_base+'/'+str(t1)+'_'+str(t2))
            
#         #     filename_figs = filename_figs_base + '/'+str(t1)+'_'+str(t2) +'/'+'j{}'.format(j)
#         #     plt.savefig(filename_figs+'_'+str(t1)+'_'+str(t2)+'.png', dpi=100)
#         #     plt.close()

#     avg_err_vit = np.mean(errs_vit, axis=0)

#     avg_err_guppy = np.mean(errs_guppy, axis=0)
#     avg_err_raw = np.mean(errs_raw, axis=0)
#     avg_err_ev = np.mean(errs_ev, axis=0)

#     filename_scores = filepath_align+'/alignment_scores.txt'
#     with open(filename_scores, "w") as txt_file:
#         txt_file.write('Viterbi: \n')
#         txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_vit), np.std(lev_vit)))
#         txt_file.write('\n')
#         txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_vit[0], avg_err_vit[1], avg_err_vit[2]))
#         txt_file.write('\n')
#         txt_file.write('\n')
#         txt_file.write('Guppy: \n')
#         txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_guppy), np.std(lev_guppy)))
#         txt_file.write('\n')
#         txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_guppy[0], avg_err_guppy[1], avg_err_guppy[2]))
#         txt_file.write('\n')
#         txt_file.write('\n')
#         txt_file.write('Raw: \n')
#         txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_raw), np.std(lev_raw)))
#         txt_file.write('\n')
#         txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_raw[0], avg_err_raw[1], avg_err_raw[2]))
#         txt_file.write('\n')
#         txt_file.write('\n')
#         txt_file.write('Event: \n')
#         txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_ev), np.std(lev_ev)))
#         txt_file.write('\n')
#         txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_ev[0], avg_err_ev[1], avg_err_ev[2]))
#         txt_file.write('\n')
#         txt_file.close() 

#     vit_mean.append(np.mean(lev_vit))
#     guppy_mean.append(np.mean(lev_guppy)) 
#     ev_mean.append(np.mean(lev_ev))

#     vit_std.append(np.std(lev_vit))
#     guppy_std.append(np.std(lev_guppy)) 
#     ev_std.append(np.std(lev_ev))

#     vit_indelsub.append(avg_err_vit)  
#     guppy_indelsub.append(avg_err_guppy)  
#     ev_indelsub.append(avg_err_ev) 

# thresh = [int(d)//100 for d in deltas]
# fig, ax = plt.subplots()
# ax.set_title('Average Levenshtein Distance')
# ax.set_ylabel('Edit Distance')
# ax.set_xlabel(r'$\delta$')
# l1 = ax.plot(thresh, [x for x in vit_mean], 'o-', color="blue", label='Viterbi')
# l1 = ax.plot(thresh, [x for x in guppy_mean], '*-', color="orange", label='Guppy')
# l1 = ax.plot(thresh, [x for x in ev_mean], '^-', color="green", label='Event')

# plt.legend()
# plt.tight_layout()

# plt.savefig('/home/work/mnt/kw5km/Research/Nanopore_Channel/Output_new/deBruijn_6mer_mean/mean_lev.eps')
# plt.close()
# np.savetxt('/home/work/mnt/kw5km/Research/Nanopore_Channel/Output_new/deBruijn_6mer_mean/means.txt', np.array([vit_mean, guppy_mean, ev_mean]))
# np.savetxt('/home/work/mnt/kw5km/Research/Nanopore_Channel/Output_new/deBruijn_6mer_mean/stds.txt', np.array([vit_std, guppy_std, ev_std]))

# with open('/home/work/mnt/kw5km/Research/Nanopore_Channel/Output_new/deBruijn_6mer_mean/errs_vit.txt', 'w') as f:
#     f.write(', '.join(str(x) for x in vit_indelsub))
# f.close()
# with open('/home/work/mnt/kw5km/Research/Nanopore_Channel/Output_new/deBruijn_6mer_mean/errs_guppy.txt', 'w') as f:
#     f.write(', '.join(str(x) for x in guppy_indelsub))
# f.close()
# with open('/home/work/mnt/kw5km/Research/Nanopore_Channel/Output_new/deBruijn_6mer_mean/errs_ev.txt', 'w') as f:
#     f.write(', '.join(str(x) for x in ev_indelsub))
# f.close()


# fig, ax = plt.subplots()
# ax.set_title('Error Types')
# ax.set_ylabel('Number of Errors')
# ax.set_xlabel(r'$\delta$')
# l1 = ax.plot(thresh, [x[0] for x in vit_indelsub], 'o-', color="blue")#, label='Insertions')
# l1 = ax.plot(thresh, [x[1] for x in vit_indelsub], 'o--', color="blue")#, label='Deletions')
# l1 = ax.plot(thresh, [x[2] for x in vit_indelsub], 'o-.', color="blue")#, label='Substitutions')

# l1 = ax.plot(thresh, [x[0] for x in guppy_indelsub], '*-', color="orange")#, label='Insertions')
# l1 = ax.plot(thresh, [x[1] for x in guppy_indelsub], '*--', color="orange")#, label='Deletions')
# l1 = ax.plot(thresh, [x[2] for x in guppy_indelsub], '*-.', color="orange")#, label='Subsitutions')

# l1 = ax.plot(thresh, [x[0] for x in ev_indelsub], '^-', color="green")#, label='Insertions')
# l1 = ax.plot(thresh, [x[1] for x in ev_indelsub], '^--', color="green")#, label='Deletions')
# l1 = ax.plot(thresh, [x[2] for x in ev_indelsub], '^-.', color="green")#, label='Subsitutions')


# f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]

# lines, markers = ['-', '--', '-.'], ['o', '*', '^']
# handles = [Line2D([0], [0], color='blue', marker='o', label='Viterbi'), Line2D([0], [0], color='orange', marker='*', label='Guppy'), Line2D([0], [0], color='green', marker='^', label='Event'), 
#            Line2D([0], [0], color='k', linestyle='-', label='Insertions'), Line2D([0], [0], color='k', linestyle='--', label='Insertions'),Line2D([0], [0], color='k', linestyle='-.', label='Insertions')]

# labels = ['Viterbi', 'Guppy', 'Event', "Insertions", "Deletions", "Substitutions"]

# plt.legend(handles, labels, loc=3, framealpha=1)

# plt.tight_layout()
# # plt.show()
# plt.savefig('/home/work/mnt/kw5km/Research/Nanopore_Channel/Output_new/deBruijn_6mer_mean/errors_lev.eps')
# plt.close()   

        
    
