# env: py27
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

n, bp, ver = 1000, 200, '2'

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
        path_var_ref = '/home/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/deBruijn_6mer_mean/n{}_bp{}_delta{}_DeepSimu/fasta_raw/'.format(n, bp, delt)
        basecall_file = 'n{}_bp{}_delta{}@v{}.fastq'.format(n, bp, delt, ver)

    else:
        path_var_ref = '/home/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/deBruijn_6mer_mean/Encoded/n{}_bp{}_delta{}_statesplit_rate{}_DeepSimu/fasta_raw/'.format(n, bp, delt, r)
        basecall_file = 'n{}_bp{}_delta{}_rate{}@v{}.fastq'.format(n, bp, delt, r, ver)

    score_f, fig_f = '/scores/', '/Figs/'
    path_var_long = '/home/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/bonito/basecalls/'
    path_var_out = path_var_long  

    bon_calls = {}
    filename_bon = path_var_long+basecall_file  
    with open(filename_bon) as bon_file:
        for record in SeqIO.parse(bon_file, "fastq"):
            bon_calls[record.id] = record.seq
            # print(record.id, record.seq)
        bon_file.close() 
    print('bonito read ids:', bon_calls.keys())
    
#     mapping = {}
#     for filename in os.listdir(path_var_long+'/fast5/'):
#         filename_list = filename.split('_')
#         if 'signal' not in filename_list:continue
#         mapping[filename_list[-2]] = filename_list[-1].replace(".fast5", "")

#     cs = guppy_calls.values()

#     scores_dwll, scores_ev , scores_guppy , scores_raw = [], [], [], []
#     lev_dwll, lev_ev , lev_guppy , lev_raw, lev_vit = [], [], [], [], []
#     errs_vit, errs_ev , errs_guppy , errs_raw = [], [], [], []

#     FP_all_vit, FN_all_vit = np.array([]), np.array([])
#     FP_all_ont, FN_all_ont = np.array([]), np.array([])


#     for j in mapping.keys():

#         # print(j)
#         ''' Guppy '''
#         mapping_id = mapping[j]

#         if mapping_id in guppy_calls.keys(): 
#             seq_guppy = guppy_calls[mapping_id]
#         else: 
#             seq_guppy = cs[0]  	
#         j = int(j)

#         ''' True '''
#         filename_true = path_var_long+'processed_genome_{}'.format(j)  
#         with open(filename_true) as true_file:
#             for line in true_file:
#                 if '>' not in line: bases_true = line
#             true_file.close()
#         seq_true = Seq(bases_true)

#         ''' Viterbi Dwell '''
#         filename_viterbi_dwll = path_var_long_out+calls_f+'j{}'.format(j)
#         # if not Path(filename_viterbi_dwll).is_file(): continue
#         if Path(filename_viterbi_dwll).is_file():
#             with open(filename_viterbi_dwll) as viterbi_file_dwll:
#                 for line in viterbi_file_dwll:
#                     viterbi_kmers_dwll = line
#                 viterbi_file_dwll.close()
            
#             ch = ['\'', '[', ']', '2', '3', ' ']
#             # Remove multiple characters from the string
#             for character in ch:
#                 viterbi_kmers_dwll = viterbi_kmers_dwll.replace(character, '')  
            
#             viterbi_kmers_dwll = viterbi_kmers_dwll.split(',')
                
#             viterbi_seq_dwll = "" + viterbi_kmers_dwll[0]
#             last_kmer = viterbi_kmers_dwll[0]
#             viterbi_kmers_proc = [last_kmer]
#             for b in viterbi_kmers_dwll[1:]:
#                 if b == last_kmer: continue
#                 else: 
#                     viterbi_kmers_proc.append(b)
#                     viterbi_seq_dwll=viterbi_seq_dwll+b[-1]
#                     last_kmer=b
#                 seq_dwll = Seq(viterbi_seq_dwll)

#             file_out = path_var_long_out+calls_f+'j{}_processed'.format(j)
#             with open(file_out, "w") as txt_file:
#                     txt_file.write(str(viterbi_kmers_proc))
#                     txt_file.close()

#         ''' Events '''
#         filename_ev = path_var_long + 'fasta/signal_{}_{}.fasta'.format(j, mapping_id)  
#         with open(filename_ev) as ev_file:
#             for line in ev_file:
#                 if '>' not in line: bases_ev = line
#             ev_file.close()
#         seq_ev = Seq(bases_ev)

#         ''' Raw '''       
#         filename_raw = path_var_long + 'fasta_raw/signal_{}_{}.fasta'.format(j, mapping_id)  
#         if not Path(filename_raw).is_file(): continue
#         with open(filename_raw) as raw_file:
#             for line in raw_file:
#                 if '>' not in line: bases_raw = line
#             raw_file.close()
#         seq_raw = Seq(bases_raw)

#         ''' Segs True'''
#         # tru_seg = [0]
#         # last = 0
#         # alig = np.loadtxt(path_var_long+'align/align_{}.ali'.format(j))
#         # for index in alig:
#         #     if index[1] == last:
#         #         continue
#         #     tru_seg.append(index[0]-1)
#         #     last = index[1]
#         # df_b = np.array(tru_seg)
#         # df_b = np.arange(bp)

#         ''' Segs Event'''
#         # filename_hdf5 = path_var_long+'hdf5/signal_{}_{}.hdf5'.format(j, mapping_id)
#         # new_dat = read_dump_seq(filename_hdf5)
            
#         # ## ONT Segmentation
#         # segs_ont=new_dat[:, 0]
#         # db_ont, fn_ont, FN_ont, FP_ont = score_seg(df_b, segs_ont)

#         ''' Segs Viterbi'''
#         ## Get Viterbi estimated kmers
#         # filename_figs_base = path_var_long_out+fig_f
                
#         # viterbi_kmers = viterbi_kmers_dwll
        
#         ## Find the segments inferred from the Viterbi estimation kmers
#         # segs_viterbi = [0]
#         # for i, b in enumerate(viterbi_kmers[1:], 1):
#         #     if b == viterbi_kmers[i-1]: 
#         #         continue
#         #     else: segs_viterbi.append(i)
#         # db_viterbi, fn_viterbi, FN_viterbi, FP_viterbi = score_seg(df_b, segs_viterbi)

#         # ''' Save Segs'''
#         # # FP_all_ont = np.append(FP_all_ont, FP_ont)
#         # # FN_all_ont = np.append(FN_all_ont, FN_ont)
#         # FP_all_vit = np.append(FP_all_vit, FP_viterbi)
#         # FN_all_vit = np.append(FN_all_vit, FN_viterbi)
        
#         ''' Align '''
#         if Path(filename_viterbi_dwll).is_file():
#             vit_align = edlib.align(str(seq_dwll), str(seq_true), task = "path")
#             vit_alignment = edlib.getNiceAlignment(vit_align, str(seq_dwll), str(seq_true))
#             lev_vit.append(vit_align["editDistance"])
#             vit_ins, vit_dels, vit_subs = format_cigar(vit_align["cigar"])
#             errs_vit.append((vit_ins, vit_dels, vit_subs))

#         guppy_align = edlib.align(str(seq_guppy), str(seq_true), task = "path")
#         guppy_alignment = edlib.getNiceAlignment(guppy_align, str(seq_guppy), str(seq_true))
#         lev_guppy.append(guppy_align["editDistance"])
#         guppy_ins, guppy_dels, guppy_subs = format_cigar(guppy_align['cigar'])
#         errs_guppy.append((guppy_ins, guppy_dels, guppy_subs))

#         ev_align = edlib.align(str(seq_ev), str(seq_true), task = "path")
#         ev_alignment = edlib.getNiceAlignment(ev_align, str(seq_ev), str(seq_true))
#         lev_ev.append(ev_align["editDistance"])
#         ev_ins, ev_dels, ev_subs = format_cigar(ev_align["cigar"])
#         errs_ev.append((ev_ins, ev_dels, ev_subs))

#         raw_align = edlib.align(str(seq_raw), str(seq_true), task = "path")
#         raw_alignment = edlib.getNiceAlignment(raw_align, str(seq_raw), str(seq_true))
#         lev_raw.append(raw_align["editDistance"])
#         raw_ins, raw_dels, raw_subs = format_cigar(raw_align["cigar"])
#         errs_raw.append((raw_ins, raw_dels, raw_subs))

#         filepath_align = path_var_long_out+score_f
#         if not os.path.exists(filepath_align):
#             os.makedirs(filepath_align)

#         filename_align = filepath_align+'align_j{}.txt'.format(j)
#         with open(filename_align, "w") as txt_file:
#             txt_file.write('(top seq<-true) \n')
#             if Path(filename_viterbi_dwll).is_file():
#                 txt_file.write('Viterbi: \n')
#                 txt_file.write("\n".join(vit_alignment.values())) 
#                 txt_file.write('\n')
#                 txt_file.write('Lev Dist: {}'.format(vit_align["editDistance"]))
#                 txt_file.write('\n')
#                 txt_file.write('Edits: {}'.format(format_cigar(vit_align["cigar"])))
#                 txt_file.write('\n')
#                 txt_file.write('\n')
#             txt_file.write('Guppy: \n')
#             txt_file.write("\n".join(guppy_alignment.values())) 
#             txt_file.write('\n')
#             txt_file.write('Lev Dist: {}'.format(guppy_align["editDistance"]))
#             txt_file.write('\n')
#             txt_file.write('Edits: {}'.format(guppy_align["cigar"]))
#             txt_file.write('\n')
#             txt_file.write('Insertions: {}, Deletions: {}, Substitutions: {}'.format(guppy_ins, guppy_dels, guppy_subs))
#             txt_file.write('\n')
#             txt_file.write('\n')
#             txt_file.write('Raw: \n')
#             txt_file.write("\n".join(raw_alignment.values()))  
#             txt_file.write('\n')
#             txt_file.write('Lev Dist: {}'.format(raw_align["editDistance"]))
#             txt_file.write('\n')
#             txt_file.write('Edits: {}'.format(raw_align["cigar"]))
#             txt_file.write('\n')
#             txt_file.write('Insertions: {}, Deletions: {}, Substitutions: {}'.format(raw_ins, raw_dels, raw_subs))
#             txt_file.write('\n')
#             txt_file.write('\n')
#             txt_file.write('Event: \n')
#             txt_file.write("\n".join(ev_alignment.values())) 
#             txt_file.write('\n')
#             txt_file.write('Lev Dist: {}'.format(ev_align["editDistance"])) 
#             txt_file.write('\n')
#             txt_file.write('Edits: {}'.format(ev_align["cigar"]))
#             txt_file.write('\n')
#             txt_file.write('\n')
#             txt_file.write('Insertions: {}, Deletions: {}, Substitutions: {}'.format(ev_ins, ev_dels, ev_subs))
#             txt_file.close()

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

        
    
