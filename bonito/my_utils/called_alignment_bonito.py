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

n, bp, ver = 1000, 200, '2'

bon_mean = []
bon_std = []
bon_indelsub = []

rates = ["106", "106", "106", "106", "106", "106", "86", "86", "86", "86", "66"]
deltas = ['000', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1000']

score_f, fig_f = 'scores/', 'Figs/'
path_var_long = '/home/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/bonito/basecalls/'

for delt, r in zip(deltas, rates):
    
    print('DELTA:', delt)
    cwd = os.getcwd()

    # 0 is entirely unconstrained so no state split
    if delt == '000':
        path_var_ref = '/home/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/deBruijn_6mer_mean/n{}_bp{}_delta{}_DeepSimu/fasta/'.format(n, bp, delt)
        basecall_file = 'n{}_bp{}_delta{}@v{}.fastq'.format(n, bp, delt, ver)
        f_out = 'n{}_bp{}_delta{}@v{}/'.format(n, bp, delt, ver)

    else:
        path_var_ref = '/home/work/mnt/kw5km/Research/Nanopore_Channel/Data_new/deBruijn_6mer_mean/Encoded/n{}_bp{}_delta{}_statesplit_rate{}_DeepSimu/fasta/'.format(n, bp, delt, r)
        basecall_file = 'n{}_bp{}_delta{}_rate{}@v{}.fastq'.format(n, bp, delt, r, ver)
        f_out = 'n{}_bp{}_delta{}_rate{}@v{}/'.format(n, bp, delt, r, ver)

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

        bon_align = edlib.align(str(seq_bon), str(seq_true), task = "path", mode="HW")
        bon_alignment = edlib.getNiceAlignment(bon_align, str(seq_bon), str(seq_true))
        lev_bon.append(bon_align["editDistance"])
        bon_ins, bon_dels, bon_subs = format_cigar(bon_align['cigar'])
        errs_bon.append((bon_ins, bon_dels, bon_subs))


        filepath_align = path_var_long_out+score_f
        if not os.path.exists(filepath_align):
            os.makedirs(filepath_align)

        filename_align = filepath_align+'align_{}.txt'.format(signal_id)
        with open(filename_align, "w") as txt_file:
            txt_file.write('(bottom seq<-true) \n')
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


    avg_err_bon = np.mean(errs_bon, axis=0)

    filename_scores = filepath_align+'alignment_scores.txt'
    with open(filename_scores, "w") as txt_file:
        txt_file.write('Bonito: \n')
        txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_bon), np.std(lev_bon)))
        txt_file.write('\n')
        txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_bon[0], avg_err_bon[1], avg_err_bon[2]))
        txt_file.write('\n')
        txt_file.close() 

    bon_mean.append(np.mean(lev_bon))
    bon_std.append(np.std(lev_bon))
    bon_indelsub.append(avg_err_bon)  

filepath_figs = path_var_long+fig_f
if not os.path.exists(filepath_figs):
            os.makedirs(filepath_figs)
            
thresh = [int(d)//100 for d in deltas]
fig, ax = plt.subplots()
ax.set_title('Average Levenshtein Distance')
ax.set_ylabel('Edit Distance')
ax.set_xlabel(r'$\delta$')
l1 = ax.plot(thresh, [x for x in bon_mean], 'o-', color="red", label='Bonito')

plt.legend()
plt.tight_layout()
plt.savefig(path_var_long+fig_f+'mean_lev_@{}.eps'.format(ver))
plt.close()

np.savetxt(path_var_long+'means_@{}.txt'.format(ver), np.array([bon_mean]))
np.savetxt(path_var_long+'stds_@{}.txt'.format(ver), np.array([bon_std]))

with open(path_var_long+'errs_@{}.txt'.format(ver), 'w') as f:
    f.write(', '.join(str(x) for x in bon_indelsub))
f.close()


fig, ax = plt.subplots()
ax.set_title('Error Types')
ax.set_ylabel('Number of Errors')
ax.set_xlabel(r'$\delta$')
l1 = ax.plot(thresh, [x[0] for x in bon_indelsub], 'o-', color="red", label='Insertions')
l1 = ax.plot(thresh, [x[1] for x in bon_indelsub], 'o--', color="red", label='Deletions')
l1 = ax.plot(thresh, [x[2] for x in bon_indelsub], 'o-.', color="red", label='Substitutions')


plt.legend()

plt.tight_layout()
# plt.show()
plt.savefig(path_var_long+fig_f+'errors_lev_@{}.eps'.format(ver))
plt.close()   

        
    
