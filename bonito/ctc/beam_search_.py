from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, List, Tuple

import numpy as np
import itertools

# from ctc_decoder.language_model import LanguageModel


def log(x: float) -> float:
    with np.errstate(divide='ignore'):
        return np.log(x)

def get_all_kmers(k, alphabet):
    l = [''.join(i) for i in itertools.product(alphabet, repeat = k)]
    return np.array(l)


def get_debruijn_edges_from_kmers(kmers):
    """
    Every possible (k-1)mer (n-1 suffix and prefix of kmers) is assigned
    to a node, and we connect one node to another if the (k-1)mer overlaps 
    another. Nodes are (k-1)mers, edges are kmers.
    """
    # store edges as tuples in a set
    edges = np.zeros((len(kmers), len(kmers)))
    edges_set = set()
    
    # compare each (k-1)mer
    for i, k1 in enumerate(kmers):
        for j, k2 in enumerate(kmers):
            if k1 != k2: 
                # if they overlap then add to edges
                if k1[1:] == k2[:-1]:
#                    print(k1, k2)
                    edges[i, j] = 1
                    edges_set.add((k1[:-1], k2[:-1]))

    return edges, edges_set

def get_scores_m(edges, kmers, score_file='/u/kw5km/Research/Nanopore_Channel/Regression_Output/6mer_means.txt'):
    """ Get the de Bruijn edge score for each kmers transition """
    means = np.loadtxt(score_file, skiprows=1, usecols=[0,1], dtype=str)

    score_vals = {k: v for k, v in zip(means[:, 0], means[:, 1].astype(float))}

    scores = np.zeros_like(edges)
    for node1 in range(edges.shape[0]):
        for node2 in range(edges.shape[1]):

            if edges[node1, node2] == 1:
                k1, k2 = kmers[node1], kmers[node2]
                val1 = score_vals[k1]
                val2 = score_vals[k2]
                scores[node1, node2] = np.abs(val1-val2)
    return scores, score_vals

def remove_nodes(edges):
    loop = True
    
    last_sums = np.zeros_like((np.sum(edges, axis=1)))
    while loop==True:
        
        edge_sums = np.sum(edges, axis = 1)
        if np.array_equal(edge_sums, last_sums): loop=False
        
        idxs_bad = [edge_sums<1][0]
        edges[:, idxs_bad] = 0
        
        last_sums=edge_sums

    return edges

def get_cs(beam, chars, kmers, edges):
    if len(list(beam.labeling))<6:
        return np.arange(1,len(chars))
    k_char_idxs = list(beam.labeling[-6:])
    k = [chars[k] for k in k_char_idxs]
    k = ''.join(k)
    k_idx = np.where(kmers==k)[0][0]
    # print('Labeling chars:', k_char_idxs, 'kmer:', k, 'kmer index:', k_idx)

    cs = np.where(edges[k_idx, :]>0)[0]
    k_cs = kmers[cs]

    k_thresh = np.unique([i[-1] for i in k_cs])
    c_thresh = [chars.index(i) for i in k_thresh]#[np.where(i==chars)[0] for i in k_thresh]
    # print('Transition idxs:', cs, 'Transition kmers:', k_cs, 'Transition char:', k_thresh, 'Transition cahr idx:', c_thresh)

    return c_thresh

@dataclass
class BeamEntry:
    """Information about one single beam at specific time-step."""
    pr_total: float = log(0)  # blank and non-blank
    pr_non_blank: float = log(0)  # non-blank
    pr_blank: float = log(0)  # blank
    pr_text: float = log(1)  # LM score
    lm_applied: bool = False  # flag if LM was already applied to this beam
    labeling: tuple = ()  # beam-labeling


class BeamList:
    """Information about all beams at specific time-step."""

    def __init__(self) -> None:
        self.entries = defaultdict(BeamEntry)

    def normalize(self) -> None:
        """Length-normalise LM score."""
        for k in self.entries.keys():
            labeling_len = len(self.entries[k].labeling)
            self.entries[k].pr_text = (1.0 / (labeling_len if labeling_len else 1.0)) * self.entries[k].pr_text

    def sort_labelings(self) -> List[Tuple[int]]:
        """Return beam-labelings, sorted by probability."""
        beams = self.entries.values()
        sorted_beams = sorted(beams, reverse=True, key=lambda x: x.pr_total + x.pr_text)
        return [x.labeling for x in sorted_beams]


# def apply_lm(parent_beam: BeamEntry, child_beam: BeamEntry, chars: str, lm: LanguageModel) -> None:
#     """Calculate LM score of child beam by taking score from parent beam and bigram probability of last two chars."""
#     if not lm or child_beam.lm_applied:
#         return

#     # take bigram if beam length at least 2
#     if len(child_beam.labeling) > 1:
#         c = chars[child_beam.labeling[-2]]
#         d = chars[child_beam.labeling[-1]]
#         ngram_prob = lm.get_char_bigram(c, d)
#     # otherwise take unigram
#     else:
#         c = chars[child_beam.labeling[-1]]
#         ngram_prob = lm.get_char_unigram(c)

#     lm_factor = 0.01  # influence of language model
#     child_beam.pr_text = parent_beam.pr_text + lm_factor * log(ngram_prob)  # probability of char sequence
#     child_beam.lm_applied = True  # only apply LM once per beam entry

def apply_lm(child_beam: BeamEntry, chars: str, edges, kmers) -> None:
    """Calculate LM score of child beam by taking score from parent beam and bigram probability of last two chars."""
    if child_beam.lm_applied or len(child_beam.labeling)<=6:
        return

    # Get kmers from seq and kmer scores
    k1_idxs = list(child_beam.labeling[-7:-1])
    k2_idxs = list(child_beam.labeling[-6:])

    k1 = [chars[k] for k in k1_idxs]
    k1 = ''.join(k1)
    k2 = [chars[k] for k in k2_idxs]
    k2 = ''.join(k2)
    
    
    k1_idx = np.where(kmers==k1)[0][0]
    k2_idx = np.where(kmers==k2)[0][0]
    # print('Labeling chars:', k_char_idxs, 'kmer:', k, 'kmer index:', k_idx)

    # not exactly 1 or 0 for the log
    if edges[k1_idx, k2_idx]>0: prob = 1 - 1e-15
    else: prob = 1e-15

    child_beam.pr_text = child_beam.pr_text + (0.2)*log(prob)  # probability of char sequence
    child_beam.lm_applied = True  # only apply LM once per beam entry

def beam_search(mat: np.ndarray, chars: str, beam_width: int = 25, threshold=1e-4, delta=0) -> str:
    """Beam search decoder.

    See the paper of Hwang et al. and the paper of Graves et al.

    Args:
        mat: Output of neural network of shape TxC.
        chars: The set of characters the neural network can recognize, excluding the CTC-blank.
        beam_width: Number of beams kept per iteration.
        lm: Character level language model if specified.

    Returns:
        The decoded text.
    """

    kmers = get_all_kmers(6, ['A', 'C', 'G', 'T'])
    edges, edges_set = get_debruijn_edges_from_kmers(kmers)
    scores, score_vals = get_scores_m(edges, kmers, score_file='/u/kw5km/Research/Nanopore_Channel/Regression_Output/6mer_means.txt')
    mask1 = [scores<delta][0]
    edges[mask1] = 0
    edges = remove_nodes(edges)

    blank_idx = 0 #len(chars) #blank is at beginning of alphabet
    max_T, max_C = mat.shape

    # initialise beam state
    last = BeamList()
    labeling = ()
    last.entries[labeling] = BeamEntry()
    last.entries[labeling].pr_blank = log(1)
    last.entries[labeling].pr_total = log(1)

    # go over all time-steps
    for t in range(max_T):
        curr = BeamList()

        # get beam-labelings of best beams
        best_labelings = last.sort_labelings()[:beam_width]

        # go over best beams
        for labeling in best_labelings:

            # probability of paths ending with a non-blank
            pr_non_blank = log(0)
            # in case of non-empty beam
            if labeling:
                # probability of paths with repeated last char at the end
                pr_non_blank = last.entries[labeling].pr_non_blank + log(mat[t, labeling[-1]])

            # probability of paths ending with a blank
            pr_blank = last.entries[labeling].pr_total + log(mat[t, blank_idx])

            # fill in data for current beam
            curr.entries[labeling].labeling = labeling
            curr.entries[labeling].pr_non_blank = np.logaddexp(curr.entries[labeling].pr_non_blank, pr_non_blank)
            curr.entries[labeling].pr_blank = np.logaddexp(curr.entries[labeling].pr_blank, pr_blank)
            curr.entries[labeling].pr_total = np.logaddexp(curr.entries[labeling].pr_total,
                                                           np.logaddexp(pr_blank, pr_non_blank))
            curr.entries[labeling].pr_text = last.entries[labeling].pr_text
            curr.entries[labeling].lm_applied = True  # LM already applied at previous time-step for this beam-labeling

            # extend current beam-labeling
            c_scored = get_cs(curr.entries[labeling], chars, kmers, edges)
            for c in c_scored:#range(1, max_C):
                # add new char to current beam-labeling
                new_labeling = labeling + (c,)

                # if new labeling contains duplicate char at the end, only consider paths ending with a blank
                if labeling and labeling[-1] == c:
                    pr_non_blank = last.entries[labeling].pr_blank + log(mat[t, c])
                else:
                    pr_non_blank = last.entries[labeling].pr_total + log(mat[t, c])

                # fill in data
                curr.entries[new_labeling].labeling = new_labeling
                curr.entries[new_labeling].pr_non_blank = np.logaddexp(curr.entries[new_labeling].pr_non_blank,
                                                                       pr_non_blank)
                curr.entries[new_labeling].pr_total = np.logaddexp(curr.entries[new_labeling].pr_total, pr_non_blank)

                # apply LM
                # apply_lm(curr.entries[new_labeling], chars, score_vals, delta=delta)

        # set new beam state
        last = curr

    # normalise LM scores according to beam-labeling-length
    last.normalize()

    # sort by probability
    best_labeling = last.sort_labelings()[0]  # get most probable labeling

    # map label string to char string
    res = ''.join([chars[label] for label in best_labeling])
    return res
