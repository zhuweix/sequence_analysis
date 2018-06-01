# -*- coding: utf-8 -*-
"""Sequence operations."""
import re
import time
import os
import numpy as np
import matplotlib.pyplot as plt
# from .biofile import simple_gff3_load
# from .biofile import simple_fasta_write


def _message(mess=""):
    """Print a message with time.
    Args:
        mess: message to be printed

    """
    time_str = time.strftime("%Y-%m-%d %H:%M:%S")
    print ("%s: %s" %(time_str, mess))

def rc_seq(seq=""):
    """Returns the reverse compliment sequence."""
    rc_nt_ls = []
    rc_dict = {
        "a": "t",
        "c": "g",
        "t": "a",
        "g": "c",
        "n": "n",
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "N": "N"
    }
    rc_nt_ls = [rc_dict[seq[i]] for i in range(len(seq)-1, -1, -1)]
    rc_seq_ = "".join(rc_nt_ls)
    return rc_seq_

def translate_exon(seq):
    """Translate exon by normal codon."""

    codon = {
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "CTT": "L",
        "CTC": "L",
        "CTG": "L",
        "CTA": "L",
        "TTA": "L",
        "TTG": "L",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "TTT": "F",
        "TTC": "F",
        "ATG": "M",
        "TGT": "C",
        "TGC": "C",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "ACT": "T",
        "ACC": "T",
        "ACG": "T",
        "ACA": "T",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "AGT": "S",
        "AGC": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TGG": "W",
        "CAA": "Q",
        "CAG": "Q",
        "AAT": "N",
        "AAC": "N",
        "CAT": "H",
        "CAC": "H",
        "GAA": "E",
        "GAG": "E",
        "GAT": "D",
        "GAC": "D",
        "AAA": "K",
        "AAG": "K",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGA": "R",
        "AGG": "R",
        "TAA": "*",
        "TAG": "*",
        "TGA": "*"
    }

    seq = seq.upper()
    size = len(seq)
    if size % 3 != 0:
        _message("Length is %d, invalid !" % size)
    protein = [codon[seq[id_: id_+3]] for id_ in range(0, size, 3)]
    protein = "".join(protein)
    return protein

# def extract_cdna(dirn: str, gff:str, prefix="out") -> tuple:
    """Extract cdna seqence."""
    os.chdir(dirn)
    gff_entries, chrom, seqs = simple_gff3_load(gff, return_fasta=True)
    chrom = [name.split()[0] for name in chrom]
    cdna_names = []
    cdna_seqs = []
    cdna_location = {}
    for entry in gff_entries:
        name = entry[8]
        name = name.replace("ID=", "")
        name = name.replace("-mRNA-1", "")
        if entry[2]  != "gene":
            continue
        cdna_names.append(name)
        chrom_id = chrom.index(entry[0])
        seq = seqs[chrom_id][entry[3] -1: entry[4]]
        if entry[6] == "-":
            seq = rc_seq(seq)
        cdna_seqs.append(seq)
        cdna_location.setdefault(entry[0], [])
        cdna_location[entry[0]].append((entry[3], entry[4], entry[8]))

    # for entry in gff_entries:
    #     name = entry[8]
    #     name.replace("-mRNA-1", "")
    #     if entry[2]  != "mRNA":
    #         continue
    #     if name in cdna_names:
    #         continue
    #     cda_name.append[name]
    #     chrom_id = chrom.index(entry[0])
    #     seq = seqs[chrom_id][entry[3] -1: entry[4]]
    #     if entry[6] == "-":
    #         seq = rc_seq(seq)
    #     cdna_seqs.append(seq)
    #     cdna_location.setdefault(entry[0], [])
    #     cdna_location[entry[0]].append((entry[3], entry[4], entry[8]))
    for chrom in cdna_location.keys():
        genes = cdna_location[chrom]
        genes.sort()
        genes = [(id_, gene[0], gene[1], gene[2]) for id_, gene in enumerate(genes, 1)]
        cdna_location[chrom] = genes

    simple_fasta_write(prefix + ".fa", cdna_names, cdna_seqs)
    return (cdna_names, cdna_seqs, cdna_location)

def load_blastp_score(score: str) -> dict:
    """Generate score matrix."""
    aa_list = []
    score_dict = {}
    with open(score, "r") as filep:
        # Omint comment
        for line in filep:
            if line[0] == "#":
                continue
            break
        # Load a.a
        aa_list = list(line.split())
        for line in filep:
            entry = line.split()
            for _id, _aa in enumerate(aa_list):
                score_dict[(entry[0], _aa)] = entry[_id + 1]
    return score_dict

def search_orf(seq:str, min_orf:int) -> list:
    """Search full orf over ceration length in 6 frames"""
    scod = "M"
    send = "*"
    orf_regions = {}
    # Load 6 reading frames
    seq1 = seq
    seq2 = seq1[1: ]
    seq3 = seq1[2: ]
    seq4 = rc_seq(seq1)
    seq5 = seq4[1: ]
    seq6 = seq4[2: ]
    # Shrink to times of 3
    seq1 = seq1[: len(seq1)//3*3]
    seq2 = seq2[: len(seq2)//3*3]
    seq3 = seq3[: len(seq3)//3*3]
    seq4 = seq4[: len(seq4)//3*3]
    seq5 = seq5[: len(seq5)//3*3]
    seq6 = seq6[: len(seq6)//3*3]
    # Translate 6 frames
    trans1 = translate_exon(seq1)
    trans2 = translate_exon(seq2)
    trans3 = translate_exon(seq3)
    trans4 = translate_exon(seq4)
    trans5 = translate_exon(seq5)
    trans6 = translate_exon(seq6)
    # All the start and stop codons
    start1 = [id_ for id_, cod in enumerate(trans1) if cod == scod]
    start2 = [id_ for id_, cod in enumerate(trans2) if cod == scod]
    start3 = [id_ for id_, cod in enumerate(trans3) if cod == scod]
    start4 = [id_ for id_, cod in enumerate(trans4) if cod == scod]
    start5 = [id_ for id_, cod in enumerate(trans5) if cod == scod]
    start6 = [id_ for id_, cod in enumerate(trans6) if cod == scod]
    end1 = [id_ for id_, cod in enumerate(trans1) if cod == send]
    end2 = [id_ for id_, cod in enumerate(trans2) if cod == send]
    end3 = [id_ for id_, cod in enumerate(trans3) if cod == send]
    end4 = [id_ for id_, cod in enumerate(trans4) if cod == send]
    end5 = [id_ for id_, cod in enumerate(trans5) if cod == send]
    end6 = [id_ for id_, cod in enumerate(trans6) if cod == send]


    if start1 and end1:
        pos1 = start1[0]
        pos2 = end1[0]
        s_i = 0
        e_i = 0
        while s_i < len(start1):
            # search for stop codon
            pos1 = start1[s_i]
            while e_i < len(end1) -1:
                if pos2 < pos1:
                    e_i += 1
                    pos2 = end1[e_i]
                else:
                    break
            # No full orf
            if pos2 < pos1:
                break
            # Identify orf
            else:
                if pos2 - pos1 >min_orf:
                    orf_regions[(pos1*3 + 1, pos2*3 + 3)] = trans1[pos1: pos2 +1]
                s_i += 1
                while s_i < len(start1):
                    pos1 = start1[s_i]
                    if pos1 < pos2:
                        s_i += 1
                        continue
                    break
    if start2 and end2:
        pos1 = start2[0]
        pos2 = end2[0]
        s_i = 0
        e_i = 0
        while s_i < len(start2):
            # search for stop codon
            pos1 = start2[s_i]
            while e_i < len(end2) -1:
                if pos2 < pos1:
                    e_i += 1
                    pos2 = end2[e_i]
                else:
                    break
            # No full orf
            if pos2 < pos1:
                break
            # Identify orf
            else:
                if pos2 - pos1 >min_orf:
                    orf_regions[(pos1*3 + 2, pos2*3 + 4)] = trans2[pos1: pos2 +1]
                s_i += 1
                while s_i < len(start2):
                    pos1 = start2[s_i]
                    if pos1 < pos2:
                        s_i += 1
                        continue
                    break
    if start3 and end3:
        pos1 = start3[0]
        pos2 = end3[0]
        s_i = 0
        e_i = 0
        while s_i < len(start3):
            # search for stop codon
            pos1 = start3[s_i]
            while e_i < len(end3) -1:
                if pos2 < pos1:
                    e_i += 1
                    pos2 = end3[e_i]
                else:
                    break
            # No full orf
            if pos2 < pos1:
                break
            # Identify orf
            else:
                if pos2 - pos1 >min_orf:
                    orf_regions[(pos1*3 + 3, pos2*3 + 5)] = trans3[pos1: pos2 +1]
                s_i += 1
                while s_i < len(start3):
                    pos1 = start3[s_i]
                    if pos1 < pos2:
                        s_i += 1
                        continue
                    break
    if start4 and end4:
        pos1 = start4[0]
        pos2 = end4[0]
        s_i = 0
        e_i = 0
        while s_i < len(start4):
            # search for stop codon
            pos1 = start4[s_i]
            while e_i < len(end4) -1:
                if pos2 < pos1:
                    e_i += 1
                    pos2 = end4[e_i]
                else:
                    break
            # No full orf
            if pos2 < pos1:
                break
            # Identify orf
            else:
                if pos2 - pos1 >min_orf:
                    orf_regions[(len(seq) - pos1*3, len(seq) - pos2*3 -2)] = trans4[pos1: pos2 +1]
                s_i += 1
                while s_i < len(start4):
                    pos1 = start4[s_i]
                    if pos1 < pos2:
                        s_i += 1
                        continue
                    break
    if start5 and end5:
        pos1 = start5[0]
        pos2 = end5[0]
        s_i = 0
        e_i = 0
        while s_i < len(start5):
            # search for stop codon
            pos1 = start5[s_i]
            while e_i < len(end5) -1:
                if pos2 < pos1:
                    e_i += 1
                    pos2 = end5[e_i]
                else:
                    break
            # No full orf
            if pos2 < pos1:
                break
            # Identify orf
            else:
                if pos2 - pos1 >min_orf:
                    orf_regions[(len(seq) - pos1*3 -1, len(seq) - pos2*3 -3,)] = trans5[pos1: pos2 +1]
                s_i += 1
                while s_i < len(start5):
                    pos1 = start5[s_i]
                    if pos1 < pos2:
                        s_i += 1
                        continue
                    break
    if start6 and end6:
        pos1 = start6[0]
        pos2 = end6[0]
        s_i = 0
        e_i = 0
        while s_i < len(start6):
            # search for stop codon
            pos1 = start6[s_i]
            while e_i < len(end6) -1:
                if pos2 < pos1:
                    e_i += 1
                    pos2 = end6[e_i]
                else:
                    break
            # No full orf
            if pos2 < pos1:
                break
            # Identify orf
            else:
                if pos2 - pos1 >min_orf:
                    orf_regions[(len(seq) - pos1*3 -2, len(seq) - pos2*3 -4)] = trans6[pos1: pos2 +1]
                s_i += 1
                while s_i < len(start6):
                    pos1 = start6[s_i]
                    if pos1 < pos2:
                        s_i += 1
                        continue
                    break

    return orf_regions

# def prot_dotplot(seq1: str, seq2:str, length=10, score="BLOSUM62.dat"):
    """Generate protein dotplot."""
    score_dict = load_blastp_score(score)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    # Generate score
    score_m1 = np.zeros((len(seq2), len(seq1)), dtype=int)
    for id_1, aa_1 in enumerate(seq1):
        for id_2, aa_2 in enumerate(seq2):
            score_m1[id_2][id_1] = score_dict[(aa_1, aa_2)]
    dot_score = np.zeros((len(seq2) - length + 1, len(seq1) - length +1), dtype=int)
    for id_1 in range(len(seq2) - length + 1):
        for id_2 in range(len(seq1) - length + 1):
            score_ = 0
            for id_ in range(length):
                score_ += score_m1[id_1 + id_][id_2 + id_]
    return dot_score

def calc_tm(seq: str) -> int:
    """Calculate the Tm for primers.

    primer: 50 nM ; Na+: 50 mM, pH = 7
    """
    seq = seq.upper()
    length = len(seq)
    try:
        assert len(seq) > 1
    except:
        _message("The sequence: %s is too short to calculate Tm" % seq)
    gc_count = seq.count("G") + seq.count("C")
    if length <= 13:
        return (length - gc_count) * 2 + gc_count *4
    else:
        return 64.9 + 41 * (gc_count - 16.4) / length

def search_pcr_primer_pair(
        sub_seq: str, u_pos: int, d_pos: int,
        tm_max: int, tm_min: int, length_min: float, length_max: float,
        strain_seqs=None, homo_mask=4, gc_end=True, trial=10) -> list:
    """Search pcr primer.
    seq: sequence to search for primer
    pos: search primer around the position
    u_flank/ d_flank: minimum up/down stream flanking region
    length_min: minimum length of primers
    length_max: maximum length of primers
    strain_seqs (list): None if not searching for unique primer, list of seqs for unique primer
    homo_mask: min length of homopolymers to be masked
    """
    # homopolymer mask
    mask = ["A" * homo_mask, "T" * homo_mask, "C" * homo_mask, "G" * homo_mask]
    # Only allow unqiue kemr
    kmer_count = {}
    if strain_seqs:
        for seq in strain_seqs:
            seq = seq.upper()
            for pos in range(len(seq) - length_min +1):
                kmer = seq[pos: pos + length_min]
                for homo in mask:
                    if homo in kmer:
                        continue
                rc_kmer = rc_seq(kmer)
                kmer = kmer if kmer < rc_kmer else rc_kmer
                kmer_count.setdefault(kmer, 0)
                kmer_count[kmer] += 1
        allow_kmer = {}
        for kmer, count in kmer_count.items():
            if count == 1:
                allow_kmer[kmer] = None
        del kmer_count
        primer_list = []
        seq = seq.upper()
        for tm in range(tm_min, tm_max +1):
            f_primer = search_forward_pcr_primer(sub_seq, u_pos, tm + 1, tm -1,
                length_min, length_max, gc_end, allow_kmer, trial)
            r_primer = search_reverse_pcr_primer(sub_seq, d_pos, tm + 1, tm -1,
                length_min, length_max, gc_end, allow_kmer, trial)
            if f_primer and r_primer:
                primer_list.append(f_primer + r_primer)
        return primer_list

def gc_ratio(seq:str) -> int:
    """Return the GC ratio"""
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return gc_count/ len(seq)

def search_forward_pcr_primer(seq: str, pos: int,
        tm_max: float, tm_min: float, length_min: float,
        length_max: float, gc_end: bool, alllow_kmer=None, trial=10):
    """Search for forward primer"""
    seq = seq.upper()
    primer_candidates = []
    if alllow_kmer:
        for end in range(pos, length_max, -1):
            kmer = seq[end - length_min: end]
            if gc_end:
                if kmer[-1] != "C" and kmer[-1] != "G":
                    continue
            rc_kmer = rc_seq(kmer)
            u_kmer = kmer if kmer < rc_kmer else rc_kmer
            if u_kmer not in alllow_kmer:
                continue
            for length in range(length_min, length_max):
                kmer = seq[end - length: end]
                # Expand the kmer to certain Tm
                tm = calc_tm(kmer)
                if tm_min < tm < tm_max:
                    gc = gc_ratio(seq)
                    if .4 < gc < .6:
                        return kmer, tm, end - length + 1
                    else:
                        primer_candidates.append((kmer, tm, end - length + 1))
                        break
            if len(primer_candidates) > trial:
                break
        # Return the best gc ratio
        # Second range 30% - 70%
        for primer in primer_candidates:
            gc = gc_ratio(primer[0])
            if .3 < gc < .7:
                return primer
        if primer_candidates:
            return primer_candidates[0]
        else:
            print(len(seq), pos)
            return

def search_reverse_pcr_primer(seq: str, pos: int,
        tm_max: float, tm_min: float, length_min: float,
        length_max: float, gc_end: bool, alllow_kmer=None, trial=10):
    """Search reverse pcr primer"""
    seq = seq.upper()
    seq = rc_seq(seq)
    pos = len(seq) - pos
    primer = search_forward_pcr_primer(seq, pos,
        tm_max, tm_min, length_min, length_max, gc_end, alllow_kmer, trial)
    if primer:
        return primer[0], primer[1], len(seq) - primer[2] + 1
    else:
        print(len(seq), pos)
        return

def kyte_doolittle_score(aa: str) -> float:
    score = {
        "A": 1.8,
        "C": 2.5,
        "D": -3.5,
        "E": -3.5,
        "F": 2.8,
        "G": -0.4,
        "H": -3.2,
        "I": 4.5,
        "K": -3.9,
        "L": 3.8,
        "M": 1.9,
        "N": -3.5,
        "P": -1.6,
        "Q": -3.5,
        "R": -4.5,
        "S": -0.8,
        "T": -0.7,
        "V": 4.2,
        "W": -0.9,
        "Y": -1.3,
        "*": 0
    }
    return score[aa.upper()]
