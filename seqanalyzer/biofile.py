# -*- coding: utf-8 -*-
"""Module for Biology File operations."""

from .sequence_lib import rc_seq
from .sequence_lib import translate_exon
import re
import time
import os
import numpy as np
import matplotlib.pyplot as plt



def _message(mess=""):
    """Print a message with time.
    Args:
        mess: message to be printed

    """
    time_str = time.strftime("%Y-%m-%d %H:%M:%S")
    print("%s: %s" % (time_str, mess))


def simple_gff3_load(fname, return_fasta=False):
    """Load gff3 files."""
    entries = []
    with open(fname, "r") as gff3:
        for line in gff3:
            line = line.strip()
            if line[0] == "#":
                if line == "##FASTA":
                    break
                continue
            entry = line.split("\t")
            if len(entry) < 9:
                _message("Error loading %s: Less than 9 items in the entry\n%s" % (fname, line))
                continue
            notes = entry.pop().split(";")
            entry.extend(notes)
            entry[3] = int(entry[3])
            entry[4] = int(entry[4])
            entries.append(entry)
        _message("Gff entries are analyzed")

    if not return_fasta:
        return entries
    else:
        names = []
        seqs = []
        with open(fname, "r") as gff3:
            line = gff3.readline()
            line = line.strip()
            while not line == "##FASTA":
                line = gff3.readline()
                line = line.strip()
            _cur_name = ""
            _cur_seq = []
            line = gff3.readline()
            while line:
                line_ = line.strip()
                line = gff3.readline()
                # Skip empty lines
                if not line_:
                    continue
                # New entry
                if line_[0] == ">":
                    # Add previous entry
                    if _cur_name:
                        names.append(_cur_name)
                        seqs.append("".join(_cur_seq))
                    _cur_name = line_[1:] # Omit >
                    _cur_seq = []
                else:
                # Update sequence of the entry
                    if not _cur_name:
                        _message("One seq without entry")
                        _message(line_)
                        return entries, [], []
                    _cur_seq.append(line_)
            # Update the final entry
            if _cur_name:
                names.append(_cur_name)
                seqs.append("".join(_cur_seq))
        return entries, names, seqs


def simple_fasta_write(fname, names, seqs, linewidth=80):
    """Write the fasta file."""
    with open(fname, "w") as fasta_p:
        total = len(names)
        for _id in range(total):
            fasta_p.write(">%s\n" %names[_id])
            _seq = seqs[_id]
            # seq with linebreak
            _seq_wb = "\n".join([_seq[i:linewidth+i]
                                 for i in range(0, len(_seq), linewidth)])
            fasta_p.write(_seq_wb)
            fasta_p.write("\n")


def simple_fastq_load(fname: str) -> tuple:
    """Load Fastq file with Phred Score in 32-ASCII code."""
    names = []
    seqs = []
    scores = []
    if not os.path.isfile(fname):
        _message("%s is not available!" % fname)
        return [], [], []
    with open(fname) as filep:
        _cur_name = ""
        _cur_seq = []
        _cur_score = []
        for count, line in enumerate(filep):
            line = line.strip()
            if not line:
                continue
            # New entry
            if count %4 == 0:
                if line[0] != '@':
                    _message("Fastq file reading error in reading %s" %line[:20])
                    return [], [], []
                if _cur_name:
                    names.append(_cur_name)
                    seqs.append(_cur_seq)
                    scores.append(_cur_score)
                _cur_name = line[1:]
                _cur_seq = []
                _cur_score = []
            elif count %4 == 1:
                _cur_seq = line
            # load Phred Score
            elif count %4 == 3:
                _cur_score = [ord(phred) - 33 for phred in line]
                if len(_cur_score) != len(_cur_seq):
                    _message("Length of sequence and score doesnt match for %s" %_cur_name)
                    return [], [], []
        names.append(_cur_name)
        seqs.append(_cur_seq)
        scores.append(_cur_score)
    return names, seqs, scores


def simple_fasta_load(fname: str) -> tuple:
    """Load fasta file.
    Args:
    fname: name of the fasta file
    Returns:
    names, seqs: name list and seqlist
    """

    names = []
    seqs = []
    if not os.path.isfile(fname):
        _message("%s is not available!" % fname)
        return [], []

    with open(fname, 'r') as fasta_p:
        _cur_name = ""
        _cur_seq = []
        for line in fasta_p:
            line_ = line.strip()
            # Skip empty lines
            if not line_:
                continue
            # New entry
            if line_[0] == ">":
                # Add previous entry
                if _cur_name:
                    names.append(_cur_name)
                    seqs.append("".join(_cur_seq))
                _cur_name = line_[1:]  # Omit >
                _cur_seq = []
            else:
                # Update sequence of the entry
                if not _cur_name:
                    _message("One seq without entry")
                    _message(line_)
                    return [], []
                _cur_seq.append(line_)

        # Update the final entry
        if _cur_name:
            names.append(_cur_name)
            seqs.append("".join(_cur_seq))
    return names, seqs


def load_fasta(fname: str) -> tuple:
    """Load fasta file.
    Args:
    fname: name of the fasta file
    Returns:
    genome = {names: seqs}
    """
    genome = {}

    with open(fname, 'r') as fasta_p:
        _cur_name = ""
        _cur_seq = []
        for line in fasta_p:
            line_ = line.strip()
            # Skip empty lines
            if not line_:
                continue
            # New entry
            if line_[0] == ">":
                # Add previous entry
                if _cur_name:
                    genome[_cur_name] = "".join(_cur_seq)
                _cur_name = line_[1:]  # Omit >
                _cur_seq = []
            else:
                # Update sequence of the entry
                if not _cur_name:
                    _message("One seq without entry")
                    _message(line_)
                    return [], []
                _cur_seq.append(line_)

        # Update the final entry
        if _cur_name:
            genome[_cur_name] = "".join(_cur_seq)
    return genome


class blast_7_entry(object):
    """Class to store the blast tab entry.
    The blast tab file generated with the following option
    is analysed.
    -outfmt "7 qseqid sseqid qstart qend sstart send \
     length  qlen slen evalue bitscore nident pident positive"
    """
    def __init__(self, line:str):
        """Inits the entry class."""
        entry = line.split()
        self.query = entry[0]
        self.ref = entry[1]
        self.qstart = int(entry[2])
        self.qend = int(entry[3])
        self.sstart = int(entry[4])
        self.send = int(entry[5])
        self.len_align = int(entry[6])
        self.len_query = int(entry[7])
        self.len_sub = int(entry[8])
        self.e = float(entry[9])
        self.bit = float(entry[10])
        self.ide = int(entry[11])
        self.p_ide = float(entry[12])
        self.pos = int(entry[13])
        if self.sstart < self.send:
            self.direction = True
        else:
            self.direction = False
        self.p_pos = self.pos / self.len_query * 100
        self.line = line.strip()
        self.line += "\t%2.2f\n" %self.p_pos


def extact_blast_tab(dirn:str, blast:str, min_p_pos=80):
    """Extract the information from blast tab file."""
    os.chdir(dirn)
    filter_entries = []
    with open(blast, "r") as filep:
        for line in filep:
            if line[0] == "#":
                continue
            entry = blast_7_entry(line)
            if entry.p_pos >= min_p_pos:
                filter_entries.append(entry)
    return filter_entries


def extract_extend_cDNA(start:int, end:int, direction:bool , ref_seq: str) -> tuple:
    """Extract the cDNA region and etend to full orf."""
    cstart = "ATG"
    cend = ["TAG", "TGA", "TAA"]
    crstart = "CAT"
    crend = ["CTA", "TCA", "TTA"]
    if direction:
        seq = ref_seq[start - 1: end]
    else:
        seq = ref_seq[end -1: start]

    if direction:
        trans = translate_exon(seq)
    else:
        trans = translate_exon(rc_seq(seq))
    if "*" in trans[:-1]:
        is_inner_stop = True
    else:
        is_inner_stop = False
    # Extend the cDNA sequence
    if direction:
        # + alignment
        left = ""
        if seq[: 3].upper() != cstart:
            for pos in range(start, 2, -3):
                code = ref_seq[pos -3: pos]
                u_code = code.upper()
                if u_code != cstart and u_code not in cend:
                    left = code + left
                elif u_code == cstart:
                    left = code + left
                    break
                else:
                    break
        seq = left + seq
        left_flank = len(left)
        right = ""
        if seq[-3: ].upper() not in cend and not is_inner_stop:
            for pos in range(end +1, len(ref_seq) -3, 3):
                code = ref_seq[pos: pos + 3]
                u_code = code.upper()
                if u_code not in cend:
                    right += code
                else:
                    right += code
                    break
        if is_inner_stop:
            stop_pos = trans.index("*")
            seq = seq[:stop_pos * 3]
        seq = seq + right
        right_flank = len(right)
    else:
        # - alignment
        left = ""
        if seq[: 3].upper() not in crend and not is_inner_stop:
            for pos in range(end, 2, -3):
                code = ref_seq[pos -3: pos]
                u_code = code.upper()
                if u_code not in crend:
                    left = code + left
                else:
                    left = code + left
                    break
        seq = left + seq
        if is_inner_stop:
            stop_pos = trans.index("*")
            seq = seq[-stop_pos*3: ]
        left_flank = len(left)
        right = ""
        if seq[-3: ].upper() != crstart:
            for pos in range(start +1, len(ref_seq) -3, 3):
                code = ref_seq[pos: pos + 3]
                u_code = code.upper()
                if u_code not in crend and u_code != crstart:
                    right += code
                elif u_code == crstart:
                    right += code
                    break
                else:
                    break
        seq = seq + right
        right_flank = len(right)
        seq = rc_seq(seq)

    if seq[:3].upper() == cstart:
        is_full = True
    else:
        is_full = False
    return (seq, left_flank, right_flank, is_full, is_inner_stop)


def extract_tblast_n_alignment(dirn:str, blast:str, prefix:str, fasta:str, min_p_pos=80):
    """Extract the tblastn alignment region.
    The aligned region is extended to have the entire cDNA region.
    The cDNA of the region is also extracted.
    """
    align_entries = extact_blast_tab(dirn, blast, min_p_pos)
    os.chdir(dirn)
    chrom, seqs = simple_fasta_load(fasta)
    chrom = [name.split()[0] for name in chrom]
    # Extract alignment region
    align_info = []
    align_seqs = {}
    full_aligns = []
    full_prots = []
    for entry in align_entries:
        query = entry.query
        sstart = entry.sstart
        send = entry.send
        ref = entry.ref
        ref_seq = seqs[chrom.index(ref)]
        seq, left_flank, right_flank, is_full, is_innerstop = extract_extend_cDNA(
            sstart, send, entry.direction, ref_seq)
        if is_innerstop:
            div_mod = abs(send - sstart) +1
            entry.pos = entry.pos * len(seq) // div_mod
            entry.len_align = entry.len_align * len(seq) // div_mod
            entry.ide = entry.ide * len(seq) // div_mod
            if entry.direction:
                # + alignment
                send = sstart + len(seq) -1
            else:
                send = sstart - len(seq) + 1
            entry.qend = entry.qstart + len(seq)//3 -1
        align_seqs.setdefault(query, [])
        id_ = len(align_seqs[query]) + 1
        align_seqs[query].append(("%s%03d" %(query, id_), seq))
        line = entry.line
        line = line.strip()
        _tmp_line = line.split()
        _tmp_line[3] = str(entry.qend)
        _tmp_line[4] = str(sstart)
        _tmp_line[5] = str(send)
        _tmp_line[6] = str(entry.len_align)
        _tmp_line[11] = str(entry.ide)
        _tmp_line[13] = str(entry.pos)
        update_pos = "%2.2f" % (entry.pos / entry.len_query * 100)
        _tmp_line[14] = update_pos
        line = "\t".join(_tmp_line)
        p_pos = entry.pos / (entry.len_align + left_flank/3 + right_flank/3) * 100
        line += "\t%d\t%d\t%d\t%2.2f\t%s\t%s\n" %(id_, left_flank/3, right_flank/3, p_pos, is_full, is_innerstop)
        align_info.append(line)
        if is_full:
            full_aligns.append(("%s %03d %d-%d" %(query, id_, sstart - left_flank, send + right_flank), seq))
            full_prots.append(("%s %03d %d-%d" %(query, id_, sstart - left_flank, send + right_flank), translate_exon(seq)))
    aligns = []
    for align in align_seqs.values():
        aligns.extend(align)
    align_names, align_seqs = zip(*aligns)
    full_names, full_seqs = zip(*full_aligns)
    prot_names, prot_seqs = zip(*full_prots)
    simple_fasta_write(prefix + ".aligns.fa", align_names, align_seqs)
    simple_fasta_write(prefix + ".aligns.full.fa", full_names, full_seqs)
    simple_fasta_write(prefix + ".aligns.full.prot.fa", prot_names, prot_seqs)
    with open(prefix + ".aligns.tab", "w") as filep:
        lines = ["query\tref\tqstart\tqend\tsstart\tsend\tlen_align\tlen_q\tlen_r\tE\tbit\tide\t%%ide\tpos\t%%pos\talign_no\tleft_flank\tright_flank\t%%pos\tis_full\tis_innerstop\n"]
        lines.extend(align_info)
        filep.write("".join(lines))
    # Filter multi full aligbments
    multi_align = {}
    for name in prot_names:
        gene = name.split()[0]
        multi_align.setdefault(gene, [])
        multi_align[gene].append(name)
    multi_names = []
    multi_seqs = []
    for gene, name in multi_align.items():
        if len(name) == 1:
            continue
        multi_names.extend(name)
        multi_seqs.extend([prot_seqs[prot_names.index(name_)] for name_ in name])
    simple_fasta_write(prefix + ".aligns.multi.prot.fa", multi_names, multi_seqs)

# class BlastAlignInfo(object):
#     """Class to store the summary of the blast result."""
#     def __init__(self, entry: list):
#         """Inits the BlastAligh object.
#         Args:
#         entry:
#             [length_query, length_ref, score, expect, ide, pos, length_align, gaps, frame, id]
#         """
#         self.score = entry[2]
#         self.expect = entry[3]
#         self.ide = entry[4]
#         self.pos = entry[5]
#         self.length_align = entry[6]
#         self.length_ref = entry[1]
#         self.gap = entry[7]
#         self.frame = entry[8]
#         self.lengt_query = entry[0]
#         self.id = entry[9]

# class BlastAlignCoverage(object):
#     """Coverage information for blast result.
#         Attributes:
#         coverage (list): per a.a alignment
#         coverage[i] = [align_id, align_status, frame]
#         align_status: "i" = Ide, "p" = Pos, "d" = Deletion, "a" = Insertion at the (i+1) th
#     """
#     def __init__(self, alignlines: list, length:int):
#         """Inits the BlastAlignCoverage object.
#             Args:
#             alignlines: lines in the blast result
#             length: length of the query
#         """
#         pass



# class tBlastnEntry(object):
#     """Class to store the tblastn Entry.

#     Attributs:
#     database (str): Name of the database
#     query (str): Name of the query
#     align (dict): Summary of aligns
#         align {(query, ref): {align_id: BlastAlignInfo}}
#     align_id (dict) : {align_id: (query, ref)}
#     coverage (list): per a.a alignment
#         coverage[i] = [align_id, align_status, frame]
#         align_status: "i" = Ide, "p" = Pos, "d" = Deletion, "a" = Insertion at the (i+1) th a.a


#     """
#     def __init__(self, entrylines=None):
#         """Inits the tBlastnEntry Class.
#         Args:
#             entrylines: lines of the entry of the alignment
#         """
#         self.database = ""
#         self.query = ""
#         self.ref = []
#         self.align = {}
#         self.align_id = {}
#         self.coverage = []
#         self.query_len = 0
#         if entrylines:
#             self.load_entrylines(entrylines)

#     def load_entrylines(self, entrylines):
#         """Load the entryline information."""
#         length = len(entrylines)
#         n = 0
#         # Load header information
#         while n < length:
#             line = entrylines[n]
#             if line.startswith("Database:"):
#                 entry = line.split()
#                 self.database = entry[1]
#                 n += 1
#             elif line.startswith("Query= "):
#                 self.query = line.split("=")[1]
#                 n += 1
#             elif line.startswith("Length= "):
#                 self.query_len = int(line.split("=")[1])
#             if line.startswith(">"):
#                 n += 1
#                 break
#         id_ = 0
#         # Load Alignment info
#         while n < length:
#             line = entrylines[n]
#             entry = line.split(",")
#             # new align
#             if entry[0].startswith("Score"):
#                 id_ += 1






class TmpTBNAlign(object):
    """Tmp"""
    def __init__(self, info):
            self.align_id = info[0]
            self.sub_s = info[1]
            self.sub_e = info[2]
            self.query_s = info[3]
            self.query_e = info[4]
            self.inter_stop = info[5]
            self.q_del = info[6]
            self.q_ins = info[7]
            self.q_omet = info[8]
            self.align_frame = info[9]
            self.align_aa = info[10]
    def info(self):
        return [self.align_id, self.sub_s, self.sub_e, self.query_s, self.query_e,
                self.inter_stop, self.q_del, self.q_ins, self.q_omet, self.align_frame]


def tmp_tblastn_analyzer(dirn, blast, rep=3, prefix="out", gap=100, align_len=60, dpi=300):
    """fast analyze tblastn result."""
    os.chdir(dirn)
    tblastn_entries = []
    tblastn_prot = []
    tblastn_summ = []
    tblastn_align = {}
    with open(blast, "r") as filep:
        current = []
        cur_query = ""
        q_len = 0
        ongoing = False
        onalign = False
        align_id = 0
        align_frame = 0
        align = []
        cur_s = 0
        cur_e = 0
        for line in filep:
            # New query
            if line.startswith("Query="):
                ongoing = True
                onalign = False
                align_id = 0
                cur_query = line.split("=")[1]
                cur_query = cur_query.strip()
                align = []
                continue
            if line.startswith("Lambda"):
                ongoing = False
                onalign = False
                if align_id != 0:
                    if align_id == 1:
                        cur_s = min([sub_s, sub_e])
                        cur_e = max([sub_s, sub_e])
                        align.append(TmpTBNAlign((align_id, sub_s, sub_e, query_s, query_e,
                                     inter_stop, q_del, q_ins, q_omet, align_frame, align_aa)))
                    else:
                        if abs(cur_s - sub_s) < gap or abs(cur_e - sub_e) < gap or abs(cur_s - sub_e) < gap or abs(cur_e - sub_s) < gap:
                            align.append(TmpTBNAlign((align_id, sub_s, sub_e, query_s, query_e,
                                         inter_stop, q_del, q_ins, q_omet, align_frame, align_aa)))
                tblastn_align[cur_query] = (q_len, align)
                continue
            if not ongoing:
                continue
            if line.startswith("Length=") and not onalign:
                q_len = int(line.split("=")[1])
                align_aa = [0] * q_len
                continue
            if line.startswith(">"):
                onalign = True
                onseq = False
                continue
            if line.startswith(" Score"):
                onseq = False
                if align_id != 0:
                    if align_id == 1:
                        cur_s = min([sub_s, sub_e])
                        cur_e = max([sub_s, sub_e])
                        align.append(TmpTBNAlign((align_id, sub_s, sub_e, query_s, query_e,
                                     inter_stop, q_del, q_ins, q_omet, align_frame, align_aa)))
                    else:
                        if abs(cur_s - sub_s) < gap or abs(cur_e - sub_e) < gap or abs(cur_s - sub_e) < gap or abs(cur_e - sub_s) < gap:
                            align.append(TmpTBNAlign((align_id, sub_s, sub_e, query_s, query_e,
                                         inter_stop, q_del, q_ins, q_omet, align_frame, align_aa)))
                            # Update region
                            cur_s = min([cur_s, sub_s, sub_e])
                            cur_e = max([cur_e, sub_s, sub_e])
                align_id += 1
                sub_s = 0
                sub_e = 0
                query_s = 0
                query_e = 0
                inter_stop = 0
                q_del = 0
                q_ins = 0
                is_end = False
                q_omet = False
                ins_id = []
                len_header = 0
                continue
            if line.startswith(" Frame"):
                align_frame = int(line.split("=")[1])
                continue
            if not onalign:
                continue
            if line.startswith("Query"):
                onseq = True
                if query_s == 0:
                    query_s = int(line.split()[1])
                    if line.split()[2][0].upper() == "M":
                        q_omet = True
                query_e = int(line.split()[-1])
                q_del += line.count("-")
                header_pat = r"Query\s+\d+\s+"
                header = re.match(header_pat, line)
                header = header.group(0)
                len_header = len(header)
                ins_id = [id_ for id_, aac in enumerate(line.split()[2]) if aac =="-"]
                if "*" in line:
                    is_end = True
                continue
            elif line.startswith("Sbjct"):
                if sub_s == 0:
                    sub_s = int(line.split()[1])
                    if q_omet:
                        if line.split()[2][0].upper() != "M":
                            q_omet = False
                sub_e = int(line.split()[-1])
                q_ins += line.count("-")
                if "*" in line and not is_end and inter_stop == 0:
                    align_ = line.split()[2]
                    inter_stop = align_.index("*") + query_e - len(align_) + 1
                elif "*" in line and is_end and inter_stop == 0:
                    align_ = line.split()[2]
                    if align_.index("*") != len(align_) -1:
                        inter_stop = align_.index("*") + query_e - len(align_) + 1
                continue
            elif onseq:
                line = line[len_header: ]
                line = line.strip()
                for id_, ami in enumerate(line):
                    if query_e > align_len:
                        pos = query_e - align_len + id_
                        if sub_e > sub_s:
                            pos_s = sub_s - align_len + id_
                        else:
                            pos_s = sub_s + align_len - id_
                    else:
                        pos = query_s + id_ -1
                        if sub_e > sub_s:
                            pos_s = sub_s + id_ -1
                        else:
                            pos_s = sub_s - id_ + 1
                    if id_ in ins_id:
                        pos -= 1
                    if ami == " ":
                        align_aa[pos] = 0
                    elif ami == "+":
                        align_aa[pos] = pos_s
                    elif ami == "*":
                        align_aa[pos] = -pos_s
                    else:
                        try:
                            align_aa[pos] = pos_s
                        except:
                            print(pos, q_len, line, cur_query)
    # Summ of alignments
    for query, aligns in tblastn_align.items():
        is_start = False
        inter_stop = []
        frame = []
        is_shift = False
        coverage = [0] * aligns[0]

        is_rep = False
        is_inter = False
        top_cover = 0
        pos_aa = [0] * aligns[0]
        for align in aligns[1]:
            if align.q_omet:
                is_start = True
            if align.inter_stop != 0:
                inter_stop.append(align.inter_stop)
            frame.append((align.query_s, align.align_frame))
            for ind in range(align.query_s -1, align.query_e):
                coverage[ind] += 1
                if align.align_aa[ind] != 0:
                    pos_aa[ind] += 1
        _pos = 0
        _count = 0
        for id_ in range(aligns[0]):
            if coverage[id_] != 0:
                pos_aa[id_] = pos_aa[id_] / coverage[id_]
                _pos += pos_aa[id_]
                _count += 1
        if _count > 0:
            pos_pr = _pos / _count *100
            pos_h = len(list(filter(lambda x: x>0, pos_aa))) / _count * 100
        else:
            pos_pr = 0
            pos_h = 0

        if max(coverage) > rep:
            is_rep = True
        top_cover = max(coverage)
        cover = (aligns[0] - coverage.count(0)) / aligns[0] * 100
        frame.sort()
        frame = [fra[1] for fra in frame]
        if len(set(frame)) > 1:
            is_shift = True
        if inter_stop:
            is_inter = True
        if pos_aa[-1] > 0:
            is_stop = True
        else:
            is_stop = False
        inter_stop.sort()
        # Summ of the ailgnment
        tblastn_summ.append((query, aligns[0], cover, is_start, is_shift, is_inter, is_rep, frame, inter_stop, top_cover, pos_pr, pos_h, is_stop))
    # Save result
    with open(prefix + ".aligns.info.tab", "w") as filep:
        lines = ["query\tlength\talign_id\tsub_s\tsub_e\tquery_s\tquery_e\tinter_stop\tq_del\tq_ins\tq_omet\talign_frame\n"]
        for query, aligns in tblastn_align.items():
            line_head = "%s\t%d\t" %(query, aligns[0])
            for align in aligns[1]:
                line = line_head + "%02d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\n" %tuple(align.info())
                lines.append(line)
        filep.write("".join(lines))
    with open(prefix + ".align.summ.tab", "w") as filep:
        lines = ["Query\tLength\tCoverage\tis_start\tis_shift\tis_inter\tis_rep\tframe\tinter_stop\ttop_coverage\tav_pos\ttop_pos\tis_stop\n"]
        for entry in tblastn_summ:
            _frame = ""
            for fra in entry[7]:
                _frame += "%s, " %str(fra)
            _frame = _frame[: -2]
            _inter = ""
            for inter in entry[8]:
                _inter += "%s, " %str(inter)
            if _inter:
                _inter = _inter[:-2]
            else:
                _inter = "NA"
            line = "%s\t%d\t%2.2f\t%s\t%s\t%s\t%s\t" %entry[:7]
            line += "%s\t%s\t%s\t%2.2f\t%2.2f\t%s\n" %(_frame, _inter, entry[9], entry[10], entry[11], entry[12])
            lines.append(line)
        filep.write("".join(lines))
    # Save protein alignment result
    # if not os.path.isdir("./" + prefix + ".protfig/"):
    #     os.mkdir("./" + prefix + ".protfig/")
    # os.chdir("./" + prefix + ".protfig/")
    # for query, aligns in tblastn_align.items():
    #     q_len = aligns[0]
    #     fname = query + ".png"
    #     fig, axes = plt.subplots(dpi=dpi)
    #     for align in aligns[1]:
    #         frame = align.align_frame
    #         print(frame, query, align.align_id)
    #         align_aa = align.align_aa
    #         align_pos = []
    #         align_ide = []
    #         align_stop = []
    #         for id_, amino in enumerate(align_aa, 1):
    #             if amino == 1:
    #                 align_pos.append(id_)
    #             elif amino == 2:
    #                 align_ide.append(id_)
    #             elif amino == -1:
    #                 align_stop.append(id_)
    #         axes.scatter(align_pos, (frame,) * len(align_pos), marker=",", color="g")
    #         axes.scatter(align_ide, (frame,) * len(align_ide), marker=",", color="b")
    #         axes.scatter(align_stop, (frame,) * len(align_stop), marker="p", color="r")
    #     axes.set_title("Coverage Map for %s" %query)
    #     axes.set_xlim([0, q_len])
    #     axes.set_ylim([-3, 3])
    #     plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    #     plt.savefig(fname, bbox_inches="tight", format="png")
    #     plt.close()



def main():
    # extract_tblast_n_alignment(
    #     dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2.canu.prot/anno/manual/",
    #     blast="bg2.man.pairwise.tab",
    #     prefix="bg2.man.2",
    #     fasta="bg2.man.interseq.fa",
    #     min_p_pos=20)
    tmp_tblastn_analyzer(
        dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2.canu.prot/anno/manual/",
        blast="bg2.man.pairwise.align",
        rep=3,
        prefix="bg2.man.3",
        gap=10000)

if __name__ == '__main__':
    main()
