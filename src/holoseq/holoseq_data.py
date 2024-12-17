# holoSeq compressed pre-computed plots

import array
from collections import OrderedDict
from functools import cmp_to_key
from pathlib import Path

import gzip
import itertools
import logging


import re

from holoseq.config import VALID_HSEQ_FORMATS
from holoseq.exceptions import HoloSeqFormatError


logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_prepare")


def getHap(contig, hap_indicator="None"):
    """
    vgp has 2 haplotypes and these need to be distinguished in tracks, so often have a suffix like "H1"

    Need a function to return suffix H1 from contig name "chrH1" - adjust to suit the specific conventions used by
    whoever assembles your genomes.
     help="None, Suffix (H[1,2]) Dashsuffix (_H...)"
    """
    if hap_indicator == "None":
        return "H1"
    elif hap_indicator == "Suffix":
        return contig[-2:]
    elif hap_indicator == "Dashsuffix":
        return contig.split("_")[-1]


def getContigs(lenFile):
    # samtools faidx will make one of these from a genome fasta
    # whitespace delimited contig names and lengths.
    contigs = []
    seen = {}
    haps = []
    with open(lenFile) as lf:
        for i, row in enumerate(lf):
            row = [x.strip() for x in row.strip().split()]
            if len(row) > 1:
                c, clen = row[:2]
                if seen.get(c, None):
                    log.debug("Contig %s seen again at row %d of %s" % (c, i, lenFile))
                else:
                    seen[c] = c
                    contigs.append((c, int(clen)))
                h = getHap(c)
                if h not in haps:
                    haps.append(h)
    return contigs, haps


def VGPsortfunc(s1, s2):
    """
    Specific sort for a set of contig naming conventions seen in VGP data
    # big fugly hack to sort super contigs before anything else
    # then by contig number or if they are the same offset
    # ('SUPER_2H1', 226668729) , ('SUPER_1H1', 284260672), ('SUPER13_unloc_5H1',..), (Scaffold_aa16H2, ...)
    # or chr10H1 etc
    # always work with uppercase
    # ^([^_]+)_(\S+)H[12]{1}$ gives group1 SUPER or Scaffold and group2 X or 22
    #     ^([^_]+)_([^_]+)_unloc_(\S+)H[12]{1}$ gives super aa x for super_aa_unloc_XH2
    """
    ssc = re.compile(r"^(CHR)_*(\d+|[a-zA-Z0-9]+)_*(\S*)$")
    # should match chrY or chr_y_paternal or chr333
    ss1 = re.compile(
        r"^([^_]+)_(\S+)H[12]{1}$"
    )  # should match VGP non-unloc - super/chr/scaffold. Always uppercase at entry.
    ss2 = re.compile(r"^([^_]+)_([^_]+)_UNLOC_(\S+)H[12]{1}$")  # for the unloc

    def intOrd(n):
        "x, y, z, w typically sex chromosomes - there be some weird beasts"
        if n:
            n = n.replace("_", "").replace("chr", "")  # in case chr_aa or something
            if n.isdigit():
                return int(n)
            else:
                return ord(n[0])
        else:
            return None

    def matchme(s):
        c1 = n1 = n2 = None
        found = ssc.search(s)
        if found:
            c1, n1, n2 = found.groups()[:3]
        else:
            found = ss2.search(s)
            if found:
                c1, n1, n2 = found.groups()[:3]
            else:
                found = ss1.search(s)
                if found:
                    c1, n1 = found.groups()[:2]
                    n2 = None
        if n1.isdigit():
            n1 = int(n1)
        else:
            n1 = ord(n1[0])
        if n2:
            if n2.isdigit():
                n2 = int(n2)
            else:
                n2 = ord(n2[0])
        return c1, n1, n2

    s1 = [s1[0].upper(), s1[1]]
    s2 = [s2[0].upper(), s2[1]]
    if s1[0] == s2[0]:  # simplest case - same contig, sort on offset
        return s1[1] - s2[1]  # neg if left sorts before
    u1 = "UNLOC" in s1[0]
    u2 = "UNLOC" in s2[0]
    if u1 and not u2:
        return 1  # u1 unloc goes after
    elif u2 and not u1:
        return -1  # u1 goes before unloc
    isSuper1 = (not u1) and (("SUPER" in s1[0]) or ("CHR" in s1[0]))
    isSuper2 = (not u2) and (("SUPER" in s2[0]) or ("CHR" in s2[0]))
    isScaff1 = (not u1) and (("SCAFFOLD" in s1[0]))
    isScaff2 = (not u2) and (("SCAFFOLD" in s2[0]))
    if isSuper1 and not isSuper2:
        return -1
    elif isSuper2 and not isSuper1:
        return 1
    # Must parse

    c1, naa, nab = matchme(s1[0].upper())
    c2, nba, nbb = matchme(s2[0].upper())
    if not c1 or not c2:
        log.debug("no match for ss1 %s and/or ss2 %s" % (ss1, ss2))
        return 0
    else:
        if isSuper1 or (
            isScaff1 and isScaff2
        ):  # if a super must both be supers or if both are scaffolds
            return naa - nba
        else:  # must both be unlocs
            if naa == nba:
                return nab - nbb
            else:
                return naa - nba


def Lengthsortfunc(s1, s2):
    """ """
    return s1[1] - s2[1]  # neg if left sorts before


def contsort(contigs, args):
    # sort and return offsets to starts of each contig
    # hstarts = list(itertools.accumulate(hlens))
    if args.contig_sort.lower() == "vgpname":
        contigs.sort(key=cmp_to_key(VGPsortfunc))
    elif args.contig_sort.lower() == "name":
        contigs.sort()
    elif args.contig_sort.lower() == "length":
        contigs.sort(key=cmp_to_key(Lengthsortfunc), reverse=True)
    clens = [x[1] for x in contigs]
    cnames = [x[0] for x in contigs]
    cstarts = list(itertools.accumulate(clens))
    cstarts.insert(0, 0)  # first one starts at 0
    scont = OrderedDict(zip(cnames, cstarts))
    maxpos = cstarts[-1] + clens[-1]
    return scont, maxpos


def is_valid_header(header: str) -> tuple[bool, int, str]:
    tokens = [token.strip() for token in header.split()]

    check = False
    hseq_format = tokens[0]
    if hseq_format not in VALID_HSEQ_FORMATS:
        msg = f"Not a valid holoSeq file. Header must start with one of: {VALID_HSEQ_FORMATS}."
        raise HoloSeqFormatError(msg)
    else:
        check = True

    num_dimensions = VALID_HSEQ_FORMATS.index(hseq_format) + 1
    plot_type = "bar"
    if num_dimensions == 1 and len(tokens) > 1:
        plot_type = tokens[1]

    return check, num_dimensions, plot_type


def load(path: Path):
    haploids = {}
    x_coords = array.array("l")
    y_coords = array.array("l")
    annotations = []
    metadata = {}
    gff_data = []
    rotated = False
    hh = []

    is_gff = False
    num_dimensions = -1
    plot_type = "bar"

    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            # Check if we have a valid holoSeq file using the first line as a header, and raise an
            # error if it is not.
            if i == 0:
                _, num_dimensions, plot_type = is_valid_header(line)
                continue

            if line[0] == "@":
                row = line[1:]
                if row[0] == "@":
                    tokens = [token.strip() for token in row[1:].split()]
                    metadata[tokens[0]] = tokens[1:]
                    if tokens[0] == "GFF":
                        is_gff = True
                    if tokens[0] == "rotated":
                        if tokens[1] == "True":
                            rotated = True

                else:
                    tokens = [token.strip() for token in row.split()]

                    if len(tokens) >= 3:
                        haploid_name, contig_name, position = tokens[:3]
                        if not haploids.get(haploid_name, None):
                            haploids[haploid_name] = {
                                "contig_names": [],
                                "positions": array.array("l"),
                            }
                            hh.append(haploid_name)

                        if num_dimensions == 2:
                            haploids[haploid_name]["contig_names"].append(contig_name)
                            haploids[haploid_name]["positions"].append(int(position))
                    else:
                        msg = (
                            "NOT A VALID holoSeq FILE.\n"
                            f"Line {i} of {str(path)} does not include a `reference name` "
                            "`contig name`, and `contig length`."
                        )
                        raise HoloSeqFormatError(msg)
            else:
                tokens = [token.strip() for token in line.split()]
                num_tokens = len(tokens)

                if num_dimensions == 2:
                    if num_tokens < 2:
                        msg = (
                            "NOT A VALID holoSeq FILE.\n"
                            f"Line {i} of {str(path)} needs at least two valid integer "
                            "coordinates to be a valid 2D holoSeq file."
                        )
                        raise HoloSeqFormatError(msg)

                    if num_tokens >= 2:
                        coords_check = tokens[0].isdigit() and tokens[1].isdigit()
                        if coords_check:
                            x_coords.append(int(tokens[0]))
                            y_coords.append(int(tokens[1]))
                            if num_tokens > 2:
                                annotations.append(tokens[2:])
                        else:
                            msg = (
                                "NOT A VALID holoSeq FILE.\n"
                                f"Line {i} of {str(path)} needs at least two valid integer "
                                "coordinates to be a valid 2D holoSeq file."
                            )
                            raise HoloSeqFormatError(msg)
                else:
                    if is_gff:
                        gff_data.append(tokens)

                    else:
                        if tokens[0].isdigit():
                            x_coords.append(int(tokens[0]))
                            if num_tokens > 1:
                                y_coords.append(int(tokens[1]))
                            if num_tokens > 2:
                                annotations.append(tokens[2:])
                        else:
                            msg = (
                                "NOT A VALID holoSeq FILE.\n"
                                f"Line {i} of {str(path)} needs at least one valid integer to be "
                                "a valid 1D holoSeq file."
                            )
                            raise HoloSeqFormatError("Not a valid holoSeq file.")
        if len(hh) < 2:
            log.debug("extending haps %s" % hh)
            hh.append(hh[0])
        hh.sort()
        return (num_dimensions, haploids, x_coords, y_coords, annotations, plot_type, metadata, gff_data, hh, rotated)
