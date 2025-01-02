# for Mashmap paf, python holoSeq_prepare_gz.py --inFile  hg002_2k99.paf --title "hg002 Mashmap" --hap_indicator None --contig_sort length
# for HiC pairs
# python holoSeq_prepare_gz.py --inFile mUroPar1H1H2.paf --xclenfile mUroPar1H1suffix.len --yclenfile mUroPar1H2suffix.len --contig_sort VGPname --hap_indicator Suffix --title "VGP mUroPar1 HiC data"
# panel serve holoSeq_display.py --show --args --inFile mUroPar1H1H2.paf_cisH1_hseq.gz mUroPar1H1H2.paf_cisH2_hseq.gz mUroPar1H1H2.paf_trans_hseq.gz  --size 1000
# or
# panel serve holoSeq_display.py --show --args --inFile mUroPar1H1H2.paf_cisH1_hseq.gz --size 1000
#
# python holoSeq_prepare_gz.py --inFile mUroPar1_protH1.gff --xclenfile mUroPar1H1suffix.len --contig_sort VGPname --title "mUroPar1 NCBI protein GFF"
#
# python holoSeq_prepare_gz.py --inFile ../hg002_bothHiC.paf --xclenfile hg002H1_suffixed.len --yclenfile hg002H2_suffixed.len --contig_sort VGPname --hap_indicator Suffix --title "T2T HG002 HiC data"
#
# Proof of concept data are Arima HiC reads from the Arctic Ground Squirrel mUroPar1 VGP genomeArk repository
# processed with Dephine's Pretext workflow using Bellerophon to remove chimeric reads
# The paired bam output is converted to PAF with an awk script (!) and that's what is read
# in this code.
# The pairs are parsed to extract the haplotype designator from the contig name
# typically a suffix like H1 extracted in getHap - rewrite that to suit your names.
# Contig ordering really matters for the plots to make any sense.
# Ideally, curators name them so they sort alphanumerically without effort.
# There's a sorthapqname function that is used here. Designed for the VGP data
# Will need to be replaced for other naming conventions.
# It's a mess.
# Sorting by contig name is based on VGP conventions - SUPER_ first, then scaffolds
# One problem to watch out for is that any differences in ordering of the X and Y contigs can make all sorts of
# artifacts appear such as the kaleidoscopic patterns seen in Pretextviewer.

# Ross Lazarus October 2024

import argparse

import logging


from holoseq import holoseq_data
from holoseq import gff
from holoseq import bigwig
from holoseq import pair2d


logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_prepare")

# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/huge.paf"



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument(
        "--inFile",
        help="PAF with paired alignments, bigwig or gff3",
       required=True,
    )
    parser.add_argument(
        "--inFtype",
        help="Only pair2d, bigwig or gff are currently supported",
       required=True,
    )
    parser.add_argument(
        "--xclenfile",
        help="X axis contig names and lengths, whitespace delimited",
        required=True,
    )
    parser.add_argument(
        "--addH1",
        help="Bigwig and gff contigs can have H1 added if that matches the supplied xclenfile contig names. Not recommended - best to map against the right fasta",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--yclenfile",
        help="Optional Y axis contig names and lengths, whitespace delimited for different reference sequences. Required for 2D plots",
        required=False,
    )
    
    parser.add_argument(
        "--title", help="Title for the plot", default="Plot title goes here"
    )
    parser.add_argument(
        "--contig_sort", help="VGPname, name, length, none", default="length"
    )
    parser.add_argument(
        "--refURI",
        help="URI for the genome reference sequence used for the coordinates for metadata",
        default="Unknown",
    )
    parser.add_argument(
        "--hap_indicator",
        help="None, Suffix (H[1,2]) Dashsuffix (_H...)",
        default="None",
    )  
    parser.add_argument(
        "--rotate",
        help="Rotate the 2D plot so the diagonal becomes the x axis",
        action="store_true",
        default=False,
    )
    parser.add_argument("--version", "-V", action="version", version="0.1")
    args = parser.parse_args()
    haps = []
    yhaps = []
    xcontigs, xhaps = holoseq_data.getContigs(args.xclenfile, args.hap_indicator)
    sxcontigs, xwidth = holoseq_data.contsort(xcontigs, args)
    if args.yclenfile:
        ycontigs, yhaps = holoseq_data.getContigs(args.yclenfile, args.hap_indicator)
        sycontigs, ywidth = holoseq_data.contsort(ycontigs, args)
    else:
        sycontigs = sxcontigs
        ywidth = xwidth
    for h in xhaps + yhaps:
        if h not in haps:
            haps.append(h)
    if len(haps) == 1:
        log.debug("extending haps %s" % haps)
        haps.append(haps[0])
    log.debug('***haps %s' % haps)
    ps = args.inFtype.lower()
    log.debug("inFile=%s, ftype = %s" % (args.inFile, ps))

    if ps == "pair2d":
        p = pair2d.pair2d(args.inFile, args, sxcontigs, sycontigs, haps, xwidth, ywidth)
        outs = p.convert()
    elif ps in ["bw", "bigwig"]:
        outf = "%s.hseq.gz" % args.inFile
        p = bigwig.bigwig(args.inFile, outf, args, sxcontigs)
        p.convert()
    elif ps in ["gff3", "gff"]:
        outf = "%s.hseq.gz" % args.inFile
        p = gff.gff(args.inFile, outf, sxcontigs, args)
        log.debug("contigs=%s" % contigs)
        p.convert()
    else:
        log.warn("%s unknown type - cannot process" % ps)
    logging.shutdown()
