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

"""oh dear
@v1HoloSeq2D
@@class pair2d
@@title VGP mUroPar1 HiC data Pairs on H1
@@datasource paf
@@datafile ./tests/mUroPar1H1H2.paf.gz
@@refURI Unknown
@@xclenfile ./tests/mUroPar1H1H2.len
@@yclenfile ./tests/mUroPar1H1H2.len
@@axes H1
@@rotated False
@H1 SUPER_1H1 0
@H2 SUPER_1H2 284260672
@H1 SUPER_2H1 568308509
@H2 SUPER_2H2 794977238

not th h2 super1 should also start wih offset 0
"""
import argparse
import copy
import logging


from holoseq import holoseq_data
from holoseq import gff
from holoseq import bigwig
from holoseq import pair2d

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

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
        help="X axis contig names and lengths, whitespace delimited - samtools faidx can generate these from the genome/haploype fasta",
        required=True,  )
    parser.add_argument(
        "--xaxis haplotype ID , - such as H1",
        required=False,
    )

    parser.add_argument(
        "--addH1",
        help="Bigwig and gff contigs can have H1 added if that matches the supplied xclenfile contig names. Not recommended - best to map against the right fasta",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--yclenfile",
        help="Optional Y axis contig names and lengths, whitespace delimited for different reference sequences. Required for 2D plots. samtools faidx can generate these from the genome/haploype fasta",
        required=False,
    ) 
    parser.add_argument(
        "--yaxishapid , - such as H2", help=' only matching contigs will be on this axis',required=False,
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
        help="None, Suffix (H[1,2]) Dashsuffix (-H...)",
        default="None",
    )  
    parser.add_argument(
        "--rotate",
        help="Rotate the 2D plot so the diagonal becomes the x axis",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--gzoutpath",
        help="Non-PAF inputs only: Path to write the gzipped precomputed plot. Default [inFile]_[inFtype]_hseq.gz",
        required=False,
    ) 
    parser.add_argument(
        "--cis1outpath",
        help="PAF inputs only: Path to write the cis hap1 gzipped precomputed plot. Default [inFile]_cis[hap]_hseq.gz",
        required=False,
    ) 
    parser.add_argument(
        "--cis2outpath",
        help="PAF inputs only: Path to write the cis hap2 gzipped precomputed plot. Default [inFile]_cis[hap]_hseq.gz",
        required=False,
    ) 
    parser.add_argument(
        "--transoutpath",
        help="PAF inputs only: Path to write the trans (both haplotypes) gzipped precomputed plot. Default [inFile]_trans_hseq.gz",
        required=False,
    ) 

    parser.add_argument("--version", "-V", action="version", version="0.1")
    args = parser.parse_args()
    haps = []
    yhaps = []
    xcontigs, xhaps = holoseq_data.getContigs(args.xclenfile, args.hap_indicator)
    sxcontigs = holoseq_data.contsort(xcontigs, args)
    if args.yclenfile:
        ycontigs, yhaps = holoseq_data.getContigs(args.yclenfile, args.hap_indicator)
        sycontigs = holoseq_data.contsort(ycontigs, args)
    else:
        sycontigs = sxcontigs
    for h in xhaps:
        if h not in haps:
            haps.append(h)
    for h in yhaps:
        if h not in haps:
            haps.append(h)
    if len(haps) == 1:
        log.debug("extending haps %s" % haps)
        haps.append(haps[0])
    log.debug('***haps %s' % haps)
    ps = args.inFtype.lower()
    log.debug("inFile=%s, ftype = %s" % (args.inFile, ps))

    if ps == "pair2d":
        p = pair2d.pair2d()
        p.inFname = args.inFile
        outs = p.convert(args)        
    elif ps in ["bw", "bigwig"]:
        if args.gzoutpath:
            outf = args.gzoutpath
        else:
            outf = "%s.hseq.gz" % args.inFile
        p = bigwig.bigwig(args.inFile)
        p.convert(outf, sxcontigs)
    elif ps in ["gff3", "gff"]:
        if args.gzoutpath:
            outf = args.gzoutpath
        else:
            outf = "%s.hseq.gz" % args.inFile
        p = gff.gff(args)
        p.convert(args.inFile, outf, sxcontigs)
    else:
        log.warn("%s unknown type - cannot process" % ps)
    logging.shutdown()
