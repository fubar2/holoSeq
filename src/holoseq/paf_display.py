import argparse
from pathlib import Path

import panel as pn

from holoseq import pair2d
from holoseq import holoseq_data

pn.extension()

parser = argparse.ArgumentParser(description="", epilog="")

parser.add_argument(
    "--size", help="Display size in pixels. Default is 1000", default=1000
)


parser = argparse.ArgumentParser(description="", epilog="")
parser.add_argument(
    "--inFile",
    help="PAF with paired alignments as text or gzip",
    default="../test-data/small.paf",
)
parser.add_argument(
    "--xclenfile",
    help="X axis genome contig names and lengths, whitespace delimited",
    required=True,
)
parser.add_argument(
    "--yclenfile",
    help="Y axis genome contig names and lengths, whitespace delimited.",
    required=True,
)
parser.add_argument(
    "--addH1",
    help="Bigwig and gff contigs can have H1 added if that matches the supplied xclenfile contig names. Not recommended - best to map against the right fasta",
    action="store_true",
    default=False,
)

parser.add_argument("--title", help="Title for the plot", default="Title")
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
parser.add_argument(
    "--show",
    help="Use panel to show the plot after creating the holoseq.gz output file",
    action="store_true",
    default=False,
)
parser.add_argument("--version", "-V", action="version", version="0.1")
args = parser.parse_args()
haps = []
yhaps = []
xcontigs, xhaps = holoseq_data.getContigs(args.xclenfile)
sxcontigs, xwidth = holoseq_data.contsort(xcontigs, args)
if args.yclenfile:
    ycontigs, yhaps = holoseq_data.getContigs(args.yclenfile)
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
haps.sort()
f = args.inFile
pair = pair2d(f, args, sxcontigs, sycontigs, haps, xwidth, ywidth)
ps = Path(f).suffix.lower()
log.debug("inFile=%s, ftype = %s" % (f, ps))
p = None
if args.show:
    counts = [pairs.nrcis1, nrcis2, nrtrans]
    for i, fname in enumerate([pair.cis0n, pair.cis1n, pair.transn]):
        if counts[i] > 0:
            thisp, t = pair.makePanel(fname, pwidth)
        if p == None:
            p = thisp
        else:
            p = p + thisp
    pn.serveable(p)


path = Path(args.inFile[0]).resolve()
width = int(args.size)
pairs = pair2d
holoseq_panel, title = maker(path=path, width=width)
pn.panel(holoseq_panel).servable(title=title)
