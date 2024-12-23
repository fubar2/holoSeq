# see https://github.com/fubar2/holoSeq
# generic pre-computed holoseq.gz display
# any number of input files. Metadata is read to
# discover the display class to call to make a stacked panel from
# each one.

import argparse

import logging


import holoviews as hv
import panel as pn

import holoseq_data
import gff
import bigwig
import pair2d


logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_display")

# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/holoSeqtest.gz"


def xportHtml(fname, hObj):
    "save a holoview object to an interactive but not adaptive scaling HTML page"
    hv.save(filename=fname, obj=hObj)


parser = argparse.ArgumentParser(description="", epilog="")
parser.add_argument(
    "--inFile",
    help="gzipped hseq coordinates and contigs",
    default="mUroPar1_cis1.hseq.gz",
    nargs="+",
)
parser.add_argument(
    "--size", help="Display size in pixels. Default is 1000", default=1000
)
parser.add_argument("--version", "-V", action="version", version="0.1")
args = parser.parse_args()
pwidth = int(args.size)
outp = None
for i, infile in enumerate(args.inFile):
    valid, metadata = holoseq_data.getMetadata(infile)
    print("Infile %s isvalid=%s " % (infile, valid))
    if valid:
        cls = metadata["class"][0]
        if cls == "bigwig":
            p1, title = bigwig.makePanel(infile, pwidth)
        elif cls == "gff":
            p1, title = gff.makeGFFPanel(infile, pwidth)
        elif cls == "pair2d":
            p1, title = pair2d.makePafPanel(infile, pwidth)
        else:
            logging.warn(
                "infile %s metadata has unknown holoseq.gz class %s so cannot display- expect pair2, bigwig or gff"
                % (infile, cls)
            )
    if outp is None:
        outp = p1
    else:
        outp = outp + p1
pn.Row(outp).servable(title=title)
