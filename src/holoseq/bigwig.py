# bigwig tracks in holoSeq
# has functions to read a bigwig and write a hseq.gz, and to prepare a panel showing that hseq.gz 
# using the generic data.load function

from bisect import bisect_left
import gzip
import logging
import numpy as np
import os
from pathlib import Path



import pybigtools
import holoviews as hv
import pandas as pd
import panel as pn



from config import VALID_HSEQ_FORMATS
from holoseq import holoseq_data


from holoviews.operation.datashader import (
    dynspread,
)
#from holoviews.operation.element import apply_when
from holoviews.operation.resample import ResampleOperation2D
from holoviews.operation import decimate


hv.extension("bokeh", "matplotlib", width=100)

# Default values suitable for this notebook
decimate.max_samples = 1000
dynspread.max_px = 8
dynspread.threshold = 0.75
ResampleOperation2D.width = 250
ResampleOperation2D.height = 250


logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("bigwig")



logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_prepare")


class bigwig:
    """
    points = [(50*i, 100+random.random()) for i in range(10000)]
    hv.Curve(points).opts(interpolation='steps-post').opts(width=1000)
    bigwig to barchart - will autoscale to line?
    """

    def __init__(self, inFname, outFname, args, contigs):
        self.inFname = inFname
        self.hsId = VALID_HSEQ_FORMATS[0]
        fakepath = "in.bw"
        if os.path.isfile(fakepath):
            os.remove(fakepath)
        p = Path(fakepath)
        p.symlink_to(inFname)  # required by pybigtools (!)
        bwf = pybigtools.open(fakepath)
        bchrlist = bwf.chroms()
        bwchrs = list(bchrlist.keys())
        bwdata = {}
        for i, bchr in enumerate(bwchrs):
            cchr = bchr
            if (not contigs.get(bchr, None)) and args.addH1:
                cchr = cchr + "H1"
            if contigs.get(cchr, None):
                cstart = contigs[cchr]
                bwdata[cchr] = {}
                bw = bwf.records(bchr)
                # Return the records of a given range on a chromosome. The result is an iterator of tuples. For BigWigs, these tuples are in the format (start: int, end: int, value: float).
                bwdata[cchr]["xstart"] = [x[0] + cstart for x in bw]
                bw = bwf.records(bchr)
                bwdata[cchr]["xend"] = [x[1] + cstart for x in bw]
                bw = bwf.records(bchr)
                bwdata[cchr]["xval"] = [x[2] for x in bw]
            else:
                log.warn(
                    "Bigwig contig %s not found in supplied X axis lengths file" % cchr
                )
        self.export_mapping(outFname, contigs, bwdata, args)

    def export(self, outFname, contigs, bwdata, args):
        """
        for bigwig
        @v1HoloSeq2D for example
        """

        def prepHeader(contigs, args):
            """
            holoSeq output format
            """
            h = ["@%s %s %d" % (holoseq_data.getHap(k), k, contigs[k]) for k in contigs.keys()]
            metah = [
                self.hsId,
                "@@bigwig 1",
                "@@title %s" % args.title,
                "@@datasource %s" % "bigwig",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % args.xclenfile,
            ]

            return metah + h

        hdr = prepHeader(contigs, args)

        with gzip.open(outFname, mode="wb") as ofn:
            ofn.write(str.encode("\n".join(hdr) + "\n"))
            for chr in bwdata.keys():
                for i in range(len(bwdata[chr]["xstart"])):
                    row = str.encode(
                        "%d %d\n" % (bwdata[chr]["xstart"][i], bwdata[chr]["xval"][i])
                    )
                    ofn.write(row)

    def makePanel(self, inFile, pwidth):
        """
        prepare a complete panel for the final display
        """

        def showX(x, y):
            if np.isnan(x):
                s = "Mouse click on image for location"
            else:
                i = bisect_left(h1starts, x)
                chrx = h1names[i - 1]
                offsx = x - h1starts[i - 1]
                s = "%s:%d" % (chrx, offsx)
            str_pane = pn.pane.Str(
                s,
                styles={
                    "font-size": "10pt",
                    "color": "darkblue",
                    "text-align": "center",
                },
                width=pwidth,
            )
            return str_pane

        (hsDims, hapsread, xcoords, ycoords, annos, plotType, metadata, gffdata, hh) = holoseq_data.load(inFile)
        self.rotated = metadata.get("rotated") == "True"
        title = " ".join(metadata["title"])
        haps = []
        print("Read nx=", len(xcoords), "ny=", len(ycoords))
        h1starts = []
        h1names = []
        for i, hap in enumerate(hapsread.keys()):
            haps.append(hap)
            for j, contig in enumerate(hapsread[hap]["cn"]):
                cstart = hapsread[hap]["startpos"][j]
                h1starts.append(cstart)
                h1names.append(contig)
        hap = haps[0]
        log.debug("h1names=%s" % h1names[:20])
        # qtic1 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
        qtic1 = [(h1starts[i], h1names[i]) for i in range(len(h1starts))]
        log.debug("qtic1=%s" % qtic1[:20])
        xax = metadata["xclenfile"][0]
        yax = metadata["yclenfile"][0] + "Bigwig value"
        pafxy = pd.DataFrame.from_dict({xax: xcoords, yax: ycoords})
        taps = hv.streams.Tap(x=0, y=0)
        showloc = pn.bind(showX, x=taps.param.x, y=taps.param.y)
        bigw = pn.pane.HoloViews(
            decimate(hv.Curve(pafxy), streams=[taps])
            .opts(interpolation="steps-pre", color="darkblue")
            .relabel("%s" % title)
            .opts(
                width=pwidth,
                height=300,
                xticks=qtic1,
                xrotation=45,
                fontsize={"xticks": 8, "yticks": 10},
                scalebar=True,
                scalebar_range="x",
                scalebar_location="top_left",
                scalebar_unit=("bp"),
                show_grid=True,
                ylim=(-0.1, 200),
                tools=[
                    "xwheel_zoom",
                    "tap",
                    "xpan",
                    "reset",
                ],
                default_tools=[],
                active_tools=["xwheel_zoom", "tap", "pan"],
            )
        )

        p1 = pn.Column(showloc, bigw)

        return p1, title

