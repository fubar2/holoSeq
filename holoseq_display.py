# see https://github.com/fubar2/holoSeq
# pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh
#
# This is a generic plot re-creator for holoSeq hseq compressed coordinate data.
# It presents interactive plots at scale using panel.
# Plot coordinates and axis metadata are prepared from genomic data such as PAF or bigwig files,
# using companion converters. They output all the information needed to re-create an interactive plot, as coordinates and axis metadata.
# The intended use case is making interactive specialised browsers for VGP and other genome assemblies and their annotation
# easily browsed and explored using a web browser.
# Annotation files with 10s of millions of data points can be zoomed smoothly down from the whole assembly scale, down to individual points,
# on a laptop.
# Ross Lazarus October 2024

import argparse
from bisect import bisect_left
from collections import OrderedDict
import gzip
import logging
import numpy as np
import os

import holoviews as hv
import panel as pn
import pandas as pd

from holoviews.operation.datashader import (
    rasterize,
    dynspread,
    datashade,
)
from holoviews.operation.element import apply_when
from holoviews.operation.resample import ResampleOperation2D
from holoviews.operation import decimate


logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_display")

from holoviews import opts
from holoviews.operation.datashader import datashade, rasterize, shade, dynspread, spread


hv.extension('bokeh','matplotlib', width=100)

# Default values suitable for this notebook
decimate.max_samples=1000
dynspread.max_px=20
dynspread.threshold=0.5
ResampleOperation2D.width=250
ResampleOperation2D.height=250


# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/holoSeqtest.gz"
holoSeqHeaders = ["@v1HoloSeq1D", "@v1HoloSeq2D"]
hv.extension("bokeh")
pn.extension()

dynspread.max_px = 7
dynspread.threshold = 0.6


def xportHtml(fname, hObj):
    "save a holoview object to an interactive but not adaptive scaling HTML page"
    hv.save(filename=fname, obj=hObj)


class holoSeq_maker:
    """
    returns panels reconstructed from hseq data - coordinates with axis and other metadata
    """

    def __init__(self, width):
        """ """
        self.pwidth = width

    def xportHtml(self, fname, hObj):
        "save a holoview object to an interactive but not adaptive scaling HTML page"
        hv.save(filename=fname, obj=hObj)

    def import_holoSeq_data(self, inFile):
        """
        reverse process of dumping the data in holoSeq format from a converter
        see https://github.com/fubar2/holoSeq/blob/main/HoloSeqOverview.md
        """
        haps = {}
        hh = []
        xcoords = []
        ycoords = []
        annos = []
        metadata = {}
        gffdata = []
        hsDims = None
        plotType = None
        with gzip.open(inFile, "rt") as f:
            for i, trow in enumerate(f):
                if i == 0:
                    hseqformat = trow.split()[0].strip()
                    if hseqformat not in holoSeqHeaders:
                        print(
                            "Supplied input",
                            inFile,
                            "has a first row =",
                            trow,
                            "so is not a valid holoSeq input file",
                        )
                        print("First row must start with one of these:", holoSeqHeaders)
                        return
                    hsDims = holoSeqHeaders.index(hseqformat) + 1
                    if hsDims == 1:
                        plotType = "bar"
                        if len(hseqformat) > 1:
                            plotType = hseqformat[1].strip()
                elif trow.startswith("@"):
                    row = trow[1:]
                    if row.startswith("@"):
                        srow = row[1:].split()
                        metadata[srow[0]] = srow[1:]
                        if srow[0] == "GFF":
                            isGFF = True
                    else:
                        srow = row.split()
                        if hsDims == 2:
                            if len(srow) >= 3:
                                hap, cname, cstart = srow[:3]
                                if not haps.get(hap, None):
                                    print("adding hap", hap)
                                    hh.append(hap)
                                    haps[hap] = {"cn": [], "startpos": []}
                                haps[hap]["cn"].append(cname.strip())
                                haps[hap]["startpos"].append(int(cstart.strip()))
                            else:
                                print(
                                    "Supplied input",
                                    inFile,
                                    "at row",
                                    i,
                                    "=",
                                    row,
                                    "lacking the required reference name, contig name and contig length. Not a valid holoSeq input file",
                                )
                                return
                        else:
                            if len(srow) >= 2:
                                hap, cname, cstart = srow[:3]
                                if not haps.get(hap, None):
                                    print("adding hap", hap)
                                    hh.append(hap)
                                    haps[hap] = {"cn": [], "startpos": []}
                                haps[hap]["cn"].append(cname.strip())
                                haps[hap]["startpos"].append(int(cstart.strip()))
                            else:
                                print(
                                    "Supplied input",
                                    inFile,
                                    "at row",
                                    i,
                                    "=",
                                    row,
                                    "lacking the required reference name, contig name and contig length. Not a valid holoSeq input file",
                                )
                                return
                else:  # not header row
                    srow = [x.strip() for x in trow.strip().split()]
                    lrow = len(srow)
                    if hsDims == 2:
                        if lrow < 2:
                            print(
                                "At row",
                                i,
                                "Supplied 2D input",
                                inFile,
                                "has",
                                trow,
                                "which does not parse into x and ycoordinates and optional annotation",
                            )
                            return
                        else:
                            if srow[0].isdigit() and srow[1].isdigit():
                                xcoords.append(int(srow[0]))
                                ycoords.append(int(srow[1]))
                                if lrow > 2:
                                    annos.append(srow[2:])
                            else:
                                print(
                                    "At row",
                                    i,
                                    "Supplied 2D input",
                                    inFile,
                                    "has",
                                    trow,
                                    "which does not parse into x and ycoordinates and optional annotation",
                                )
                                return
                    else:
                        if isGFF:
                            gffdata.append(
                                srow
                            )  # mrna XP_026254554.1 SUPER_5 1006698964 1006703995 100 100 + 3226 1006703993

                        else:
                            if srow[0].isdigit():
                                xcoords.append(int(srow[0]))
                                if lrow > 1:
                                    ycoords.append(int(srow[1]))
                                if lrow > 2:
                                    annos.append(srow[2:])
                            else:
                                print(
                                    "At row",
                                    i,
                                    "Supplied 1D input",
                                    inFile,
                                    "has",
                                    trow,
                                    "which does not parse into x coordinate and optional annotation",
                                )
                                return

        return (hsDims, haps, xcoords, ycoords, annos, plotType, metadata, gffdata)

    def makePafPanel(self, inFile, pwidth):
        """
        prepare a complete panel for the final display
        """

        def showH1(x, y):
            if np.isnan(x) or np.isnan(y):
                s = "Mouse click on image for location"
            else:
                i = bisect_left(h1starts, x)
                chrx = h1names[i - 1]
                offsx = x - h1starts[i - 1]
                i = bisect_left(h1starts, y)
                chry = h1names[i - 1]
                offsy = y - h1starts[i - 1]
                s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
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

        def showH2(x, y):
            if np.isnan(x) or np.isnan(y):
                s = "Mouse click on image for location"
            else:
                i = bisect_left(self.h2starts, x)
                chrx = h2names[i - 1]
                offsx = x - h2starts[i - 1]
                i = bisect_left(h2starts, y)
                chry = h2names[i - 1]
                offsy = y - h2starts[i - 1]
                s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
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

        def showTrans(x, y):
            if np.isnan(x) or np.isnan(y):
                s = "Mouse click on image for location"
            else:
                i = bisect_left(h1starts, x)
                chrx = h1names[i - 1]
                offsx = x - h1starts[i - 1]
                i = bisect_left(h2starts, y)
                chry = h2names[i - 1]
                offsy = y - h2starts[i - 1]
                s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
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

        (hsDims, hapsread, xcoords, ycoords, annos, plotType, metadata, gffdata) = (
            self.import_holoSeq_data(inFile)
        )
        title = " ".join(metadata["title"])
        hqstarts = OrderedDict()
        haps = []
        print("Read nx=", len(xcoords), "ny=", len(ycoords))
        h1starts = []
        h1names = []
        for i, hap in enumerate(hapsread.keys()):
            haps.append(hap)
            hqstarts[hap] = OrderedDict()
            for j, contig in enumerate(hapsread[hap]["cn"]):
                cstart = hapsread[hap]["startpos"][j]
                hqstarts[hap][contig] = cstart
                if i == 0:
                    h1starts.append(cstart)
                    h1names.append(contig)
        hap = haps[0]

        qtic1 = [(h1starts[i], h1names[i]) for i in range(len(h1starts))]
        qtic2 = qtic1
        isTrans = False
        if (
            hsDims == "2"
        ):  # may be one or two distinct haplotype identifiers - H1 +/- H2
            if len(haps) > 1:
                hap = haps[1]
                isTrans = True
                h2starts = [hqstarts[hap][x] for x in hqstarts[hap].keys()]
                h2names = list(hqstarts[hap].keys())
                qtic2 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
                print("h2names=", h1names[:20], len(h2names))
        if not isTrans:
            print(len(h1names), "h1names=", h1names[:20])
        # can take the np.tril or filter the upper triangle while processing pairs
        # and rotate so the diagonal becomes the x axis but need some kind of
        # sideways scroller to work right
        # xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
        # pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
        # --------------------cut here---------------------------------
        # once the pairs have been read and mapped into a grid, the code
        # below does the plotting.
        # it can be copied, edited to suit your needs and
        # run repeatedly without waiting for the data to be mapped.
        xcf = os.path.splitext(metadata["xclenfile"][0])[0]
        ycf = os.path.splitext(metadata["yclenfile"][0])[0]
        print("xcf", xcf, "ycf", ycf)
        pafxy = pd.DataFrame.from_dict({xcf: xcoords, ycf: ycoords})
        pafp = hv.Points(pafxy, kdims=[xcf, ycf])

        stream = hv.streams.Tap(x=0, y=0)
        if isTrans:
            showloc = pn.bind(showTrans, x=stream.param.x, y=stream.param.y)
        else:
            showloc = pn.bind(showH1, x=stream.param.x, y=stream.param.y)
        # to rotate so the diagonal becomes the x axis
        # xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
        # pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
        # bisect.bisect_left(a, x, lo=0, hi=len(a), *, key=None)
        # prepare and show the 3 plots

        p1 = pn.Column(
            showloc,
            pn.pane.HoloViews(
                dynspread(rasterize(pafp), streams=[stream])
                .relabel("%s" % title)
                .opts(
                    shared_axes=False,
                    cmap="inferno",
                    cnorm="log",
                    colorbar=True,
                    width=self.pwidth,
                    height=self.pwidth,
                    xticks=qtic1,
                    yticks=qtic2,
                    xrotation=45,
                    fontsize={"xticks": 5, "yticks": 5},
                    tools=["tap"],
                    scalebar=True,
                    scalebar_range="x",
                    scalebar_location="top_left",
                    scalebar_unit=("bp"),
                    show_grid=True,
                )
            ),
        )
        return p1, title

    def makeBWPanel(self, inFile, pwidth):
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

        (hsDims, hapsread, xcoords, ycoords, annos, plotType, metadata, gffdata) = (
            self.import_holoSeq_data(inFile)
        )
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
        print("h1names=", h1names[:20])
        # qtic1 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
        qtic1 = [(h1starts[i], h1names[i]) for i in range(len(h1starts))]
        print("qtic1=", qtic1[:20])
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

    def makeGFFPanel(self, inFile, pwidth):
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
                xi = bisect_left(segs['x1'], x)
                xtarget = segs['target'][xi]
                yreal = segs['y1'][xi]
                s += ' x %s xy %d real y %d' % (xtarget, yreal, y )

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

        (hsDims, hapsread, xcoords, ycoords, annos, plotType, metadata, gffdata) = self.import_holoSeq_data(inFile)
        segs = {
            "x1": [],
            "x2": [],
            "y1": [],
            "y2": [],
            "target": [],
            "colour": [],
            "thickness": [],
            "alpha": []
        }
        """
cds XP_026238700.1 1401967516 1401967635 100 100 - 204
mrna XP_026248570.1 SUPER_3H1 531341254 531595863 100 100 + 1102 -1
cds XP_026248570.1 531341254 531341334 100 100 + 134
        """
        mthick = 4
        cdthick = 50
        for i, rows in enumerate(gffdata):
            if rows[0].lower() == "mrna":
                (kind, targ,  contig, startp, endp, y1, y2, strand, score) = rows[:10]
                startp = int(startp)
                endp = int(endp)
                y1 = int(y1)
                y2 = int(y2)
                colr = "black"
                segs["target"].append(targ)
                segs["x1"].append(startp)
                segs["x2"].append(endp)
                segs["y1"].append(y1)
                segs["y2"].append(y2)
                segs["colour"].append(colr)
                segs["thickness"].append(mthick)
                segs["alpha"].append(1.0)
                #segs["stopc"].append(stopc)
            elif rows[0].lower() == "cds": #  f"cds {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                (kind, targ, contig, startp, endp , y1, y2, strand, score) = rows[:10]
                startp = int(startp)
                endp = int(endp)
                y = int(y1)
                colr = "blue"
                if strand == "-":
                    colr = "maroon"        
                segs["target"].append(targ)
                segs["x1"].append(startp)
                segs["x2"].append(endp)
                segs["y1"].append(y)
                segs["y2"].append(y)
                segs["colour"].append(colr)
                segs["thickness"].append(cdthick)
                segs["alpha"].append(0.5)
        xmin = min(segs["x1"])
        xmax = max(segs["x2"])
        ymin = min(segs["y1"])
        ymax = max(segs["y1"])
        title = " ".join(metadata["title"])
        haps = []
        print("GFF rows read =", len(gffdata))
        h1starts = []
        h1names = []
        qtic1 = []
        for i, hap in enumerate(hapsread.keys()):
            haps.append(hap)
            for j, contig in enumerate(hapsread[hap]["cn"]):
                cstart = hapsread[hap]["startpos"][j]
                h1starts.append(cstart)
                h1names.append(contig)
                qtic1.append((cstart, contig))
        hap = haps[0]
        # print("h1names=", h1names[:20])
        # qtic1 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
        # print("qtic1=", qtic1[:20])
        gffp = hv.Segments(
            segs, ['x1', 'y1', 'x2', 'y2'], vdims=["target", "colour", "thickness", "alpha"]
            )

        gffp.opts(color="colour",
            line_width="thickness",
            alpha="alpha",
            width=pwidth,
            height=300,
            xticks=qtic1,
            xrotation=45,
            fontsize={"xticks": 8, "yticks": 10},
            show_grid=True,
            #ylim=(50 , 105),
            tools=[
                "xwheel_zoom",
                "box_zoom",
                "tap",
                "xpan",
                "reset",
            ],
            default_tools=[],
            active_tools=["xwheel_zoom", "tap", "xpan"],
            )
        apply_when(gffp, operation=datashade, predicate=lambda x: x > 1000)

        taps = hv.streams.Tap(source=gffp, x=0, y=0)
        showloc = pn.bind(showX, x=taps.param.x, y=taps.param.y)
        gp = pn.pane.HoloViews(gffp)
        p = pn.Column(showloc, gp)

        return p, title


parser = argparse.ArgumentParser(description="", epilog="")
parser.add_argument(
    "--inFile",
    help="gzipped hseq coordinates and contigs",
    default="mUroPar1_cis1.hseq.gz",
    nargs="+",
)
parser.add_argument(
    "--size", help="Display size in pixels. Default is 800", default=1000
)
parser.add_argument("--version", "-V", action="version", version="0.1")
args = parser.parse_args()
pwidth = int(args.size)
hsm = holoSeq_maker(pwidth)
for i, infile in enumerate(args.inFile):
    print("Infile = ", infile)
    if "bw.hseq.gz" in infile:
        p1, title = hsm.makeBWPanel(infile, pwidth)
    elif "gff.hseq.gz" in infile:
        p1, title = hsm.makeGFFPanel(infile, pwidth)
    else:
        p1, title = hsm.makePafPanel(infile, pwidth)
    if i == 0:
        outp = p1
    else:
        outp = outp + p1
pn.Row(outp).servable(title=title)
