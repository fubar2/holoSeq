# 2d heatmaps from paf or hic tracks in holoSeq
# has functions to read a paf and write a hseq.gz, and to prepare a panel showing that hseq.gz 
# using the generic  holoseq_.load function
# need to add the hic code here too

from bisect import bisect_left
from collections import OrderedDict
import gzip
import io
import logging
import numpy as np
import os

import sys



import holoviews as hv
import pandas as pd
import panel as pn



from config import VALID_HSEQ_FORMATS

from holoseq import holoseq_data

from rotater import rotater


from holoviews.operation.datashader import (
    rasterize,
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
log = logging.getLogger("pair2d")


class pair2d:
    """
    updated to stream row at a time -
    slower but no room for the entire output if the input is 60GB+

    paf to xy and axis metadata
    Assumes pairs of points representing HiC contact pairs, or Mashmap sequence similarity hits. HiC  holoseq_ typically comes from pairs of haplotypes and is used to help assemble all
    the contigs into chromosomes

    These coordinate pairs are of 3 types - both ends on one or the other reference sequence, or one end on each for HiC pairs.
    Let's call these cis if the same reference and trans if between the two references. The ones in cis are likely to be consistent-ish between the two haplotypes, but the trans ones turn
    out to be very interesting...

    For HiC  holoseq_, there should really be little difference - whether a haplotype contacts another haplotype or itself depends on the 3D folding and since the haplotypes are wound together into
    a helix, then the whole total contig length - about 2 meters for mammals are folded up into a 10μ^3 ball.

    python holoSeq_prepare_gz.py --inFile mUroPar1H1H2.paf --xclenfile mUroPar1H1suffix.len --yclenfile mUroPar1H2suffix.len --contig_sort VGPname --hap_indicator Suffix
    """

    def __init__(self, inFname, args, xcontigs, ycontigs, haps, xwidth, ywidth):
        """
        if rotating paf line at a time, better to pre-calculate the constants once rather than once for each point
        """
        self.inFname = inFname
        self.args = args
        self.rot = rotater(xwidth, ywidth)
        self.rotate = args.rotate
        self.rot.onepointRot = True
        self.rot.origin = (0, ywidth)
        self.outFprefix = inFname
        self.xcontigs = xcontigs
        self.ycontigs = ycontigs
        self.haps = haps
        # have the axes set up so prepare the three plot x/y vectors
        # for a second pass to calculate all the coordinates.
        # adding tooltips just does not scale so abandoned - see the tooltip old version
        hsId = VALID_HSEQ_FORMATS[1]
        self.inFname = inFname
        self.prepPafGZ(hsId, haps, xcontigs, ycontigs)
        if self.isGzip(inFname):
            with gzip.open(inFname, "rt") as f:
                self.readPAF(f)
        else:
            with open(inFname) as f:
                self.readPAF(f)
        self.cis1f.close()
        self.cis2f.close()
        self.transf.close()


    def makePanel(self, inFile, pwidth):
        """
        prepare a complete panel for the final display
        """

        def showTap(x, y, rot, cxstarts, cxnames, cystarts, cynames, rotated):
            if np.isnan(x) or np.isnan(y):
                s = "Mouse click on image for location"
            else:
                chrx = "Out of range"
                offsx = 0
                chry = "Out of range"
                offsy = 0
                xur = 0
                yur = 0
                if rotated == "True":
                    xur, yur = rot.unrotatecoords(xr=x, yr=y)
                    i = bisect_left(cxstarts, xur)
                    if i > 0 and i <= len(cxnames):
                        chrx = cxnames[i - 1]
                        offsx = xur - cxstarts[i - 1]
                    i = bisect_left(cystarts, yur)
                    if i > 0 and i <= len(cynames):
                        chry = cynames[i - 1]
                        offsy = yur - cystarts[i - 1]
                    s = "Rotated X axis genome %s:%d Rotated Y axis genome %s:%d x %d y %d xur %d yur %d" % (
                        chrx,
                        offsx,
                        chry,
                        offsy,
                        x,
                        y,
                        xur,
                        yur,
                    )
                else:
                    i = bisect_left(cxstarts, x)
                    if i > 0 and i <= len(cxnames):
                        chrx = cxnames[i - 1]
                        offsx = x - cxstarts[i - 1]
                    i = bisect_left(cystarts, y)
                    if i > 0 and i <= len(cynames):
                        chry = cynames[i - 1]
                        offsy = y - cystarts[i - 1]
                    s = "X axis genome %s:%d Y axis genome %s:%d x %d y %d xur %d yur %d" % (
                chrx,
                offsx,
                chry,
                offsy,
                x,
                y,
                xur,
                yur,
            )
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

        (hsDims, hapsread, xcoords, ycoords, annos, plotType, metadata, gffdata, hh) = (
            holoseq_data.load(inFile)
        )
        self.rotated = metadata.get("rotated", [False,])[0]
        print('rotated', self.rotated)
        rot = rotater(max(xcoords), max(ycoords))
        title = " ".join(metadata["title"])
        hqstarts = OrderedDict()
        haps = []
        print("Read nx=", len(xcoords), "ny=", len(ycoords))
        h1starts = []
        h1names = []
        h2starts = []
        h2names = []
        for i, hap in enumerate(hapsread.keys()):
            haps.append(hap)
            hqstarts[hap] = OrderedDict()
            for j, contig in enumerate(hapsread[hap]["cn"]):
                cstart = hapsread[hap]["startpos"][j]
                hqstarts[hap][contig] = cstart
                if i == 0:
                    h1starts.append(cstart)
                    h1names.append(contig)
                else:
                    h2starts.append(cstart)
                    h2names.append(contig)
        hap = hh[0]
        if len(h2starts) == 0:
            h2starts = h1starts
            h2names = h1names
            log.warn("only one haplotype read for %s" % title)
        qtic1 = [(h1starts[i], h1names[i]) for i in range(len(h1starts))]
        hap = hh[1]
        if self.rotated: # yaxis makes no real sense
            qtic2 = [(0,'')]
        else:
            qtic2 = [(h2starts[i], h2names[i]) for i in range(len(h2starts))]
        # once the pairs have been read and mapped into a grid, the code
        # below does the plotting.
        # it can be copied, edited to suit your needs and
        # run repeatedly without waiting for the  holoseq_ to be mapped.
        xcf = os.path.splitext(metadata["xclenfile"][0])[0]
        ycf = "Y:" + os.path.splitext(metadata["yclenfile"][0])[0]
        print("xcf", xcf, "ycf", ycf)
        pafxy = pd.DataFrame.from_dict({xcf: xcoords, ycf: ycoords})
        pafp = hv.Points(pafxy, kdims=[xcf, ycf])

        # apply_when(pafp, operation=rasterize, predicate=lambda x: len(x) > 5000)
        stream = hv.streams.Tap(x=0, y=0)
        ax = metadata.get("axes", [None])[0]
        log.debug("axes = %s" % ax)
        if ax == "BOTH":
            showloc = pn.bind(
                showTap,
                x=stream.param.x,
                y=stream.param.y,
                rot=rot,
                cxnames=h1names,
                cxstarts=h1starts,
                cynames=h2names,
                cystarts=h2starts,
                rotated=self.rotated,
            )
        elif ax == haps[0]:
            showloc = pn.bind(
                showTap,
                x=stream.param.x,
                y=stream.param.y,
                rot=rot,
                cxnames=h1names,
                cxstarts=h1starts,
                cynames=h1names,
                cystarts=h1starts,
                rotated=self.rotated,
            )
        elif ax == haps[1]:
            showloc = pn.bind(
                showTap,
                x=stream.param.x,
                y=stream.param.y,
                rot=rot,
                cxnames=h2names,
                cxstarts=h2starts,
                cynames=h2names,
                cystarts=h2starts,
                rotated=self.rotated,
            )
        else:
            log.warn("ax = %s for title = %s - cannot assign axes" % (ax, title))
            sys.exit(2)
        # an alternative but can't get a stream in there..nice to have control over the resample_when but.
        # dat.hvplot(kind="scatter", x="x", y="y", color="maroon", rasterize=True, resample_when=200, cnorm='log', padding=(0, 0.1), cmap="inferno",
        #   min_height=700, autorange='y', title="Datashader Rasterize", colorbar=True, line_width=2 ,marker="x" )

        p1 = pn.Column(
            showloc,
            pn.pane.HoloViews(
                dynspread(rasterize(pafp), streams=[stream])
                .relabel("%s" % title)
                .opts(
                    cmap="inferno",
                    cnorm="log",
                    colorbar=True,
                    shared_axes=False,
                    width=self.pwidth,
                    height=self.pwidth,
                    xticks=qtic1,
                    yticks=qtic2,
                    xrotation=45,
                    fontsize={"xticks": 5, "yticks": 5},
                    tools=["tap"],
                    scalebar=True,
                    scalebar_range="x",
                    scalebar_location="bottom_left",
                    scalebar_unit=("bp"),
                    show_grid=True,
                )
            ),
        )
        return p1, title

    def readPAF(self, f):
        # might be a gz
        self.nrcis1 = self.nrcis2 = self.nrtrans = 0
        for rowi, rows in enumerate(f):
            row = rows.strip().split()
            if len(row) > 7:
                c1 = row[0]
                c2 = row[5]
                n1 = int(row[2])
                n2 = int(row[7])
                H1 =  holoseq_data.getHap(c1)
                H2 =  holoseq_data.getHap(c2)
                if H1 != H2:  # trans
                    if H1 == self.haps[0]:  # x is h1 for trans - otherwise ignore
                        x = self.xcontigs[c1] + n1
                        y = self.ycontigs[c2] + n2
                        if self.rotate:
                            if y <= x:  # lower triangle
                                x, y = self.rot.rotatecoords(x, y, True)
                            else:
                                x = None
                        if x is not None:
                            row = str.encode("%d %d\n" % (x, y))
                            self.transf.write(row)
                            self.nrtrans += 1
                    else:
                        x = self.xcontigs[c2] + n2
                        y = self.ycontigs[c1] + n1
                        if self.rotate:
                            if y <= x:  # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row = str.encode("%d %d\n" % (x, y))
                            self.transf.write(row)
                            self.nrtrans += 1
                else:  # cis
                    if H1 == self.haps[0]:
                        x = self.xcontigs[c1] + n1
                        y = self.xcontigs[c2] + n2
                        if self.rotate:
                            if y <= x:  # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row = str.encode("%d %d\n" % (x, y))
                            self.cis1f.write(row)
                            self.nrcis1 += 1
                    else:
                        x = self.ycontigs[c1] + n1
                        y = self.ycontigs[c2] + n2
                        if self.rotate:
                            if y <= x:  # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row = str.encode("%d %d\n" % (x, y))
                            self.cis2f.write(row)
                            self.nrcis2 += 1
        log.debug("nrcis1=%d, nris2=%d, nrtrans=%d" % (self.nrcis1, self.nrcis2, self.nrtrans))

    def isGzip(self, inFname):
        with gzip.open(inFname, "r") as fh:
            try:
                fh.read(1)
                return True
            except gzip.BadGzipFile:
                log.info("inFname %s is not a gzip so will read as text" % inFname)
                return False

    def prepPafGZ(self, hsId, haps, xcontigs, ycontigs):
        """
        A pairwise 2d plot involves a genome for each axis.
        A paf or HiC track file contains pairs of coordinates (contig, offset).
        Both coordinates might be on the same genome (termed "cis" in this 2 haplotype context).
        It is also possible that the coordinates are on two separate genomes (termed "trans" here). 
        Pairs are typically binned to downscale the plot density for HiC pairs, and this distinction is usually ignored - all pairs starting or ending in the window 
        are counted, no matter which haplotype is involved. This may help smooth the  holoseq_ and make it easier to visualise, but it throws away a lot of information
        For individual points, it may be interesting to see the raw  holoseq_ but that requires plot axes 
        to correspond precisely to the genome lengths file, for the coordinates to be correct.
        """

        def prepHeader(
            haps,
            hsId,
            xcontigs,
            ycontigs,
            args,
            outf,
            subtitle,
            xclenfile,
            yclenfile,
            ax,
        ):
            """
            holoSeq output format - prepare gzip output channels
            """
            h = [
                "@%s %s %d" % ( holoseq_data.getHap(k), k, xcontigs[k]) for k in xcontigs.keys()
            ]
            if len(haps) > 1:
                h += [
                    "@%s %s %d" % ( holoseq_data.getHap(k), k, ycontigs[k])
                    for k in ycontigs.keys()
                ]
            metah = [
                hsId,
                "@@class heatmap",
                "@@title %s" % self.args.title + subtitle,
                "@@datasource %s" % "paf",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % self.args.refURI,
                "@@xclenfile %s" % xclenfile,
                "@@yclenfile %s" % yclenfile,
                "@@axes %s" % ax,
                "@@rotated %s" % self.rotate,
            ]

            outs = "\n".join(metah + h) + "\n"
            outf.write(str.encode(outs))

        fn1 = "%s_cis%s_hseq.gz" % (self.inFname, haps[0])
        if self.rotate:
            fn1 = "%s_rotated_cis%s_hseq.gz" % (self.inFname, haps[0])
        self.cis0n = fn1
        f1 = gzip.open(fn1, mode="wb")
        self.cis1f = io.BufferedWriter(f1, buffer_size=1024 * 1024)
        prepHeader(
            haps[0],
            hsId,
            xcontigs,
            xcontigs,
            self.args,
            self.cis1f,
            " Pairs on %s" % haps[0],
            xclenfile=self.args.xclenfile,
            yclenfile=self.args.xclenfile,
            ax=haps[0],
        )
        fn2 = "%s_cis%s_hseq.gz" % (self.inFname, haps[1])
        if self.rotate:
            fn2 = "%s_rotated_cis%s_hseq.gz" % (self.inFname, haps[1])
        self.cis1n = fn2
        f2 = gzip.open(fn2, mode="wb")
        self.cis2f = io.BufferedWriter(f2, buffer_size=1024 * 1024)
        prepHeader(
            haps[1],
            hsId,
            ycontigs,
            ycontigs,
            self.args,
            self.cis2f,
            " Pairs on %s" % haps[1],
            xclenfile=self.args.yclenfile,
            yclenfile=self.args.yclenfile,
            ax=haps[1],
        )
        fn3 = "%s_trans_hseq.gz" % (self.inFname)
        self.transn = fn1
        if self.rotate:
            fn3 = "%s_rotated_trans_hseq.gz" % (self.inFname)
        f3 = gzip.open(fn3, mode="wb")
        self.transf = io.BufferedWriter(f3, buffer_size=1024 * 1024)
        prepHeader(
            haps,
            hsId,
            xcontigs,
            ycontigs,
            self.args,
            self.transf,
            " Pairs on different haplotypes",
            xclenfile=self.args.xclenfile,
            yclenfile=self.args.yclenfile,
            ax="BOTH",
        )
