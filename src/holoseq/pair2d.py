# 2d heatmaps from paf or hic tracks in holoSeq
# has functions to read a paf and write a hseq.gz, and to prepare a panel showing that hseq.gz
# using the generic  holoseq_.load function
# need to add the hic code here too

from bisect import bisect_left
from collections import OrderedDict
import copy
import gzip
import logging

import os

import sys


import holoviews as hv
import numpy as np
import pandas as pd
import panel as pn


from holoseq.config import VALID_HSEQ_FORMATS

from holoseq import holoseq_data

from holoseq.rotater import rotater


from holoviews.operation.datashader import (
    rasterize,
    dynspread,
)

# from holoviews.operation.element import apply_when
from holoviews.operation.resample import ResampleOperation2D
from holoviews.operation import decimate

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


hv.extension("bokeh", "matplotlib", width=100)

# Default values suitable for this notebook
decimate.max_samples = 1000
dynspread.max_px = 8
dynspread.threshold = 0.75
ResampleOperation2D.width = 250
ResampleOperation2D.height = 250


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

    def __init__(self):
        """
        
        """
        self.Hsid = VALID_HSEQ_FORMATS[1]

    def convert(self,args):

        haps = []
        yhaps = []
        xcontigs, xhaps = holoseq_data.getContigs(args.xclenfile, args.hap_indicator, 'H1')
        sxcontigs = holoseq_data.contsort(xcontigs, args)
        if args.yclenfile:
            ycontigs, yhaps = holoseq_data.getContigs(args.yclenfile, args.hap_indicator, 'H2')
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
        log.debug("inFile=%s, ftype = %s" % (self.inFname, ps))
        print("sycontigs###",list(sycontigs.keys())[:30], list(sycontigs.values())[:30], )
        self.args = args
        self.haps = haps
        self.xcontigs = sxcontigs
        self.ycontigs = sycontigs
        self.infname = self.inFname
        self.yclenfile = args.yclenfile
        self.xclenfile = args.xclenfile
        self.prepPafGZ()
        if self.isGzip(self.inFname):
            with gzip.open(self.inFname, "rt") as f:
                self.readPAF(f)
        else:
            with open(self.inFname) as f:
                self.readPAF(f)
        self.cis1fio.close()
        self.cis2fio.close()
        self.transfio.close()
        return self.cis0n, self.cis1n, self.trans1n

    def readPAF(self, f):
        # file is text even if gz
        # parse a paf
        # into metadata and coordinates
        self.nrcis1 = self.nrcis2 = self.nrtrans = 0
        for rowi, rows in enumerate(f):
            row = rows.strip().split()
            if len(row) >= 7:
                c1 = row[0]
                c2 = row[5]
                n1 = int(row[2])
                n2 = int(row[7])
                H1 = holoseq_data.getHap(c1, hap_indicator= self.args.hap_indicator)
                H2 = holoseq_data.getHap(c2, hap_indicator= self.args.hap_indicator)
                if H1 != H2:  # trans
                    if H1 == self.haps[0]:  # x is h1 for trans - otherwise ignore
                        x = self.xcontigs[c1] + n1
                        y = self.ycontigs[c2] + n2
                        if self.args.rotate:
                            if y <= x:  # lower triangle
                                x, y = self.rot.rotatecoords(x, y, True)
                            else:
                                x = None
                        if x is not None:
                            row = "%d %d\n" % (x, y)
                            self.transfio.write(row)
                            self.nrtrans += 1
                    else:
                        x = self.xcontigs[c2] + n2
                        y = self.ycontigs[c1] + n1
                        if self.args.rotate:
                            if y <= x:  # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row ="%d %d\n" % (x, y)
                            self.transfio.write(row)
                            self.nrtrans += 1
                else:  # cis.
                    if H1 == self.haps[0]:
                        x = self.xcontigs[c1] + n1
                        y = self.xcontigs[c2] + n2
                        if self.args.rotate:
                            if y <= x:  # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row = "%d %d\n" % (x, y)
                            self.cis1fio.write(row)
                            self.nrcis1 += 1
                    else:
                        x = self.ycontigs[c1] + n1
                        y = self.ycontigs[c2] + n2
                        if self.args.rotate:
                            if y <= x:  # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row = "%d %d\n" % (x, y)
                            self.cis2fio.write(row)
                            self.nrcis2 += 1
        log.debug(
            "nrcis1=%d, nrcis2=%d, nrtrans=%d" % (self.nrcis1, self.nrcis2, self.nrtrans)
        )

    def isGzip(self, inFname):
        with gzip.open(inFname, "r") as fh:
            try:
                fh.read(1)
                return True
            except gzip.BadGzipFile:
                log.info("inFname %s is not a gzip so will read as text" % inFname)
                return False

    def prepPafGZ(self):
        """
        write a holoseq paired gzip metadata and coordinates
        for plotting by the makepanel function
        """

        def prepHeader(
            haps,
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
            Depending on the configuration,Mashmap sequence similarity paf output pairs may both be coordinates on contigs from the same haplotype, or from different ones.
            HiC data usually has two possible haplotypes.
            So typically 3 kinds of pairs - same haplotype termed "cis" and different termed "trans"
            That's not how those words are usually used but it seems a natural extension for HiC

            prepare the three potentially needed gzip output channels
            """
            h = [
                "@%s %s %d" % (holoseq_data.getHap(k, hap_indicator= self.args.hap_indicator), k, xcontigs[k])
                for k in xcontigs
            ]
            if len(haps) > 1:
                h += [
                    "@%s %s %d" % (holoseq_data.getHap(k, hap_indicator= self.args.hap_indicator), k, ycontigs[k])
                    for k in ycontigs.keys()
                ]
            metah = [self.Hsid,
                "@@class pair2d",
                "@@title %s" % self.args.title + subtitle,
                "@@datasource %s" % "paf",
                "@@datafile %s" % self.infname,
                "@@refURI %s" % self.args.refURI,
                "@@xclenfile %s" % self.xclenfile,
                "@@yclenfile %s" % self.yclenfile,
                "@@axes %s" % ax,
                "@@rotated %s" % self.args.rotate,
            ]
            menc = metah + h
            outs = '\n'.join(menc)
            outf.write(outs)

        if self.args.cis1outpath:
            fn1 = self.args.cis1outpath
        else:
            fn1 = "%s_cis%s_hseq.gz" % (self.inFname, self.haps[0])
            if self.args.rotate:
                fn1 = "%s_rotated_cis%s_hseq.gz" % (self.inFname, self.haps[0])
        self.cis0n = fn1
        self.cis1fio = gzip.open(fn1, mode="wt")
        prepHeader(
            self.haps[0],
            self.xcontigs,
            self.xcontigs,
            self.args,
            self.cis1fio,
            " Pairs on %s" % self.haps[0],
            xclenfile=self.args.xclenfile,
            yclenfile=self.args.xclenfile,
            ax=self.haps[0],
        )
        
        if self.args.cis2outpath:
            fn2 = self.args.cis2outpath
        else:
            fn2 = "%s_cis%s_hseq.gz" % (self.inFname, self.haps[1])
            if  self.args.rotate:
                fn2 = "%s_rotated_cis%s_hseq.gz" % (self.inFname, self.haps[1])
        self.cis2fio = gzip.open(fn2, mode="wt")
        self.cis1n = fn2
        prepHeader(
            self.haps[1],
            self.ycontigs,
            self.ycontigs,
            self.args,
            self.cis2fio,
            " Pairs on %s" % self.haps[1],
            xclenfile=self.args.yclenfile,
            yclenfile=self.args.yclenfile,
            ax=self.haps[1],
        )
        if self.args.transoutpath:
            fn3 = self.args.transoutpath
        else:
            fn3 = "%s_trans_hseq.gz" % (self.inFname)
            if self.args.rotate:
                fn3 = "%s_rotated_trans_hseq.gz" % (self.inFname)
        self.transfio = gzip.open(fn3, mode="wt")
        self.trans1n = fn3
        
        prepHeader(
            self.haps,
            self.xcontigs,
            self.ycontigs,
            self.args,
            self.transfio,
            " Pairs on different haplotypes",
            xclenfile=self.args.xclenfile,
            yclenfile=self.args.yclenfile,
            ax="BOTH",
        )


    def makePanel(self,infname, pwidth=1000):
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
                    s = (
                        "Rotated X axis genome %s:%d Rotated Y axis genome %s:%d x %d y %d xur %d yur %d"
                        % (
                            chrx,
                            offsx,
                            chry,
                            offsy,
                            x,
                            y,
                            xur,
                            yur,
                        )
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
                    s = (
                        "X axis genome %s:%d Y axis genome %s:%d x %d y %d xur %d yur %d"
                        % (
                            chrx,
                            offsx,
                            chry,
                            offsy,
                            x,
                            y,
                            xur,
                            yur,
                        )
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


        self.pwidth = pwidth
        (num_dimensions,
            haploids,
            x_coords,
            y_coords,
            annotations,
            plot_type,
            metadata,
            gff_data,
            hh,
            rotated) = holoseq_data.load(infname)
        self.rotated = rotated
        rot = rotater(max(x_coords), max(y_coords))
        title = " ".join(metadata["title"])
        hqstarts = OrderedDict()
        haps = []
        h1starts = []
        h1names = []
        h2starts = []
        h2names = []

        for i, hap in enumerate(haploids.keys()):
            haps.append(hap)
            cnames = haploids[hap]["contig_names"]
            cstarts = haploids[hap]['starts']
            hqstarts[hap] = OrderedDict(zip(cnames,cstarts))
            if i == 0: # first haplotype
                h1starts = copy.copy(cstarts)
                h1names = copy.copy(cnames)
            elif i == 1:
                h2starts = copy.copy(cstarts)
                h2names = copy.copy(cnames)
            else:
                log.debug('i=%d but code only expects 2 max.' % i)
        hap = hh[0]
        if len(h2starts) == 0:
            h2starts = copy.copy(h1starts)
            h2names = copy.copy(h1names)
            log.warn("only one haplotype read for %s" % title)
        qtic1 = [(h1starts[i], h1names[i]) for i in range(len(h1starts))]
        hap = hh[1]
        if self.rotated:  # yaxis makes no real sense
            qtic2 = [(0, "")]
        else:
            qtic2 = [(h2starts[i], h2names[i]) for i in range(len(h2starts))]
        #print('*****qtic2: %s' % qtic2)
        # once the pairs have been read and mapped into a grid, the code
        # below does the plotting.
        # it can be copied, edited to suit your needs and
        # run repeatedly without waiting for the  holoseq_ to be mapped.
        xcf = os.path.splitext(metadata["xclenfile"][0])[0]
        ycf = "Y:" + os.path.splitext(metadata["yclenfile"][0])[0]
        pafxy = pd.DataFrame.from_dict({xcf: x_coords, ycf: y_coords})
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
                    xticks=qtic1,
                    yticks=qtic2,
                    xrotation=45,
                    width=pwidth,
                    height=pwidth,
                    # fontsize={"xticks": 8, "yticks": 8},
                    fontscale=1,
                    scalebar=True,
                    scalebar_range="x",
                    scalebar_location="top_left",
                    scalebar_unit=("bp"),
                    show_grid=True,
                )
            ),
        )
        return p1, title