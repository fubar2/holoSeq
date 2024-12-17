import gzip
import io
import logging

from config import VALID_HSEQ_FORMATS
import data


from rotater import rotater


logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("pair2d")


class pafConvert:
    """
    updated to stream row at a time -
    slower but no room for the entire output if the input is 60GB+

    paf to xy and axis metadata
    Assumes pairs of points representing HiC contact pairs, or Mashmap sequence similarity hits. HiC data typically comes from pairs of haplotypes and is used to help assemble all
    the contigs into chromosomes

    These coordinate pairs are of 3 types - both ends on one or the other reference sequence, or one end on each for HiC pairs.
    Let's call these cis if the same reference and trans if between the two references. The ones in cis are likely to be consistent-ish between the two haplotypes, but the trans ones turn
    out to be very interesting...

    For HiC data, there should really be little difference - whether a haplotype contacts another haplotype or itself depends on the 3D folding and since the haplotypes are wound together into
    a helix, then the whole total contig length - about 2 meters for mammals are folded up into a 10μ^3 ball.

    python holoSeq_prepare_gz.py --inFile mUroPar1H1H2.paf --xclenfile mUroPar1H1suffix.len --yclenfile mUroPar1H2suffix.len --contig_sort VGPname --hap_indicator Suffix
    """

    def __init__(self, inFname, args, xcontigs, ycontigs, haps, xwidth, ywidth):
        """
        if rotating paf line at a time, better to pre-calculate the constants once rather than once for each point
        """
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
        self.prepPafGZ(hsId, haps, xcontigs, ycontigs, args)
        if self.isGzip(inFname):
            with gzip.open(inFname, "rt") as f:
                self.readPAF(f)
        else:
            with open(inFname) as f:
                self.readPAF(f)
        self.cis1f.close()
        self.cis2f.close()
        self.transf.close()

    def readPAF(self, f):
        # might be a gz
        ncis1 = ncis2 = ntrans = 0
        for rowi, rows in enumerate(f):
            row = rows.strip().split()
            if len(row) > 7:
                c1 = row[0]
                c2 = row[5]
                n1 = int(row[2])
                n2 = int(row[7])
                H1 = data.getHap(c1)
                H2 = data.getHap(c2)
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
                            ntrans += 1
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
                            ntrans += 1
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
                            ncis1 += 1
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
                            ncis2 += 1
        log.debug("ncis1=%d, ncis2=%d, ntrans=%d" % (ncis1, ncis2, ntrans))

    def isGzip(self, inFname):
        with gzip.open(inFname, "r") as fh:
            try:
                fh.read(1)
                return True
            except gzip.BadGzipFile:
                log.info("inFname %s is not a gzip so will read as text" % inFname)
                return False

    def prepPafGZ(self, hsId, haps, xcontigs, ycontigs, args):
        """
        @v1HoloSeq2D for example
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
                "@%s %s %d" % (data.getHap(k), k, xcontigs[k]) for k in xcontigs.keys()
            ]
            if len(haps) > 1:
                h += [
                    "@%s %s %d" % (data.getHap(k), k, ycontigs[k])
                    for k in ycontigs.keys()
                ]
            metah = [
                hsId,
                "@@heatmap",
                "@@title %s" % args.title + subtitle,
                "@@datasource %s" % "paf",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
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
        f1 = gzip.open(fn1, mode="wb")
        self.cis1f = io.BufferedWriter(f1, buffer_size=1024 * 1024)
        prepHeader(
            haps[0],
            hsId,
            xcontigs,
            xcontigs,
            args,
            self.cis1f,
            " Pairs on %s" % haps[0],
            xclenfile=args.xclenfile,
            yclenfile=args.xclenfile,
            ax=haps[0],
        )
        fn2 = "%s_cis%s_hseq.gz" % (self.inFname, haps[1])
        if self.rotate:
            fn2 = "%s_rotated_cis%s_hseq.gz" % (self.inFname, haps[1])
        f2 = gzip.open(fn2, mode="wb")
        self.cis2f = io.BufferedWriter(f2, buffer_size=1024 * 1024)
        prepHeader(
            haps[1],
            hsId,
            ycontigs,
            ycontigs,
            args,
            self.cis2f,
            " Pairs on %s" % haps[1],
            xclenfile=args.yclenfile,
            yclenfile=args.yclenfile,
            ax=haps[1],
        )
        fn3 = "%s_trans_hseq.gz" % (self.inFname)
        if self.rotate:
            fn3 = "%s_rotated_trans_hseq.gz" % (self.inFname)
        f3 = gzip.open(fn3, mode="wb")
        self.transf = io.BufferedWriter(f3, buffer_size=1024 * 1024)
        prepHeader(
            haps,
            hsId,
            xcontigs,
            ycontigs,
            args,
            self.transf,
            " Pairs on different haplotypes",
            xclenfile=args.xclenfile,
            yclenfile=args.yclenfile,
            ax="BOTH",
        )
