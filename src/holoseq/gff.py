from bisect import bisect_left
import gzip
import html
import logging
import numpy as np
import os

import urllib.request

from config import VALID_HSEQ_FORMATS
import data

import holoviews as hv
from holoviews.operation.datashader import (
    rasterize,
)

from holoviews.operation.element import apply_when
import panel as pn

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_prepare")


class gff:
    """
            Only care about mRNA cds and stop codons initally. Turn into segments. Filter so only data in contigs is retained from input.
    SUPER_1 miniprot        mRNA    139006290       139072696       22660   -       .       ID=MP000006;Rank=1;Identity=0.9979;Positive=0.9984;Target=XP_026244093.1 1 4350
        fix positions as we go with a lookup contig -> cumulated offset
        can pop open https://www.ncbi.nlm.nih.gov/protein/XP_026244093.1
    """

    def __init__(self, gff, outFname, contigs, args):
        mrnaseen = {}
        segs = {}
        comment = "#"
        self.hsId = VALID_HSEQ_FORMATS[0]
        self.inFname = gff
        log.debug("contigs=%s" % str(contigs)[:1000])
        with open(gff) as g:
            for i, row in enumerate(g):
                if not row.startswith(comment):
                    (id, name, kind, startp, endp, score, strand, phase, text) = [
                        x.strip() for x in row.split()[:9]
                    ]
                    if not segs.get(id, None):
                        segs[id] = []
                    if kind.lower() in ["cds", "mrna"]:
                        anno = text.split(";")
                        tanno = [
                            x.strip()[7:]
                            for x in anno
                            if x.lower().startswith("target=")
                        ]
                        target = tanno[0]
                    startp = int(startp)
                    endp = int(endp)
                    offset = contigs.get(id, -1)
                    if args.addH1 and offset < 0:
                        offset = contigs.get(id + "H1", -1)
                        id = id + "H1"
                    if offset < 0:
                        log.warn(
                            "Ignored gff3 id %s missing from supplied xcontigs, in row %d %s of %s with addH1=%s"
                            % (id, i, row, gff, args.addH1)
                        )
                    else:
                        if kind.lower() == "mrna":
                            if target:
                                if mrnaseen.get(target, None):
                                    log.debug(
                                        "Seeing mrna target %s again at row %d"
                                        % (target, i)
                                    )
                                else:
                                    mrnaseen[target] = target
                                segs[id].append(
                                    (
                                        startp + offset,
                                        endp + offset,
                                        strand,
                                        score,
                                        target,
                                        "mrna",
                                    )
                                )
                            else:
                                log.warn("no target found in %s at row %d" % (text, i))
                        elif kind.lower() == "stop_codon":
                            segs[id].append((startp + offset, target, "stopc"))
                        elif kind.lower() == "cds":
                            segs[id].append(
                                (
                                    startp + offset,
                                    endp + offset,
                                    strand,
                                    score,
                                    target,
                                    "cds",
                                )
                            )

        self.export_mapping(outFname, contigs, segs, args)

    def export_mapping(self, outFname, contigs, segs, args):
        """
        for GFF
        @v1HoloSeq2D for example
        A default  Y value of 100 is set for each mRNA's extent, but often there are dozens of different named sequences in the databases that will overlap.
        Assuming the GFF is sorted, overlapping mRNA is identified as ending or starting in the previous extent and the y value is decremented to avoid overlap
        Y is reset if no overlap
        """

        def prepHeader(contigs, args):
            """
            holoSeq output format
            """
            h = ["@%s %s %d" % (data.getHap(k), k, contigs[k]) for k in contigs.keys()]
            metah = [
                self.hsId,
                "@@GFF 1",
                "@@title %s" % args.title,
                "@@datasource GFF",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % args.xclenfile,
            ]

            return metah + h

        def ranges_overlap(x1, x2, y1, y2):
            # buried in https://stackoverflow.com/questions/6821156/how-to-find-range-overlap-in-python
            if x1 == x2 or y1 == y2:
                return False
            return x1 <= y2 and y1 <= x2

        hdr = prepHeader(contigs, args)

        with gzip.open(outFname, mode="wb") as ofn:
            ofn.write(str.encode("\n".join(hdr) + "\n"))
            y = 100
            for con in contigs.keys():
                lastseg = ("", 0, 0)
                subs = segs.get(con, [])
                if len(subs) > 0:
                    subs.sort(key=lambda x: x[0])
                    for i, m in enumerate(subs):
                        kind = m[-1]
                        if kind == "mrna":
                            (startp, endp, strand, score, targ, _) = m
                            if ranges_overlap(lastseg[1], lastseg[2], startp, endp):
                                y -= 1
                            else:
                                y = 100
                                lastseg = (targ, startp, endp)
                            row = str.encode(
                                f"mrna {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                            )
                            ofn.write(row)
                        elif kind == "cds":
                            (startp, endp, strand, score, targ, _) = m
                            row = str.encode(
                                f"cds {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                            )
                            ofn.write(row)
                        elif kind == "stopc":
                            (startp, targ, _) = m
                            row = str.encode(f"stopc {targ} {con} {startp}\n")
                            ofn.write(row)

    def makeGFFPanel(self, inFile, pwidth):
        """
                prepare a complete panel for the final display
        https://www.ncbi.nlm.nih.gov/gene/?term=XP_026235740.1

        import urllib.request
        xpuri = 'https://www.ncbi.nlm.nih.gov/gene/?term=XP_026235740.1'
        req = urllib.request.Request(xpuri)
        with urllib.request.urlopen(req) as response:
           apage = response.read()
        escaped_html = html.escape(apage)

        # Create iframe embedding the escaped HTML and display it
        iframe_html = f'<iframe srcdoc="{escaped_html}" style="height:100%; width:100%" frameborder="0"></iframe>'

        # Display iframe in a Panel HTML pane
        pn.pane.HTML(iframe_html, height=350, sizing_mode="stretch_width")
        """

        def get_ncbi(target):

            xpuri = 'https://www.ncbi.nlm.nih.gov/gene/?term=%s' % target
            req = urllib.request.Request(xpuri)
            with urllib.request.urlopen(req) as response:
                apage = response.read()
            escaped_html = html.escape(apage)
            # Create iframe embedding the escaped HTML and display it
            iframe_html = f'<iframe srcdoc="{escaped_html}" style="height:100%; width:100%" frameborder="0"></iframe>'
            # Display iframe in a Panel HTML pane
            pn.pane.HTML(iframe_html, height=350, sizing_mode="stretch_width")


        def showX(x, y):
            if np.isnan(x):
                s = "Mouse click on image for location"
            else:
                i = bisect_left(h1starts, x)
                chrx = h1names[i - 1]
                offsx = x - h1starts[i - 1]
                s = "%s:%d" % (chrx, offsx)
                xi = bisect_left(segs[xcf], x)
                xtarget = segs["target"][xi]
                s += " x %s" % (xtarget)

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
            self.import_holoSeq_data(inFile)
        )
        xcf = os.path.splitext(metadata["xclenfile"][0])[0]
        segs = {
            xcf: [],
            "x2": [],
            "wy1": [],
            "y2": [],
            "target": [],
            "colour": [],
            "thickness": [],
            "alpha": [],
        }
        """
cds XP_026238700.1 1401967516 1401967635 100 100 - 204
mrna XP_026248570.1 SUPER_3H1 531341254 531595863 100 100 + 1102 -1
cds XP_026248570.1 531341254 531341334 100 100 + 134
        """
        mthick = 3
        cdthick = 50
        for i, rows in enumerate(gffdata):
            if rows[0].lower() == "mrna":
                (kind, targ, contig, startp, endp, y1, y2, strand, score) = rows[:10]
                startp = int(startp)
                endp = int(endp)
                y1 = int(y1)
                y2 = int(y2)
                colr = "blue"
                if strand == "-":
                    colr = "maroon"
                segs["target"].append(targ)
                segs[xcf].append(startp)
                segs["x2"].append(endp)
                segs["wy1"].append(y1)
                segs["y2"].append(y2)
                segs["colour"].append(colr)
                segs["thickness"].append(mthick)
                segs["alpha"].append(1.0)
                # segs["stopc"].append(stopc)
            elif (
                rows[0].lower() == "cds"
            ):  #  f"cds {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                (kind, targ, contig, startp, endp, y1, y2, strand, score) = rows[:10]
                startp = int(startp)
                endp = int(endp)
                y = int(y1)
                colr = "blue"
                if strand == "-":
                    colr = "maroon"
                segs["target"].append(targ)
                segs[xcf].append(startp)
                segs["x2"].append(endp)
                segs["wy1"].append(y)
                segs["y2"].append(y)
                segs["colour"].append(colr)
                segs["thickness"].append(cdthick)
                segs["alpha"].append(1.0)
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
            segs,
            [xcf, "wy1", "x2", "y2"],
            vdims=["target", "colour", "thickness", "alpha"],
        )

        gffp.opts(
            title="title",
            color="colour",
            line_width="thickness",
            alpha="alpha",
            width=pwidth,
            height=300,
            xticks=qtic1,
            xrotation=45,
            scalebar=True,
            scalebar_range="x",
            scalebar_location="bottom_left",
            scalebar_unit=("bp"),
            fontsize={"xticks": 8, "yticks": 10},
            show_grid=True,
            autorange="y",
            tools=[
                "xwheel_zoom",
                "box_zoom",
                "tap",
                "xpan",
                "reset",
            ],
            default_tools=[],
            active_tools=["xwheel_zoom", "tap", "xpan"],
            shared_axes=True,
        )
        apply_when(gffp, operation=rasterize, predicate=lambda x: len(x) > 5000)
        taps = hv.streams.Tap(source=gffp, x=0, y=0)
        showloc = pn.bind(showX, x=taps.param.x, y=taps.param.y)
        gp = pn.pane.HoloViews(gffp)
        p = pn.Column(showloc, gp)

        return p, title


