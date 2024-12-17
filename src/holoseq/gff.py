import gzip
import io
import logging

from config import VALID_HSEQ_FORMATS
import data



logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_prepare")



class gffConvert:
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
        self.hsId = config.VALID_HSEQ_FORMATS[0]
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
            h = ["@%s %s %d" % (getHap(k), k, contigs[k]) for k in contigs.keys()]
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
