
## HoloSeq axis contig names and lengths: Avoiding coordinate slop.

One of HoloSeq’s design goals is to allow multiple feature tracks to be reliably compared visually. That implies tracks being compared are all derived from, and displayed, in a consistent axis coordinate system. Preventing users from creating misaligned displays is an important aspect of the application architecture. For each species, the historical record of genomes reflect improvements, as the reference sequence is progressively refined from the initial assembly. Without explicit metadata, there is a substantial risk that users will try to compare feature tracks that are not directly comparable, leading to misleading visualisations.

For some tracks, such as hic and bigwig, if the generating software faithfully stores all contigs and their lengths from the reference fasta as recoverable metadata, there is no problem - so Charles’ hic data probably has correct internally consistent metadata for example. For a single track that will only ever be viewed by itself, deriving the contigs and their lengths to create the linear axes, from the feature data itself, is one option. 

For paf, bed, gff3 or vcf, there is no guarantee that every contig will be represented in the data and no explicit contig metadata is stored in those formats. If one small contig happens not to have any features, it will be missing from the inferred axis. The axes rely on cumulated contig lengths, so  "slop" may be introduced in terms of the relationship between axis position and precise genomic location, if trying to derive the axis from the contigs found in a paf. The viewer will assume that stacked tracks are all on the same backbone assembly, but features will be misaligned because of cumulated coordinate slop, or worse, from different axis assembly versions or even species.

An unambiguous source of truth is needed - unless the tracks will only ever be viewed in complete isolation, which limits their utility. The HoloSeq display code explicitly allows multiple processed coordinate gzips to be stacked into a final display - a paf heatmap can be shown together with a 1D coverage bigwig for example. Fidelity of that combined display will depend on having uniform and accurate metadata for the axes. The best source of truth is an explicit contig lengths file, so each individual track has metadata to help prevent mix-ups - for example to ensure that every feature track being compared was mapped with the same hg38 reference fasta.  A suitable lengths file can be made using samtools faidx on the original fasta used for feature mapping. Only the first two columns are needed, containing the contig names and lengths.

Updated holoSeq code will now require a chromosome lengths file from the same genome fasta used to generate each track input to the converter. This will avoid any risk of slop, or in the case of hic, mismatches between assemblies used to derive feature tracks. The effect of coordinate slop will usually be small, but the risk of presenting locations of features on two different assemblies as if they were the same, should be minimised by good design. Each lengths file describes a consistent axis for simultaneous display of any arbitrary group of pre-computed coordinates and their metadata.

The downside for users is that the input track data (gff, bigwig, vcf, paf) will need to be accompanied by a contig lengths file derived from the same fasta that was used to generate the track.

The benefit is that the contig lengths file name is stored in the pre-computed gzip as part of the metadata, so when a user combines multiple hseq.gz files as inputs for a display, the code can determine which feature tracks can be directly compared on the same axis. Where the axis metadata matches, the two tracks share the same fasta “backbone”, so they can reliably be stacked and visually compared, and axes can be linked without slop.