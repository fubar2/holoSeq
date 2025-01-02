TESTDAT="./tests"
#python scripts/holoseq_prepare_gz.py --inFile $TESTDAT/mUroPar1H1H2.paf.gz --inFtype pair2d --xclenfile $TESTDAT/mUroPar1H1H2.len --yclenfile $TESTDAT/mUroPar1H1H2.len --contig_sort VGPname --hap_indicator Suffix --title "VGP mUroPar1 HiC data" 

#python scripts/holoseq_prepare_gz.py --inFile $TESTDAT/mUroPar1_protH1.gff --inFtype gff --xclenfile $TESTDAT/mUroPar1H1suffix.len --yclenfile $TESTDAT/mUroPar1H1suffix.len --contig_sort VGPname --hap_indicator Suffix --title "VGP mUroPar1 NCBIprotein alignments H1" 

python scripts/holoseq_prepare_gz.py --inFile $TESTDAT/h1coverage.bw --inFtype bigwig --xclenfile $TESTDAT/mUroPar1HC.len --contig_sort VGPname --hap_indicator Suffix --title "VGP mUroPar1 H1 PacBio depth of coverage" 
