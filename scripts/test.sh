$TESTDAT=testdata
python holoseq_prepare_gz.py --inFile $TESTDAT/mUroPar1H1H2.paf.gz  --xclenfile $TESTDAT/mUroPar1H1H2.len --yclenfile $TESTDAT/mUroPar1H1H2.len --contig_sort VGPname --hap_indicator Suffix --title "VGP mUroPar1 HiC data" --inFtype pair2d
