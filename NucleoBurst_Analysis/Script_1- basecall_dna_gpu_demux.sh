
# Info: 
# https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revt_14dec2018/barcoding-demultiplexing

guppy_basecaller -r -i ./fast5 -s ./fastq -c dna_r9.4.1_450bps_sup.cfg -x "cuda:0" --barcode_kits EXP-PBC096
