
DAT_DIR="/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/prbDB"

OUT_DIR="${DAT_DIR}/merged"
OUT_TSV="${OUT_DIR}/GRCh36-GRCh37-GRCh38-GRCm10.u50-cgn.seq-sorted.tsv.gz"

mkdir -p ${OUT_DIR}

gzip -dc ${DAT_DIR}/*u50*.tsv.gz | \
    # head -n 1000 | \
    cut -f 1,2,4,5 | \
    sort -u | \
    gzip -c -> ${OUT_TSV}

# End of file
