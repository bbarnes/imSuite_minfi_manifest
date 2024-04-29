
import simplesam

barcodes = {}
with open('read_id_barcode_umi.txt') as barcodes_file:
  for line in barcodes_file:
    # should check the delimiter in this file. If it's ' ' or \t or ','
    read_id, umi, barcode = line.rstrip().split()
    barcode[read_id] = (umi, barcode)
    # reading this entire file could use a TON of memory if
    # if you have lots of reads

# set the tag names - take a look at SAM spec to pick an appropriate one
barcode_tag = 'ZB'
umi_tag = 'ZU'

with simplesam.Reader(open('in.bam')) as in_bam:
  with simplesam.Writer(open('out.sam', 'w'), in_bam.header) as out_sam:
    for read in in_bam:
      read[umi_tag] = barcodes[read.qname][0]
      read[barcode_tag] = barcodes[read.qname][1]
      out_sam.write(read)

## End of file
