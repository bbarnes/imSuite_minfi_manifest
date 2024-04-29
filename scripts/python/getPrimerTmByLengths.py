
import primer3
import sys
import csv

manifest = sys.argv[1];

data = []
with open(manifest) as f:
    reader = csv.DictReader(f)
    data = [r for r in reader]

print("IlmnID\tTrimLength\tTrimmed_AlleleA_ProbeSeq\tOriginalTM\tPrimer3_TM")
for idx in range(len(data)):
    orgSeq = data[idx]['AlleleA_ProbeSeq']
    name   = data[idx]['IlmnID']
    orgTm  = data[idx]['Tm']
    for jj in range(0,25):
        subSeq = orgSeq[jj:]
        tm  = primer3.calcTm(subSeq)
        print("%s\t%d\t%s\t%s\t%s" % (name, jj, subSeq, orgTm, tm))




#lines = [line.rstrip('\n') for line in open('filename')]

#print primer3.calcTm('GTAAAACGACGGCCAGT')

## End of file
