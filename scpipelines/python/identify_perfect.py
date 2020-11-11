import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("extract_cells_from_loom.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--whitelist", default=None, type=str,
                    help='a file barcodes extracted using umi whitelist')
parser.add_argument("--read1", default=None, type=str,
                    help='read1 fastq  file')
parser.add_argument("--read2", default=None, type=str,
                    help='read2 fastq file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')
args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Import whitelist    ############################## #
# ########################################################################### #


def perfect_barcode_umi(string):
    i = 0
    n = 0
    umis = []
    for x in range(0, len(string)):
        

        substr = string[x:x+2]
        if i % 2:
            pass
        else:
            if ("CC" in substr or "GG" in substr or "AA" in substr or "TT" in substr):
                n+=1
                if n == len(string)/2:
                    umi = string[::2]
                    return(umi)
            else:
                break
        i +=1



read1 = iotools.open_file(paste0(args.outname,".fastq.1.gz"),"w")
read2 = iotools.open_file(paste0(args.outname,".fastq.2.gz"),"w")



perfect = 0
total = 0
with pysam.FastxFile(args.read1) as fh, pysam.FastxFile(args.read2) as fh2:
    for record_fh, record_fh2  in zip(fh, fh2):

        total += 1
        barcode_umi = perfect_barcode_umi(record_fh.sequence)
        if barcode_umi:
            perfect += 1
            read1.write("@%s\n%s\n+\n%s\n" % (record_fh.name, barcode_umi, record_fh.quality[20:]))
            read2.write("@%s\n%s\n+\n%s\n" % (record_fh2.name, record_fh2.sequence, record_fh2.quality))
        else:
            pass

L.info("Number of perfect read pairs:")
L.info(perfect)
L.info("Total number of read pairs:")
L.info(total)

outf.close()
outf2.close()
