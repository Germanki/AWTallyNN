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
L = logging.getLogger("complement_ployA.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='infile fastq  file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')

args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #


tab = str.maketrans("ACTG", "TGAC")

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]

forward_strand = "GACAGCCTTC" #must be caps ttcctgcaagaaatcaagccg
reverse_strand = "CTTGCAGGAA" #must be caps cggcttgaTTTCTTGCAGGAA

outfile = open(args.outname, "w")
log =  iotools.open_file(args.outname + ".log","w")
n = 0
y = 0
number_of_reads_below_150 = 0
with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        y +=1
        #Check length of read
        if len(record.sequence) < 150:
            number_of_reads_below_150 +=1
        # Check for forward strand
        if forward_strand in record.sequence:
            n +=1
            outfile.write("@%s\n%s\n+\n%s\n" % (record.name, record.sequence, record.quality))
        else:
            # Check for reverse strand and reverse complement if found
            if reverse_strand in record.sequence:
                n +=1
                sequence = reverse_complement_table(str(record.sequence))
                quality = str(record.quality)[::-1]
                outfile.write("@%s\n%s\n+\n%s\n" % (record.name, sequence, quality))

        

log.write("The number of total reads with polyA: %s\n" %(n))
log.write("The number of total reads is: %s\n" %(y))
log.write("The number of total recovered percent is: %s\n" %((n/y)*100))
log.write("The number of total reads with length < 150: %s\n" %(number_of_reads_below_150))

log.close()
outfile.close()
