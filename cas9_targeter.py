#!/usr/bin/python

##################################################
## Daniel Pass | github.com/passdan | June 2016 ##
##  ~Following specifications of IGEM_Cardiff~  ##
##################################################

import argparse
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import re

description = """---------------
This cas9_targeter.py script takes a region or genome, deconstructs it and tests for paired target regions that satisfy the guidelines below. This is then passed to BLASTn to test for simple alignment against the reference dataset which could be the rest of this species' genome, or multiple cross-reactive species. Final checks should be perforemd using additional online tools for certainty. The tool was designed specifically for targetting species not catered for with online tools.

Reqired input is a fasta formatted file of the region you want to find targets within. Can be multiple fasta lines in one file (multiple regions).
Outputs are a fasta table of the generated probes, and a table.txt file of the same information in a graphical representation.

Following the guidelines from: http://www.clontech.com/GB/Products/Genome_Editing/CRISPR_Cas9/Resources/Designing_sgRNA
    Use the following guidelines to select a genomic DNA region that corresponds to the crRNA sequence of the sgRNA:
        - The 3' end of the DNA target sequence must have a proto-spacer adjacent motif (PAM) sequence (5'-NGG-3').
        - The 20 nucleotides upstream of the PAM sequence will be your targeting sequence (crRNA) and Cas9 nuclease will cleave approximately 3 bases upstream of the PAM.
        - The PAM sequence itself is absolutely required for cleavage, but it is NOT part of the sgRNA sequence and therefore should not be included in the sgRNA.
        - The target sequence can be on either DNA strand.
    Tips for designing sgRNAs: Through our own experience, we have identified additional tips for designing sgRNAs.
    We have found that the best sgRNAs for several tested genes have a G at position 1 and an A or T at position 17.
---------------"""

if __name__ == '__main__':
    # Define the input arguments
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', help='Input fasta format file')
    parser.add_argument('-o', default="cas9-output", help='Output handle')
    parser.add_argument('-B', action='store_true', help='Pass argument to Blast the results')
    parser.add_argument('-d', default='ecoli-k12-genome', help='Database to Blast the results against')

    # Read the arguments in
    args = parser.parse_args()

def main():
    # Generate output files
    outfasta = open(args.o + ".fasta", 'w')
    outtable = open(args.o + "-table.txt", 'w')

    # Run the main target finding function
    casTargeter(outfasta,outtable)

    # If you passed the -B flag to the script, it'll do a blast now of the data.
    # Got to make sure blastn is installed and findable (in the PATH) on your machine or it will fail, maybe with unhelpful errors, maybe silently...
    if args.B is True:
        blastFunc(args.o + ".fasta")

    # All done!
    print "All done, see files:", args.o + "-table.txt,", args.o + ".fasta and", args.o + ".bln"


def casTargeter(outfasta,outtable):
    # Load input
    print "Reading input fasta file:", args.i
    for seq_record in SeqIO.parse(args.i, "fasta"):
        i = 0
        # Create empty dictionaries
        targets =[]
        targetDict = {}
        fwdDict = {}
        revDict = {}

        # Deconstructs the input into 23bp lengths (20bp target and the 3bp PAM sequence) and puts it in a dictionary
        print "Building 20bp dictionary."
        while i < len(seq_record):
            targets.append(seq_record.seq[i:(i+23)])
            targetDict[i] = seq_record.seq[i:(i+23)]
            i += 1

        # Making the graphical table output
        print "Testing if regions match cas9 design guidlines."
        outtable.write(' '.join(("Gene:", seq_record.description, "\n")))
        outtable.write("Pos  : 5' - 20bp sgRNA      | PAM | Spacer | PAM | RevComp 5'-20bp sgRNA\n")
        outtable.write("========================================================================\n")

        # Finding forward strand targets
        for pos in sorted(targetDict.keys()):
            t = targetDict[pos]
            # The important part: This regex says that the 20bp must | Start with a G | have 15bp or anything | A or a T | 3bp of anything | any character then GG |
            if re.match('^G.{15}[A|T].{3}.GG$', str(t)) is not None:
                # Make sure sequences aren't full of N's
                if re.search('N', str(t)) is not None:
                    print '[NB: ambiguous bases present]'
                else:
                    # Making a dictionary of the table outputs
                    a = t[0:20], "|", t[20:23]
                    fwdDict[pos] = ' '.join(str(i) for i in a)

        # The same thing for the reverse
        for pos in sorted(targetDict.keys()):
            t = targetDict[pos].complement()
            # Reversed version of the regex above
            if re.match('^GG..{3}[A|T].{15}G', str(t)) is not None:
                if re.search('N', str(t)) is not None:
                    print '[NB: ambiguous bases present]'
                else:
                    a = t[0:3], "|", t[3:23]
                    revDict[pos] = ' '.join(str(i) for i in a)

        matchInt = 0

        # Outputting the data but only if the forward and reverse pairs are within a certain range
        for fwdpos in sorted(fwdDict.keys()):
            for revpos in sorted(revDict.keys()):
                spacer = (revpos + 3) - (fwdpos + 20)
                # If the distance between the pairs is between 5bp and 100bp
                if (spacer < 100) and (spacer > 5):
                    # Output the table
                    outtable.write(' '.join((str(fwdpos).zfill(4),":", fwdDict[fwdpos], "--", str(spacer), "bp--", revDict[revpos], "\n")))
                    matchInt += 1
                    # output the fasta file
                    outfasta.write(''.join((">", str(i + matchInt) + "-forward | FwdPos:", str(fwdpos).zfill(4), "\n")))
                    outfasta.write(''.join((fwdDict[fwdpos][0:20], "\n")))
                    outfasta.write(''.join((">", str(i + matchInt) + "-reverse | FwdPos:", str(fwdpos).zfill(4), "\n")))
                    outfasta.write(''.join((revDict[revpos][6:26], "\n")))

    # Close the outputs so they can be read by the blast program
    outfasta.close()
    outtable.close()

def blastFunc(blastfasta):
    print "Running blastn for all generated sequences."
    blastn_cline = NcbiblastnCommandline(query=blastfasta, db=args.d, outfmt=7, out=args.o + ".bln")()
    #stdout, stderr = blastn_cline()
    #print blastn_cline

if __name__ == "__main__":
    main()
