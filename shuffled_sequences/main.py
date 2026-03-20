# Shuffle all sequences in network using altschulEriksonDinuclShuffle.py

import altschulEriksonDinuclShuffle as diShuffl


def formatFasta(inputFasta, outputFasta):
    f = open(inputFasta, "r")
    out = open(outputFasta + ".fasta", "w")

    count = 0
    for line in f:

        if count % 2 == 0:
            # first line >RefSeqID - print  >ID_Shuffled
            out.write(line.rstrip() + "_Shuffled\n")
            print(line)
        else:
            # second line = sequence - print shuffled (call external script)
            seq = line.rstrip()
            print("sequence: ", len(seq))
            print("CDS: ", len(seq[0:100]))
            print("UTR: ", len(seq[100:]))
            out.write(diShuffl.dinuclShuffle(seq[0:100]) + diShuffl.dinuclShuffle(seq[100:]) + "\n")

        count += 1
        print(count)

    f.close()
    out.close()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Path to fasta with sequences
    wd = "/Users/rnadeau2/Documents/Structures/hcm/"
    inputFasta = wd + "corrNet2-400_3utr_w100cds.fasta"

    for i in range(11,16):
        outputFasta = wd + "corrNet2-400_3utr_w100cds_diNuclShuffled_rand" + str(i)
        formatFasta(inputFasta, outputFasta)
