import os
from itertools import combinations
import subprocess
import statistics


def combine_data(inputPrefix, outputFile, fdr_output, comb, rescale, significantMotifsFile, tpd):

    for i in comb:
        # Create the file
        combinedResultsFile = f"{outputFile}{'_'.join(map(str, i))}.tsv"
        with open(combinedResultsFile, 'w') as out:
            for fileIdx in i:

                combinedFile = f"{inputPrefix}{fileIdx}_struct_bp4_significantScores.tsv"
                if not tpd:
                    combinedFile = f"{inputPrefix}{fileIdx}_struct_bp4_r{rescale}_significantScores.tsv"

                with open(combinedFile, 'r') as infile:
                    out.write(infile.read())

        # Execute FDR check
        params = "/Users/rnadeau2/Documents/Structures/enrichmentAnalysis/LOCAL_motifs_params_file.txt"
        jScript = "/Users/rnadeau2/Documents/GitHub/Structures/fdr/ MainFDR"
        oFile = f"{fdr_output}{'_'.join(map(str, i))}.tsv"
        cmd = f"java -cp {jScript} {params} {combinedResultsFile} {oFile} {significantMotifsFile}"
        print(cmd)
        subprocess.run(cmd, shell=True, capture_output=True, text=True)


def summarize_fdrs(fdr_input_file_prefix, summary_file, comb):
    with open(summary_file, 'w') as out:
        s = '\t'.join(f'S{x}' for x in range(1, len(samplings[0]) + 1))
        out.write(s + "\tFDR\tPval\t#Motifs\t#ShuflledMotifs\n")

        for i in comb:
            if os.path.exists(f"{fdr_input_file_prefix}{'_'.join(map(str, i))}.tsv"):
                with open(f"{fdr_input_file_prefix}{'_'.join(map(str, i))}.tsv", 'r') as inFile:
                    next(inFile)  # skip header
                    for line in inFile:
                        l = line.split("\t")
                        if float(l[0]) > 0 and int(l[2]) > 0:
                            out.write('\t'.join(map(str, i)) + "\t" + line)
                            break


def evaluateFDRCoefficientOfVariation(summaryFile, outputFile, samples):

    # create header
    if not os.path.exists(outputFile):
        # Open the file in write mode and write the initial content
        with open(outputFile, 'w') as file:
            file.write("nSamples\tmean(fdrs)\tstDev(fdrs)\tcoeff(fdrs)\n")

    fdrs = []
    with open(summary_file, 'r') as inFile:
        next(inFile)  # skip header
        for line in inFile:
            fdrs.append(float(line.split("\t")[samples]))

    fdr_mean = statistics.mean(fdrs)
    fdr_stdev = statistics.stdev(fdrs)

    fdr_coeff = fdr_stdev / fdr_mean

    with open(outputFile, 'a') as out:
        out.write(f"{samples}\t{fdr_mean}\t{fdr_stdev}\t{fdr_coeff}\n")


if __name__ == '__main__':

    tpd = True
    #wd = "/Users/rnadeau2/Documents/Structures/enrichmentAnalysis/fdr/bp0/tpd/"
    #inputFilePrefix = f"{wd}corrNet2-400_tpd_rand"
    #signifcantMotifs = f"{wd}corrNet2-400_tpd_fwd_struct_bp0_significantScores.tsv"

    wd = "/Users/rnadeau2/Documents/Structures/enrichmentAnalysis/fdr/robust/corrNet2_400/nwTPD2/"
    inputFilePrefix = f"{wd}corrNet2-400_nwTPD2_rand"
    signifcantMotifs = f"{wd}corrNet2-400_nwTPD2_fwd_struct_bp4_significantScores.tsv"

    coefficient_file = f"{wd}/coefficient_summary.tsv"

    #rescales = [0, 0.5, 0.85]
    rescales = [0]

    for r in rescales:

        if not tpd:
            wd = f"/Users/rnadeau2/Documents/Structures/enrichmentAnalysis/fdr/upsampling/nwTPD2/r{r}/"
            inputFilePrefix = f"{wd}corrNet2-400_nwTPD2_rand"
            signifcantMotifs = f"{wd}corrNet2-400_nwTPD2_fwd_struct_bp4_r{r}_significantScores.tsv"
            coefficient_file = f"{wd}/coefficient_summary_r{r}.tsv"

        #elements = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        elements = [1, 2, 3, 4, 5]
        for i in range(3, len(elements) +1):
        #for i in 1:
        #i = 1
        # create directory
            newpath = f"{wd}/n{i}/"
            if not os.path.exists(newpath):
                os.makedirs(newpath)

            # Use itertools.combinations to generate combinations
            samplings = list(combinations(elements, i))

            #combine the upsamplings and perform FDR
            combined_scores_prefix = f"{wd}/n{i}/significantScores_n{i}_"
            fdr_output_prefix = f"{wd}/n{i}/fdr_"

            combine_data(inputFilePrefix, combined_scores_prefix, fdr_output_prefix, samplings, r, signifcantMotifs, tpd)

            # summarize FDRs
            summary_file = f"{wd}/n{i}/summary_FDR_{i}.tsv"
            summarize_fdrs(fdr_output_prefix, summary_file, samplings)

            # perform coefficient analysis
            if(i < len(elements)):
                evaluateFDRCoefficientOfVariation(summary_file, coefficient_file, i)
