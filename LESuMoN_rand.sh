#!/bin/bash 

# dependencies: java, viennaRNA package, python 10+, RNAbayespairing(+dependencies)

# input: fasta for CDS and 3'UTR
cds_fasta=$0
utr_fasta=$1
mapping_file=$2
output="fasta/seq_"

java sequences/FormatSequencesForRNAfold "struct/${output}" $struct

# RNA fold steps
# create a fasta file per sequence in text file
split -l 2 -d -a 5 --additional-suffix=".fasta" $input $output

nFiles=$(ls ${output}*.fasta | wc -l)

# fold with RNAfold
echo "fold RNA"
RNAfold -o${output} --noPS ${input}
mv ${output} "struct/"

## format structure file for BayesPairing
echo "format RNAfold output"
java sequences/FormatRNAfoldOutput "struct/${output}" $struct

# execute bayespairing
echo "Running BP"
for (( i=0; i<$nFiles; i++ ))
do
        echo $i
        if [ $i -lt 10 ]
        then
        inputI="${output}0${i}.fasta"
        else
        inputI="${output}${i}.fasta"
        fi

        outputI="${output}${i}"

        python3 /bayespairing/src/parse_sequences.py -seq $inputI -ss infile -d RELIABLE -o $outputI
done

# remove temporary files 
rm ${output}*.pickle

# evaluate BP output and generate annotation file
jsonIdx="json_idx.tsv"
assessBayesPairing/IndexJSONFromFasta $output $outputI $jsonIdx
assessBayesPairing/GenerateAnnotationFile $wd $threshol $condition $shuffled
assessBayesPairing/FilterAnnotationFile $annotationFile $filteredAnnotations

#local enrichment in network
localEnrichment/Main $params
