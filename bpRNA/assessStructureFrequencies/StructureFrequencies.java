package assessStructureFrequencies;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class StructureFrequencies {

	public static void main(String[] args) {


		/* load structural elements */
		String out = "/Users/rnadeau2/Documents/Structures/bpRNAseqs/sequenceFrequencies_5testSeqs.tsv";

		List<Sequence> seqs = new ArrayList<>();


		for(int i=1; i<=5; i++) {

			String file = "/Users/rnadeau2/Documents/Structures/bpRNAseqs/seq"+ i +".st";
			seqs.add(new Sequence(file));
		}
		
		printFrequencies(seqs, out);
	}

	public static void printFrequencies(List<Sequence> seqs, String outputFile) {

		// get set of structures
		HashSet<String> structures = new HashSet<>();
		for(Sequence s: seqs) {
			structures.addAll(s.frequencyMap.keySet());
		}

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			// header 
			out.write("Structure\t");
			for(int i=0; i<seqs.size(); i++) {
				out.write("Seq" + (i+1) + "\t");
			}
			out.write("\n");

			// body
			for(String struct:  structures) {
				
				out.write(struct + "\t");
				
				for(int i=0; i<seqs.size(); i++) {
					Sequence s = seqs.get(i);
					if(s.frequencyMap.containsKey(struct)) {
						out.write(s.frequencyMap.get(struct) + "\t");
					} else {
						out.write("0\t");
					}
				}
				out.write("\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


	}

}
