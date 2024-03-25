package assessProteinsByModule;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

public class mapProteinsByModules {

	public static void main(String[] args) {

		String modulesDetailsFile = "/Users/rnadeau2/Documents/Structures/enrichmentAnalysis/fdr/bp4_seq/corrNet2-400_nwTPD2_rand1Seq_struct_bp4_clusteringDetails.tsv";
		String annotationFile = "/Users/rnadeau2/Documents/Structures/Annotations/modules/annotations/corrNetTop2-400_structureModules_seqRand1_4.0.tsv";

		String outputFile = "/Users/rnadeau2/Documents/Structures/Annotations/modules/membership/corrNet2-400_seqRand1_bp4_proteinsByModule.tsv";
		double threshold = 1.0;
		
		/* Get list of significant motifs */
		HashMap<String, Module> modules = getSignificantModules(modulesDetailsFile, threshold);

		/* Parse annotation file for protein information */
		setProteins(annotationFile, modules);

		/* combine info */
		printCombinedInfo(outputFile, modules);
	}

	public static HashMap<String, Module> getSignificantModules(String inputFile, double threshold){

		HashMap<String, Module> modules = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine(); // header
			
			while(line!=null) {
				String[] col = line.split("\t");

				if(Double.parseDouble(col[3]) <= threshold) {
					modules.put(col[0], new Module(col[0], Integer.parseInt(col[1])));
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return modules;
	}

	public static void setProteins(String annotationFile, HashMap<String, Module> modules) {

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));

			String line = in.readLine();
			line = in.readLine();

			while(line!=null) {
				String[] col = line.split("\t");
				if(modules.containsKey(col[0])){
					modules.get(col[0]).setProteins(col[2].split("\\|"));
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void printCombinedInfo(String outputFile, HashMap<String, Module> modules) {

		/* get list of all proteins */
		HashSet<String> proteins = new HashSet<>();
		for(Module m : modules.values()) {
			proteins.addAll(m.getProteins());
		}

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			// header
			out.write("Protein\t");
			for(Module m : modules.values()) {
				out.write(m.getProteinName() + "\t");
			}
			out.write("\n");

			for(String prot: proteins) {
				out.write(prot + "\t");

				for(Module m : modules.values()) {
					if(m.containsProtein(prot)) {
						out.write(1 + "\t");
					} else {
						out.write(0 + "\t");
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
