package assessMembership;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

public class DetermineModulesMembership {

	public static void main(String[] args) {

		String annotationFile = "/Users/rnadeau2/Documents/Structures/Annotations/modules/annotations/corrNetTop2-400_structureModules_fwdSeq_4.0.tsv";

		String outputFile = "/Users/rnadeau2/Documents/Structures/Annotations/modules/membership/corrNet2-400_fwdSeq_bp4_proteinsByModule.tsv";
		String outputFile2 = "/Users/rnadeau2/Documents/Structures/Annotations/modules/membership/corrNet2-400_fwdSeq_bp4_proteinsByModule_wScores.tsv";

		/* Get list of significant motifs */
		HashMap<String, Module> modules = getSignificantModules(annotationFile);

		/* combine info */
		printCombinedInfo(outputFile, modules);
		printCombinedInfoWithScores(outputFile2, modules);
	}

	public static HashMap<String, Module> getSignificantModules(String inputFile){

		HashMap<String, Module> modules = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line!=null) {
				String[] col = line.split("\t");
				if(col.length>2) {
					modules.put(col[0], new Module(col[0], col[2].split("\\|")));
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return modules;
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
	
	public static void printCombinedInfoWithScores(String outputFile, HashMap<String, Module> modules) {

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
					
					/* get index map */
					if(m.containsProtein(prot)) {
						out.write(m.getProteinScores().get(m.getProteinIdxMap().get(prot)) + "\t");
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
