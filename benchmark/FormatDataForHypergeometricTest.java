import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

public class FormatDataForHypergeometricTest {

	public static void main(String[] args) {

		String annotationFile = "/Users/rnadeau2/Documents/Structures/benchmark/corrNetTop2-400_structureModules_4.0.tsv";
		String mclFile = "/Users/rnadeau2/Documents/Structures/benchmark/mclOutput_i6.txt";
		String outputPrefix = "/Users/rnadeau2/Documents/Structures/benchmark/i6/moduleInfo_";
		
		formatHypergeometricInput(mclFile, annotationFile, outputPrefix);
	}
	
	public static void formatHypergeometricInput(String mclFile, String annotationFile, String outputPrefix) {
		
		/* map proteins to MCL clusters */
		HashMap<String, Integer> mclClusterMap = loadMCLclusters(mclFile);

		/* read annotation file - generate 1 file for HyperGeometric test per module */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(annotationFile)));

			String line = in.readLine(); // header 
			line = in.readLine();
			int module = 0;
			
			while(line!=null) {
				
				/* get proteins in cluster */
				HashSet<String> proteins = new HashSet<>();
				if(line.split("\t").length >= 3) {
					proteins.addAll(Arrays.asList(line.split("\t")[2].split("\\|")));
				}
				
				/* print module info */
				printDataForHypergeometricTest(outputPrefix, module, proteins, mclClusterMap);
				
				line = in.readLine();
				module++;
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Load MCL clusters as a map of {protein = cluster number}
	 * @param mclFile
	 * @return 
	 */
	public static HashMap<String, Integer> loadMCLclusters(String mclFile){
		
		HashMap<String, Integer> mclClustersMap = new HashMap<>(); // key = protein, value = cluster number
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(mclFile)));

			String line = in.readLine();
			int cluster = 1;
			
			while(line!=null) {
			
				String[] prots = line.split("\t");
				for(String p : prots) {
					mclClustersMap.put(p, cluster);
				}
				line = in.readLine();
				cluster++;
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return mclClustersMap;
	}
	
	public static void printDataForHypergeometricTest(String outputFile, int cluster, HashSet<String> proteinsInClusterSet, HashMap<String, Integer> clusterMap) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + cluster + ".tsv")));

			out.write("Protein\tGroup\tCluster\n");
			for(Entry<String, Integer> e : clusterMap.entrySet()) {
				
				
				String prot = e.getKey();
				out.write(prot + "\t");
				
				// check module
				if(proteinsInClusterSet.contains(prot)) {
					out.write("Module\t");
				} else {
					out.write("Not\t");
				}
				
				// write cluster
				out.write(e.getValue() + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
