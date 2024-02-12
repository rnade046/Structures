package info;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Properties;

public class SummarizeSignificantModules {

	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));	

		String wd = params.getProperty("working_directory");

		String networkType = params.getProperty("network_name");
		
		/* input file for module information */
		String moduleDetailsFile = wd + networkType + "_structure_clusteringDetails.tsv";
		String moduleMetaFile = wd + "reliable-info2.tsv";
		
		/* input files for protein information */
		String annotationFile = wd + params.getProperty("annotationFile"); 
		String proteinInfoFile = wd + "ioFiles/" + networkType + "_proteinsInNetwork_info.tsv";
		String jsonIdxFile = wd + args[1];
		
		double threshold = Double.parseDouble(args[2]); // 0.000902548470845385
		String summaryOutputFile = wd + "significantModules_summary.tsv";
		
		/* load significant modules (clustering details file) */
		HashMap<String, Module> modules = getSignificantModules(moduleDetailsFile, threshold);
 
		/* get module info */
		setModuleMetaData(moduleMetaFile, modules);
		
		/* load proteins (annotation file) */
		HashMap<String, Protein> proteins = loadProteins(annotationFile, modules);
	
		/* load refSeq IDs (protein info file) */
		setRefSeqIDs(proteinInfoFile, proteins);
		
		/* load JSON paths (JSON index file) */
		setJSONpaths(jsonIdxFile, proteins);
		
		/* print combined info */
		printModuleSummary(summaryOutputFile, modules.values(), proteins);
	}
	
	public static HashMap<String, Module> getSignificantModules(String inputFile, double threshold){

		HashMap<String, Module> modules = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine(); // header
			
			while(line!=null) {
				String[] col = line.split("\t");

				if(Double.parseDouble(col[3]) <= threshold) {
					/* constructor = Module(String n, int nProts, double p) */
					modules.put(col[0], new Module(col[0], Integer.parseInt(col[1]), Double.parseDouble(col[3])));
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return modules;
	}
	
	public static void setModuleMetaData(String moduleMetaDataFile, HashMap<String, Module> modules) {

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(moduleMetaDataFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line!=null) {
				String[] col = line.split("\t"); // [0] = module, [1] = atlasID
				
				if(modules.containsKey(col[0])) {
					modules.get(col[0]).setMetaData(col[1], col[1].split("\\_")[0]);
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static HashMap<String, Protein> loadProteins(String annotationFile, HashMap<String, Module> modules) {
		
		HashMap<String, Protein> proteinMap = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));

			String line = in.readLine();
			line = in.readLine();

			while(line!=null) {
				String[] col = line.split("\t");
				if(modules.containsKey(col[0])){
					
					String[] proteins = col[2].split("\\|");
					
					/* initialize protein in list */
					for(String p: proteins) {
						String name = p.split("\\_")[0];
						
						if(!proteinMap.containsKey(name)) {
							proteinMap.put(name, new Protein(name));
						}
					}
					/* set protein indexes in module */
					modules.get(col[0]).setProteins(proteins);
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinMap;
	}
	
	public static void setRefSeqIDs(String proteinInfoFile, HashMap<String, Protein> proteinMap) {
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfoFile))));

			String line = in.readLine();
			line = in.readLine();

			while(line!=null) {
				String[] col = line.split("\t");
				
				if(col.length > 1) { // i.e. contains at least (1) RefSeqID = [1]
					if(proteinMap.containsKey(col[0])) { // [0] = protein
						proteinMap.get(col[0]).setIds(col[1].split("\\|"));
					}
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void setJSONpaths(String jsonIdxFile, HashMap<String, Protein> proteinMap) {
		
		/* load paths */
		HashMap<String, String> idMap = new HashMap<>(); // RefSeqId = JSONpath
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(jsonIdxFile))));

			String line = in.readLine();
			line = in.readLine();

			while(line!=null) {
				String[] col = line.split("\t"); // [0] = refSeqId, [1] = jsonPath
				idMap.put(col[0], col[1]);
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* associate paths to ids in ProteinMap */
		for(Protein p : proteinMap.values()) {
			for(String id: p.getRefSeqIds()) {
				if(idMap.containsKey(id)) {
					p.setJSONPath(id, idMap.get(id));
				}
			}
		}
	}
	
	public static void printModuleSummary(String outputFile, Collection<Module> modules, HashMap<String, Protein> proteinMap) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("Module\tType\tAtlasId\tPval\t#Prots\tProteinList\n");
			for(Module m: modules) {
				
				/* module meta data */
				for(String info : m.getModuleSummary()) {
					out.write(info + "\t");
				}
				
				/* protein names */
				for(String protein : m.getProteins()) {
					Protein p = proteinMap.get(protein);
					out.write(p.getName() + "_[");
					
					for(Entry<String, String> e : p.getIds()) {
						out.write(e.getKey() + "=" + e.getValue() + ",");
					}
					out.write("]|");
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
