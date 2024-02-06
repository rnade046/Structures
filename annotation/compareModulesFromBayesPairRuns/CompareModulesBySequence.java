package compareModulesFromBayesPairRuns;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

public class CompareModulesBySequence {

	public static void main(String[] args) {

		String jsonMappingFile1 = "jsonIdxOfRefSeqIds_fwd.tsv";
		String path1 = "/home/rnade046/scratch/bp/";

		String jsonMappingFile2 = "corrNet2-400_nwTPD_fwd_bp4_idxOfFasta.tsv";
		String path2 = "/home/rnade046/scratch/bp_nwTPD/bp_corrNet2-400_3utr_w100cds_";

		String sequenceInfoOutputFile = "corrNet2-400_nwTPD_fwd_bp4_bpComp_seqInfo.tsv";
		String moduleInfoOuputFile = "corrNet2-400_nwTPD_fwd_bp4_bpComp_moduleInfo.tsv";

		/* get list of proteins to RefSeqIds from mapping2 file */
		HashMap<String, HashSet<String>> proteinIdMap = mapProteinsToRefSeqIds(jsonMappingFile2);

		/* get file names for original run */
		HashMap<String, String> idFileMapping1 = getFileNamesFromOriginalRun(jsonMappingFile1, path1);

		/* get file names for new run */
		HashMap<String, String> idFileMapping2 = getFileNamesFromSecondRun(jsonMappingFile2, path2);

		HashMap<String, Sequence> seqMap = compareSequenceModules(proteinIdMap, idFileMapping1, idFileMapping2);

		printSequenceInfo(sequenceInfoOutputFile, seqMap);
		
		HashMap<String, Module> moduleMap = assessModuleOverlap(proteinIdMap, seqMap);
		
		printModuleInfo(moduleInfoOuputFile, moduleMap);
		
	}

	public static HashMap<String, HashSet<String>> mapProteinsToRefSeqIds(String mappingFile){

		HashMap<String, HashSet<String>> idMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(mappingFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while (line != null) {
				/* idx	ID	Length	Protein */
				String[] col = line.split("\t");

				if(idMap.containsKey(col[3])) {
					idMap.get(col[3]).add(col[1]);
				}else {
					HashSet<String> ids =  new HashSet<>();
					ids.add(col[1]);
					idMap.put(col[3], ids);
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return idMap;
	}

	public static HashMap<String, String> getFileNamesFromOriginalRun(String jsonMapping, String path){

		HashMap<String, String> fileMapping = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(jsonMapping))));

			String line = in.readLine(); // header

			while (line != null) {
				String[] col = line.split("\t");
				fileMapping.put(col[0], path + col[1]);

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return fileMapping;
	}

	public static HashMap<String, String> getFileNamesFromSecondRun(String jsonMapping, String path){

		HashMap<String, String> fileMapping = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(jsonMapping))));

			String line = in.readLine(); // header

			while (line != null) {
				String[] col = line.split("\t");	/* idx	ID	Length	Protein */

				fileMapping.put(col[1], path + col[0] + ".json");

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return fileMapping;
	}

	public static HashMap<String, Sequence> compareSequenceModules(HashMap<String, HashSet<String>> proteinIdMap, HashMap<String, String> idFileMapping1, HashMap<String, String> idFileMapping2){

		HashMap<String, Sequence> seqMap = new HashMap<>();

		for(Entry<String, HashSet<String>> idInfo : proteinIdMap.entrySet()) {
			for(String id: idInfo.getValue()) {

				if(idFileMapping1.containsKey(id) && idFileMapping2.containsKey(id)) {
					System.out.println(id);
					seqMap.put(id, new Sequence(id, idInfo.getKey(), idFileMapping1.get(id), idFileMapping2.get(id)));
				}
			}
		}
		return seqMap;
	}

	public static void printSequenceInfo(String outputFile, HashMap<String, Sequence> sequenceMap) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("Protein\tID\t#ModulesKept\t#ModulesGained\t#ModulesLost\tListModuleKept\tListModuleGained\tListModuleLost\n");

			for(Entry<String, Sequence> seq : sequenceMap.entrySet()) {

				Sequence value = seq.getValue();

				out.write(value.getProteinName() + "\t" + seq.getKey() + "\t" + value.getKeptModules().size() + "\t" + value.getNewModules().size() + "\t" + value.getLostModules().size() + "\t");

				/* modules kept */
				for(String m : value.getKeptModules()) {
					out.write(m + "|");
				}
				out.write("\t");

				/* modules gained */
				for(String m : value.getNewModules()) {
					out.write(m + "|");
				}
				out.write("\t");

				/* modules lost */
				for(String m : value.getLostModules()) {
					out.write(m + "|");
				}
				out.write("\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static HashMap<String, Module> assessModuleOverlap(HashMap<String, HashSet<String>> proteinIdMap, HashMap<String, Sequence> seqMap){

		HashMap<String, Module> moduleMap = new HashMap<>();

		/* iterate for each protein */
		for(Entry<String, HashSet<String>> protEntry : proteinIdMap.entrySet()) {

			String protein = protEntry.getKey();
			System.out.println("Protein: " + protein);

			/* iterate for each sequence */
			for(String id: protEntry.getValue()) {

				Sequence seq = seqMap.get(id);
				System.out.println(id + "|");
				
				/* kept modules = kept proteins */
				for(String module : seq.getKeptModules()) {
					if(!moduleMap.containsKey(module)) {
						moduleMap.put(module, new Module(module));
					}
					moduleMap.get(module).getKeptProteins().add(protein);
				}

				/* new modules = new proteins */
				for(String module : seq.getNewModules()) {
					if(!moduleMap.containsKey(module)) {
						moduleMap.put(module, new Module(module));
					}
					moduleMap.get(module).getNewProteins().add(protein);
				}

				/* lost modules = lost proteins */
				for(String module : seq.getLostModules()) {
					if(!moduleMap.containsKey(module)) {
						moduleMap.put(module, new Module(module));
					}
					moduleMap.get(module).getLostProteins().add(protein);
				}
			}
			System.out.println();
		}

		return moduleMap;
	}

	public static void printModuleInfo(String outputFile, HashMap<String, Module> moduleMap) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("Module\t#ProteinsKept\t#ProteinsGained\t#ProteinsLost\tListProteinsKept\tListProteinsGained\tListProteinsLost\n");

			for(Module module : moduleMap.values()) {

				out.write(module.getModuleId() + "\t" + module.getKeptProteins().size() + "\t" + module.getNewProteins().size() + "\t" + module.getLostProteins().size() + "\t");

				/* proteins kept */
				for(String m : module.getKeptProteins()) {
					out.write(m + "|");
				}
				out.write("\t");

				/* proteins gained */
				for(String m : module.getNewProteins()) {
					out.write(m + "|");
				}
				out.write("\t");

				/* proteins lost */
				for(String m : module.getLostProteins()) {
					out.write(m + "|");
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
