package investigateModules;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class InvestigateGraphPropertiesOfModules {

	public static void main(String[] args) throws IOException {

		String wd = "/Users/rnadeau2/Documents/Structures/";
		String condition = args[0];
		String threshold = args[1];

		/* network specific files */
		String networkFile = wd + "enrichmentAnalysis/inputFiles/corrNetTop2-400_formattedNetwork.tsv";
		String proteinFile = wd + "enrichmentAnalysis/ioFiles/corrNetTop2-400_proteinsInNetwork_info.tsv";
		String connectedDMFile = wd + "enrichmentAnalysis/ioFiles/corrNetTop2-400_connectedComponent_distanceMatrix.tsv";

		/* annotation file */
		String annotationFile = wd + "Annotations/modules/annotations/corrNetTop2-400_structureModules_" + condition + "_" + threshold + ".tsv";
		String annotationPercentileFile = wd + "Annotations/modules/annotations/corrNetTop2-400_percentileAnnotations_r0.88_"+ condition +"_" + threshold + ".tsv";

		/* make directories */ 
		Files.createDirectories(Paths.get(wd + "Annotations/modules/" + condition + "/"));

		/* determine protein info: index in matrix & #degrees - Map<String = Protein> */
		System.out.println("loading proteins");
		HashMap<String, Protein> proteinMap = obtainProteinsInNetwork(proteinFile);

		/* load distance matrix as double[][] */
		System.out.println("loading dm");
		double[][] dm = loadDistanceMatrix(connectedDMFile, proteinMap.size());

		/* load edges as Map<Prot1_Prot2 = weight> */
		System.out.println("loading edge list");
		HashMap<String, Double> edgeList = obtainEdgeList(networkFile, proteinMap);

		/* Search annotations */
		try {
			BufferedReader in1 = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			BufferedReader in2 = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationPercentileFile))));

			/* headers */
			String line1 = in1.readLine();
			String line2 = in2.readLine();

			line1 = in1.readLine();
			line2 = in2.readLine();

			while(line1!=null) {

				String[] col = line1.split("\t");
				String[] col2 = line2.split("\t");

				if(!col[0].equals(col2[0])) {
					while(!col[0].equals(col2[0])) {
						line1 = in1.readLine();
						col = line1.split("\t");
					} 
				}
				if(col.length > 2) {

					/* format annotation */
					System.out.println("Module: " + col[0]);
					Annotation module = new Annotation(line1.split("\t")[2].split("\\|"), line2.split("\t")[2].split("\\|"));

					if(module.getProteinList().size()>1 && module.getProteinList().size() < 2000) {
						/* output files */
						Files.createDirectories(Paths.get(wd +"Annotations/modules/bp4_" + condition + "/m" + col[0] + "/"));

						String degreeDistributionFile = wd +"Annotations/modules/bp4_"+ condition + "/m" + col[0] + "/degreeDistribution_l2000_m" + col[0] + ".tsv";
						String bpScoreDistributionFile = wd +"Annotations/modules/bp4_"+ condition + "/m" + col[0] + "/bpScoreDistribution_l2000_m" + col[0] + ".tsv";
						String percentileDistributionFile = wd +"Annotations/modules/bp4_"+ condition + "/m" + col[0] + "/percentileDistribution_l2000_m" + col[0] + ".tsv";

						String edgeWeightDistributionFile =  wd +"Annotations/modules/bp4_"+ condition + "/m" + col[0] + "/edgeWeightDistribution_l2000_m" + col[0] + ".tsv";
						String shortestPathDistributionFile = wd +"Annotations/modules/bp4_"+ condition + "/m" + col[0] + "/shPathDistribution_l2000_m" + col[0] + ".tsv";
						String nwShortestPathDistributionFile = wd +"Annotations/modules/bp4_"+ condition + "/m" + col[0] + "/nwShPathDistribution_l2000_r0.88_m" + col[0] + ".tsv";

						/* assess degree distribution */ 
						module.assessDegreeDistribution(proteinMap, degreeDistributionFile);

						/* assess score distribution */
						module.assessBPscoreDistribution(bpScoreDistributionFile);

						/* assess percentile distribution */
						module.assessPercentileScoreDistribution(percentileDistributionFile);

						/* for each pair of protein */
						List<String> proteins = module.getProteinList();
						List<Double> percentiles = module.getPercentileScores();

						List<Double> edges = new ArrayList<>();
						List<Double> shortestPaths = new ArrayList<>();
						List<Double> nwShortestPaths = new ArrayList<>();

						for(int i=0; i<proteins.size(); i++) {
							for(int j=i+1; j<proteins.size(); j++) {

								/* assess edge distribution */
								String interactors1 = proteins.get(i) + "_" + proteins.get(j);
								String interactors2 = proteins.get(j) + "_" + proteins.get(i);

								if(edgeList.containsKey(interactors1)) {
									edges.add(edgeList.get(interactors1));
								}

								if(edgeList.containsKey(interactors2)) {
									edges.add(edgeList.get(interactors2));
								}

								Protein prot1 = proteinMap.get(proteins.get(i));
								Protein prot2 = proteinMap.get(proteins.get(j));

								/* assess shortest path distribution */ 
								double shPath = dm[prot1.getIdx()][prot2.getIdx()];
								shortestPaths.add(shPath);

								/* assess node weighted shortest path distribution */ 
								double nodeAverage = (percentiles.get(i) + percentiles.get(j)) / (double) 2;
								nwShortestPaths.add(shPath / nodeAverage);
							}
						}
						assessEdgeDistribution(edges, edgeWeightDistributionFile);
						assessShortestPathDistribution(shortestPaths, shortestPathDistributionFile);
						assessShortestPathDistribution(nwShortestPaths, nwShortestPathDistributionFile);
					}
				}
				line1 = in1.readLine();
				line2 = in2.readLine();
			}
			in1.close();
			in2.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public static HashMap<String, Protein> obtainProteinsInNetwork(String proteinInfoFile){

		HashMap<String, Protein> proteinMap = new HashMap<>();

		/* obtain index of proteins in DM */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfoFile))));
			String line = in.readLine();
			int lineCount=0;

			while(line!=null) {
				String prot = line.split("\t")[0];
				proteinMap.put(prot, new Protein(prot, lineCount));

				line = in.readLine();
				lineCount++;
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinMap;
	}

	public static double[][] loadDistanceMatrix(String distance_matrixFile, int numProteins) {
		/* Import distance matrix from text file, it's dimensions are based on the size of the 
		 * proteinNetwork List initially used to build the distance Matrix file */

		double[][] distanceMatrix = new double[numProteins][numProteins]; // Initialize distance matrix

		try {

			InputStream in = new FileInputStream(new File(distance_matrixFile));				
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // read first line

			int x = 0; // initialize row counter

			while (line != null) {

				String[] col = line.split("\t"); // split columns
				int y = 0; // initialize column counter (resets at the end of every row)

				for (String str : col) { // str is the index to go through all element of col

					distanceMatrix[x][y] = Double.parseDouble(str);;// set the value (distance) at the appropriate coordinates
					y++; // adds one to the value of y (change column)
				}
				x++; // adds one to the vale of x (change row)
				line = input.readLine(); // read next line

			}
			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		return distanceMatrix;
	} // end import distance matrix

	public static HashMap<String, Double> obtainEdgeList(String networkFile, HashMap<String, Protein> proteinMap){

		HashMap<String, Double> edges = new HashMap<>();

		/* determine #degrees from network file */ 
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(networkFile))));
			String line = in.readLine();

			while(line!=null) {
				String[] col = line.split("\t");

				/* if proteins are in the network - update degree count */
				if(proteinMap.containsKey(col[0]) && proteinMap.containsKey(col[1])) {

					if(!edges.containsKey(col[0] + "_" + col[1]) && !edges.containsKey(col[1] + "_" + col[0])) {
						proteinMap.get(col[0]).increaseDegrees();
						proteinMap.get(col[1]).increaseDegrees();

						edges.put(col[0] + "_" + col[1], (1-Double.parseDouble(col[2])));
					}
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return edges;
	}

	public static void assessEdgeDistribution(List<Double> edges, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("EdgeWeights\tcount\tnormalizeCount\n");

			for(double i=0; i<0.6; i=i+0.1) {

				int count = 0;
				for(double s: edges) {
					if(s >= i && s < (i+0.1)) {
						count++;
					}
				}
				out.write(i + "\t" + count + "\t" + (count / (double) edges.size()) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void assessShortestPathDistribution(List<Double> shortestPaths, String outputFile) {

		double min = Math.floor(Collections.min(shortestPaths));
		double max = Math.ceil(Collections.max(shortestPaths));

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("EdgeWeights\tcount\tnormalizeCount\n");

			for(double i=min; i<max; i= i+0.5) {

				int count = 0;
				for(double s: shortestPaths) {
					if(s >= i && s < (i+0.5)) {
						count++;
					}
				}
				out.write(i + "\t" + count + "\t" + (count / (double) shortestPaths.size()) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
