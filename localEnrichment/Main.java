import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Properties;

import graph.Interaction;
import graph.Protein;
import network.Calculator;
import network.CorrelationGraphLoader;
import network.DistanceMatrix;
import network.NetworkProteins;
import sampling.ApproximateNormalDistribuiton;
import sampling.MotifSampling;
import sampling.MotifSamplingPerModule;
import sampling.AssessProteinAnnotations;
import utils.AnnotationCompanionFile;
import utils.AssessEnrichment;
import utils.MotifEnrichment;
import utils.Scale;
import utils.ScorePercentiles;

public class Main {
	public static void main(String[] args) throws FileNotFoundException, IOException {

		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));	

		String wd = params.getProperty("working_directory");

		String networkType = params.getProperty("network_name");

		String networkFile = wd + params.getProperty("networkRepositoryFile"); 		// set params
		String proteinMappingFile = wd + params.getProperty("proteinInfoMappingFile");	// set params
		String proteinMappingFile2 = wd +  "ioFiles/" + networkType + "_proteinsInNetwork_info.tsv";

		String originalDistMatrixFile = wd +  "ioFiles/" + networkType + "_initial_distanceMatrix.tsv";
		String connectedComponentDistMatrixFile = wd +  "ioFiles/" + networkType + "_connectedComponent_distanceMatrix.tsv";

		String annotationFile = wd + params.getProperty("annotationFile"); 	 	// set params 
		String proteinAnnotationFrequencyFile = wd +  "ioFiles/" + networkType + "_proteinAnnotationFrequency.tsv";
		String annotationCompanionFile = wd + "ioFiles/" + networkType + "_annotationCompanionFile.tsv";

		String mcSamplingPrefix = wd + "mcSampling/" + networkType + "_mcSamplingDistribution_";
		String normalDistributionParamsFile = wd + "ioFiles/" + networkType + "_normalDistributionParams.tsv";

		String annotationOutputFile = wd + networkType + "_structure_clusteringDetails.tsv";
		String significanceScoreFile = wd + networkType + "_structure_significantScores.tsv";

		int clusteringMeasure = Integer.parseInt(params.getProperty("clusteringMeasure"));
		double clusteringThreshold = Double.parseDouble(params.getProperty("clusteringThreshold"));

		/* Load interaction network - from ABC format */ 
		System.out.println("** Loading interaction repository **");
		ArrayList<Interaction> interactionList = CorrelationGraphLoader.loadGraphFromCorrelationNetwork(networkFile, proteinMappingFile);
		System.out.println("Number of interactions:" + interactionList.size() + "\n");

		/* Determine proteins in network  */ 
		System.out.println("** Getting list of proteins in network **");
		ArrayList<Protein> proteinList = NetworkProteins.getProteinsInNetwork(interactionList);
		System.out.println("Number of Proteins: " + proteinList.size() + "\n");

		/* Generate distance matrix of connected component */
		double[][] distanceMatrix = new double[1][1];
		ArrayList<Protein> proteinList2 = new ArrayList<>();

		File f = new File(originalDistMatrixFile);
		File f1 = new File(connectedComponentDistMatrixFile);
		if(!f1.exists() && !f1.isDirectory()) {

			if(!f.exists() && !f.isDirectory()) {
				System.out.println("** Generating initial distance matrix **");
				distanceMatrix = DistanceMatrix.computeDistanceMatrix(interactionList, proteinList, originalDistMatrixFile);

			} else {
				System.out.println("** Loading distance matrix **");
				distanceMatrix = DistanceMatrix.loadDistanceMatrix(originalDistMatrixFile, proteinList); 
			}

			/* Determine which proteins are disconnected*/ 
			boolean[] proteinsToKeep = Calculator.determineConnectedProteins(distanceMatrix);

			/* Update proteins to keep in the network (ie. those that contribute to the most connected component) */
			proteinList2 = NetworkProteins.modifyNetworkProteinsList(proteinList, proteinsToKeep, proteinMappingFile2);

			/* Update distance matrix to only include proteins that contribute to the most connected component of the graph */
			System.out.println("** Updating distance matrix for connected component **");
			distanceMatrix = DistanceMatrix.updateDistanceMatrix(proteinsToKeep, distanceMatrix, connectedComponentDistMatrixFile);

		} else {
			/* Load distance matrix representing fully connected component */
			System.out.println("** Connected Component distance matrix **\n");
			proteinList2 = NetworkProteins.loadProteinsInNetwork(proteinMappingFile2);
			distanceMatrix = DistanceMatrix.loadDistanceMatrix(connectedComponentDistMatrixFile, proteinList2);
		}

		/* Parameters for enrichment */
		int lowerBound = Integer.parseInt(params.getProperty("lowerBoundToSample", "3"));
		int upperBound = Integer.parseInt(params.getProperty("upperBoundToSample", "2000"));
		int numOfSamplings = Integer.parseInt(params.getProperty("numberOfSamplings"));

		/* accessory set */
		HashSet<String> proteinSet = NetworkProteins.getProteinSet(proteinList2);

		/* scale scores */
		if(clusteringMeasure == 3) {
			String scaledAnnotationFile = wd + "ioFiles/" + networkType + "_scaledAnnotations.tsv";
			
			Scale.scaleScores(annotationFile, scaledAnnotationFile);
			annotationFile = scaledAnnotationFile;
		}
		
		if(clusteringMeasure == 4) {
			String percentileAnnotationFile = wd + "ioFiles/" + networkType + "_percentileAnnotations.tsv";
			String percentilesFile = wd + params.getProperty("percentilesFile");
			
			double threshold = 0.75;
			
			ScorePercentiles.assignScorePercentile(annotationFile, percentilesFile, percentileAnnotationFile, threshold);
			annotationFile = percentileAnnotationFile;
		}

		f = new File(annotationCompanionFile);
		f1 = new File(proteinAnnotationFrequencyFile);

		if(!f1.exists() && !f1.isDirectory()) {

			if(!f.exists() && !f.isDirectory()) {
				System.out.println("** Generating companion file **\n");
				AnnotationCompanionFile.determineAnnotatedProteinsInNetwork(annotationFile, annotationCompanionFile, proteinSet, lowerBound, upperBound);
			}

			/* Calculate protein annotation frequency */
			System.out.println("** Enumerating protein annotation frequency files **\n");
			AssessProteinAnnotations freq = new AssessProteinAnnotations(lowerBound, upperBound, proteinSet);
			freq.computeFrequencyOfProteins(annotationFile, annotationCompanionFile, proteinAnnotationFrequencyFile);
		}

		/* Perform Monte Carlo Sampling procedure */
		if(Boolean.parseBoolean(params.getProperty("performMCprocedure"))) {

			System.out.println("**Performing Monte Carlo Sampling Procedure**");

			/* compute distribution by number of proteins to sample */
			if(Integer.parseInt(params.getProperty("samplingMethod"))==0) {
				MotifSampling sampling = new MotifSampling(proteinAnnotationFrequencyFile, proteinList2, distanceMatrix, clusteringMeasure, clusteringThreshold); // 2 - Initialize sampling
				sampling.computeMultipleDistributions(Integer.parseInt(args[1]), Integer.parseInt(args[2]), numOfSamplings, mcSamplingPrefix, annotationCompanionFile); // 3 - Perform sampling for n proteins

			} else { 
				/* compute distribution per annotation */ 
				MotifSamplingPerModule sampling = new MotifSamplingPerModule(proteinAnnotationFrequencyFile, proteinList2, distanceMatrix, clusteringMeasure, clusteringThreshold);
				sampling.computeMCdistributionsForAllModules(annotationCompanionFile, annotationFile, numOfSamplings, mcSamplingPrefix, Integer.parseInt(args[1]));
			}
		}

		if(Boolean.parseBoolean(params.getProperty("calculateNormalDistributionParams"))) {
			if(Integer.parseInt(params.getProperty("samplingMethod"))==0) {
				ApproximateNormalDistribuiton.getNormalDistributionParams(mcSamplingPrefix, lowerBound, upperBound, numOfSamplings, normalDistributionParamsFile);
			} else { 
				ApproximateNormalDistribuiton.getNormalDistributionParamsByModule(annotationCompanionFile, mcSamplingPrefix, Integer.parseInt(args[1]), Integer.parseInt(args[2]), numOfSamplings, normalDistributionParamsFile);
			}
		}

		/* Load and test significance annotations */
		if(Boolean.parseBoolean(params.getProperty("testMotifs"))) {
			System.out.println("**Assessing motif clustering**");
			MotifEnrichment m = new MotifEnrichment(distanceMatrix, proteinList2, normalDistributionParamsFile, lowerBound, upperBound,  clusteringMeasure, clusteringThreshold);
			m.testMotifClustering(annotationFile, annotationCompanionFile, annotationOutputFile, Integer.parseInt(args[1]), Integer.parseInt(args[2]), Integer.parseInt(params.getProperty("samplingMethod")));
		}

		/* Look at the overall distribution of significance scores once all annotations have been tested */
		if(Boolean.parseBoolean(params.getProperty("assessSignificanceScores"))) {
			System.out.println("**Assessing significance scores**");
			AssessEnrichment.assessSignificanceScores(annotationOutputFile, significanceScoreFile);
		}
	}
}
