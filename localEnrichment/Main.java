import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

import graph.Interaction;
import graph.Protein;
import network.Calculator;
import network.CorrelationGraphLoader;
import network.DistanceMatrix;
import network.NetworkProteins;

public class Main {
	public static void main(String[] args) {
		
		String networkFile = "";
		String fastaFile = "";
		String biomartMappingFile = "";
		String proteinMappingFile = "";
		
		String originalDistMatrixFile = "";
		String connectedComponentDistMatrixFile = "";
		
		/* Load interaction network - from ABC format */ 
		System.out.println("** Loading interaction repository **");
		ArrayList<Interaction> interactionList = CorrelationGraphLoader.loadGraphFromCorrelationNetwork(networkFile, fastaFile, biomartMappingFile);
		System.out.println("Number of interactions:" + interactionList.size() + "\n");

		/* Determine proteins in network  */ 
		System.out.println("** Getting list of proteins in network **");
		ArrayList<Protein> proteinList = NetworkProteins.getProteinsInNetwork(interactionList);
		System.out.println("Number of Proteins: " + proteinList.size() + "\n");
		
		
	
		File f = new File(originalDistMatrixFile);
		if(!f.exists() && !f.isDirectory()) {
			System.out.println("** Generating initial distance matrix **");
			DistanceMatrix.computeDistanceMatrix(interactionList, proteinList, originalDistMatrixFile);
		}

		/* Perform motif enumeration around here */

		System.out.println("**Loading distance matrix**");
		double[][] distanceMatrix = DistanceMatrix.loadDistanceMatrix(originalDistMatrixFile, proteinList); 

		/* Determine which proteins are disconnected*/ 
		boolean[] proteinsToKeep = Calculator.determineConnectedProteins(distanceMatrix);
		/* Update proteins to keep in the network (ie. those that contribute to the most connected component) */
		ArrayList<Protein> proteinList2 = NetworkProteins.modifyNetworkProteinsList(proteinList, proteinsToKeep, proteinMappingFile);
		HashSet<String> proteinSet = NetworkProteins.getProteinSet(proteinList2);

		//CheckDegreeDistributions.assessDegreeDistribution(proteinList2, interactionList, (wd + projectName + "_degreesInNetwork.tsv"));

		
		if(proteinList.size() != proteinList2.size()) {
			System.out.println("**Checking for disconnected components**");
			File f1 = new File(connectedComponentDistMatrixFile);
			if(!f1.exists() && !f1.isDirectory()) {
				/* Update distance matrix to only include proteins that contribute to the most connected component of the graph */
				System.out.println("**Updating distance matrix**");
				DistanceMatrix.updateDistanceMatrix(proteinsToKeep, distanceMatrix, connectedComponentDistMatrixFile);
			}

			/* Load distance matrix representing fully connected component */
			System.out.println("**Loading updated distance matrix**\n");
			distanceMatrix = DistanceMatrix.loadDistanceMatrix(connectedComponentDistMatrixFile, proteinList2);
		} 
		
		
		
	}
}
