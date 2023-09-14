package sampling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import graph.Protein;
import utils.ClusteringMeasure;

public class MotifSamplingPerModule extends MotifSampling {

	public MotifSamplingPerModule(String annotationFreqFile, ArrayList<Protein> protList, double[][] dm,
			int clustering_measure, double percent_threshold) {
		super(annotationFreqFile, protList, dm, clustering_measure, percent_threshold);
	}

	public void assessMCdistributions(String annotationCompanionFile, String annotationFile, int samplingFreq,  String mcFilePrefix) { 

		/* determine annotations that are tested - and therefore will need a MC distribution */
		HashSet<String> motifsToTest = loadMotifsToTest(annotationCompanionFile);

		/* iterate over annotations & get list of scores */
		try {
			InputStream in = new FileInputStream(new File(annotationFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header
			int moduleCount = 0;
			while(line!=null) {

				String col[] = line.split("\t");

				if(motifsToTest.contains(col[0])) {
					
					/* obtain scores */
					String[] elements = col[2].split("\\|");
					List<Double> scores = new ArrayList<>();
					
					for(String e : elements) { // e = ProteinName_Score
						scores.add(Double.parseDouble(e.split("\\_")[1]));
					}

					/* obtain MC sampling distribution */
					HashMap<Double, Double> distribution = computeMCdistributionForCurrentAnnotation(scores, samplingFreq);
					String mcFile = mcFilePrefix + "s" + samplingFreq + "_n" + moduleCount;
					printDistribution(mcFile, distribution, moduleCount);
				}
				line = input.readLine();
				moduleCount++;
			}
			System.out.println("Done\n");
			input.close();
		} catch (IOException e) {

		}

	}

	private HashMap<Double, Double> computeMCdistributionForCurrentAnnotation(List<Double> scores, int samplingFrequency) {
		
		
		HashMap<Double, Double> distribution;
		double[] tpdSampleList = new double[samplingFrequency]; // array to store the results of sampling {list of TPDs}
		
		for(int i=0; i<samplingFrequency; i++) {
			
			/* Select random proteins in network */
			ArrayList<Integer> randomProteins = getRandomWeightedProteinsWithBinarySearch(scores.size());
			
			/* Assign scores to random proteins by shuffling list order */
			Collections.shuffle(scores);
			
			/* compute the clustering measure */ 
			tpdSampleList[i] = ClusteringMeasure.computeWNodeTPD(randomProteins, scores, distanceMatrix);
		}
	
		/* Measure the frequencies of sampled total pairwise distances from x amount of proteins */
		distribution = Sampling.computeFrequenciesOfSampledTPDs(tpdSampleList, samplingFrequency);

		/* Sanity check to ensure the sum of TPD frequencies equals 1 */
		Sampling.checkFrequencyTotal(distribution, samplingFrequency);
		
		return distribution;
	}
	
	private HashSet<String> loadMotifsToTest(String annotationCompanionFile){

		HashSet<String> motifSet = new HashSet<>();

		try {
			InputStream in = new FileInputStream(new File(annotationCompanionFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				motifSet.add(line.split("\t")[0]);
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}
}