package utils;

import java.util.HashMap;

import sampling.ApproximateNormalDistribuiton;

public class MotifEnrichmentFromMCdist {

	public static double getPvaluesFromMonteCarloDistribution(String distributionInputFile, double clusteringMeasure) {
		
		double pValue = 1;	

		/* load monte carlo distribution */
		HashMap<Double, Double> distributionMap = ApproximateNormalDistribuiton.loadMonteCarloDistributions(distributionInputFile);

		/* if line is empty; confirm that distribution totals 1; compute mean and standard deviation of distribution= */
		if(verifySumOfFrequencies(distributionMap)) {

			double pval = compute_pValue(distributionMap, clusteringMeasure);
			if(pval < pValue && pval != 0 ) {
				pValue = pval;
			} 
		}
		return pValue;
	}

	/**
	 * Compute the sum of frequencies for a given distribution, if it isn't equal to 1, make error. 
	 * 
	 * @param distributionMap	
	 * @return frequenciesVerified	 boolean; if true, sum of frequencies is 1
	 */
	public static boolean verifySumOfFrequencies(HashMap<Double, Double> distributionMap) { 

		double countFrequencies = 0;
		boolean frequenciesVerified = true;

		/* compute the sum of frequencies in the distribution */
		for(Double tpd:distributionMap.keySet()) {
			countFrequencies += distributionMap.get(tpd);
		}

		/* check if sum of frequencies equals 1 */
		if(Math.rint(countFrequencies)!=1) {
			//System.out.println("ERROR loading distribution; frequency = " + countFrequencies);
			frequenciesVerified = false;
		}
		return frequenciesVerified;
	}

	public static double compute_pValue(HashMap<Double, Double> probOfTPDmap, double TPD) {

		double p_value1 = 0; // initialize p-value

		for (double distance : probOfTPDmap.keySet()) {
			if (distance <= TPD) { // frequencies of distances smaller then the obtained TPD contribute to the p-value
				p_value1 += probOfTPDmap.get(distance);
			}

		}
		return p_value1;
	} // end compute p_value

}
