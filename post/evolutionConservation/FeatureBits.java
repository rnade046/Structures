package evolutionConservation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

public class FeatureBits {

	private double trials;
	private double probSuccess;
	private double motifCons;
	private double mean;
	private double stDev;
	private double foldChange;
	private double pval;

	public FeatureBits(String moduleFBFile) {

		/* load file - set trials, probSuccess and utr conservation */
		loadModuleFile(moduleFBFile);
		this.foldChange = this.motifCons / this.trials;
		
		/* compute mean and stDev */
		approximateValuesForNormalDistribution();
		
		/* compute conservation probability */
		this.pval = computeNormPvalue(this.mean, this.stDev, this.motifCons, this.trials);
	}

	private void loadModuleFile(String file) {

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(file))));
			
			String line = in.readLine();
			int count = 1;
			
			double utr = 0, utrCons = 0;
			
			while(line!=null) {
				
				if(count%2==0) {
					switch(count){
					case 2: 	// 3utr
						utr = obtainFeatureBitsValues(line);
						break;
					case 4:		// 3utr-conservation
						utrCons = obtainFeatureBitsValues(line);
						break;
					case 6:		// motif-conservation
						this.motifCons = obtainFeatureBitsValues(line);
						break;
					case 8:		// motif = trials
						this.trials = obtainFeatureBitsValues(line);
						break; 
					}
				}
				line = in.readLine();
				count++;
			}
			/* probability of success */ 
			if(utr != 0) {
				this.probSuccess = utrCons / utr;
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private double obtainFeatureBitsValues(String line) {
		
		String[] values = line.split("\\s+");
		double fBits = Integer.parseInt(values[0]); 
		return fBits;
	}
	
	private void approximateValuesForNormalDistribution() {
		
		this.mean = this.trials * probSuccess;
		this.stDev = this.mean * (1-this.probSuccess);
	}
	
	private double computeNormPvalue(double mean, double stdev, double value, double max){
		double prob = 0.0;
		
		for(double j = value; j <= max; j=j+0.0001){
			
			double s2Pi = Math.sqrt(2*Math.PI);
			double sd = Math.sqrt(stdev);
			double exp = (-((j-mean)*(j-mean)))/(2*stdev);
			double add = ((1/(sd*s2Pi))*Math.pow(Math.E,exp)*0.0001);
			prob = prob + add;
		}
		return prob;
	}
	
	public String getConservationProperties() {
		
		return String.valueOf(this.motifCons) +"\t" + String.valueOf(this.trials) + "\t"+ String.valueOf(this.probSuccess) +
				"\t" + String.valueOf(this.foldChange) + "\t" + String.valueOf(this.pval);
	}
}
