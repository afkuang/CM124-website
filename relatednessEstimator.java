package cs124project3;

import java.io.IOException;
import java.util.Random;
import java.util.Vector;
 
public class relatednessEstimator {
     
    public static boolean baselineRE(String person1, String person2, int snp_size) {
        int numMatches = 0;
        for (int curAllele = 0; curAllele < snp_size; curAllele++) {
            if (person1.charAt(curAllele) == person2.charAt(curAllele))
                numMatches++;
        }
         
        if (((double) numMatches/ (double) snp_size) > 0.50){
            return true;
        }
        else{
            return false;
        }
    }
    
    public static void createRelationMatrices(double maf, double[][] unrelated, double[][] related) {
    	
    	unrelated[0][0] = Math.pow((1-maf), 4.0);
		unrelated[0][1] = 2 * Math.pow((1-maf), 3.0) * maf;
		unrelated[0][2] = Math.pow((1-maf), 2.0) * Math.pow(maf, 2.0);
    	unrelated[1][0] = 2 * Math.pow((1-maf), 3.0) * maf;;
		unrelated[1][1] = 4 * Math.pow((1-maf), 2.0) * Math.pow(maf, 2.0);
		unrelated[1][2] = 2 * Math.pow((maf), 3.0) * (1-maf);
    	unrelated[2][0] = Math.pow((1-maf), 2.0) * Math.pow(maf, 2.0);
		unrelated[2][1] = 2 * Math.pow(maf, 3.0) * (1-maf);
		unrelated[2][2] = Math.pow(maf, 4.0);
		
    	related[0][0] = 1.0 * 1.0 * unrelated[0][0] +
    					0.5 * 0.5 * unrelated[0][1] +
    					0.5 * 0.5 * unrelated[1][0] +
    					0.25 * 0.25 * unrelated[1][1];
    	
		related[0][1] = 0.5 * 0.5 * unrelated[0][1] +
						0.5 * 0.5 * unrelated[1][0] +
						0.5 * 0.25 + unrelated[1][1];
		
		related[0][2] = 0.25 * 0.25 * unrelated[1][1];
    	
		related[1][0] = 0.5 * 0.5 * unrelated[0][1] +
    					0.5 * 0.5 * unrelated[1][0] +
    					0.5 * 0.25 * unrelated[1][1];
		
		related[1][1] = 0.5 * 0.5 * unrelated[0][1] +
						0.5 * 0.5 * unrelated[1][0] +
						0.5 * 0.5 * unrelated[1][1] +
						1.0 * 1.0 * unrelated[0][2] +
						1.0 * 1.0 * unrelated[2][0] +
						0.5 * 0.5 * unrelated[1][2] +
						0.5 * 0.5 * unrelated[2][1];
		
		related[1][2] = 0.5 * 0.5 * unrelated[1][2] +
						0.5 * 0.5 * unrelated[2][1] +
						0.5 * 0.25 * unrelated[1][1];
    	
		related[2][0] = 0.25 * 0.25 * unrelated[1][1];
		
		related[2][1] = 0.5 * 0.5 * unrelated[1][2] +
						0.5 * 0.5 * unrelated[2][1] +
						0.5 * 0.25 * unrelated[1][1];
		
		related[2][2] = 0.5 * 0.5 * unrelated[1][2] +
						0.5 * 0.5 * unrelated[2][1] +
						1.0 * 1.0 * unrelated[2][2] +
						0.25 * 0.25 * unrelated[1][1];
    }
    
    public static double computeRelatedness(String person1, String person2, double[][] relationMatrix, int snp_size) {
    	double runningProb = 0.0;
    	
    	for (int curAllele = 0; curAllele < snp_size; curAllele++) {
    		int index1 = Integer.parseInt(String.valueOf(person1.charAt(curAllele)));
    		int index2 = Integer.parseInt(String.valueOf(person2.charAt(curAllele)));
    		runningProb += relationMatrix[index1][index2];
    	}
    	
    	double runningProbAvg = runningProb/((double) snp_size);
    	return runningProbAvg;
    }
     
    public static void main(String[] args) throws IOException {
         
        int num_trials = 100;
        int num_correct_base = 0;
        int num_correct_improved = 0;
        double base_accuracy = 0;
        double improved_accuracy = 0;
        int snp_size = 100;
        double maf = 0.4;           //minor allele frequency
         
        Vector<String> child_data = new Vector<String>();
        createData cd = new createData();
//      cd.createUnrelatedChildren(maf, snp_size, sample_size, child_data);
//      
//      cd.breed(maf, snp_size, 2, child_data);
//      
//      for(int i = 0; i < child_data.size(); i++) {
//          System.out.println((i+1) + ": " + child_data.elementAt(i));
//      }
         
         
        /*---------------------------------------------------------------------------------------------*/
        // Calculate baseline accuracy
         
        Random prob = new Random();
        
        double averageAccuracy = 0.0;
             
        for (int curTest = 0; curTest < 50; curTest++) {
        	
        	base_accuracy = 0.0;
        	num_correct_base = 0;
        	
	        for (int curTrial = 0; curTrial < num_trials; curTrial++) {
	            int prob_related = prob.nextInt(2);
	            if (prob_related == 0) { //Unrelated
	                //System.out.println("Children should be unrelated");
	                cd.createChildren(maf, snp_size, 4, child_data, false);
	                if (!baselineRE(child_data.elementAt(0), child_data.elementAt(1), snp_size))
	                    num_correct_base++;
	                child_data.removeAllElements();
	            }
	            else { //Related
	                //System.out.println("Children should be related");
	                cd.createChildren(maf, snp_size, 4, child_data, true);
	                if (baselineRE(child_data.elementAt(0), child_data.elementAt(1), snp_size))
	                    num_correct_base++;
	                child_data.removeAllElements();
	            }
	        }
	        
	        base_accuracy = ((double) num_correct_base/(double) num_trials);
//	        System.out.println("Baseline Accuracy: " + base_accuracy*100 + "%");
	        
	        averageAccuracy += base_accuracy;
	        averageAccuracy /= 2;
        }
                
        /*---------------------------------------------------------------------------------------------*/
        
        /*---------------------------------------------------------------------------------------------*/
        // Calculate improved algorithm accuracy
        
        //Initialize two "empty" 3x3 matrices
        double[][] unrelated = {{0.0, 0.0, 0.0}, 
        						{0.0, 0.0, 0.0}, 
        						{0.0, 0.0, 0.0}};
        double[][] related = {	{0.0, 0.0, 0.0}, 
								{0.0, 0.0, 0.0}, 
								{0.0, 0.0, 0.0}};
				        
        //Fill matrices with probabilities for each event given the MAF for the trial
        createRelationMatrices(maf, unrelated, related);
        
        double avgRelatedProb = 0.0;
        double avgUnrelatedProb = 0.0;
        
        double averageImprovedAccuracy = 0.0;
        
        for (int curTest = 0; curTest < 50; curTest++) {
        	
        	improved_accuracy = 0.0;
        	num_correct_improved = 0;
       
	        for (int curTrial = 0; curTrial < num_trials; curTrial++) {
	            int prob_related = prob.nextInt(2);
	            if (prob_related == 0) { //Unrelated
	                //System.out.println("Children should be unrelated");
	                cd.createChildren(maf, snp_size, 4, child_data, false);
	                avgRelatedProb = computeRelatedness(child_data.elementAt(0), child_data.elementAt(1), related, snp_size);
	                avgUnrelatedProb = computeRelatedness(child_data.elementAt(0), child_data.elementAt(1), unrelated, snp_size);
	                if (avgRelatedProb < avgUnrelatedProb)
	                    num_correct_improved++;
	                child_data.removeAllElements();
	            }
	            else { //Related
	                //System.out.println("Children should be related");
	                cd.createChildren(maf, snp_size, 4, child_data, true);
	                avgRelatedProb = computeRelatedness(child_data.elementAt(0), child_data.elementAt(1), related, snp_size);
	                avgUnrelatedProb = computeRelatedness(child_data.elementAt(0), child_data.elementAt(1), unrelated, snp_size);
	                if (avgUnrelatedProb < avgRelatedProb)
	                    num_correct_improved++;
	                child_data.removeAllElements();
	            }
	        }
	        
	        improved_accuracy = ((double) num_correct_improved/(double) num_trials);
		      
	        averageImprovedAccuracy += improved_accuracy;
	        averageImprovedAccuracy /= 2;
	        
        }
        

        //Compare the two accuracies.
        System.out.println("Baseline Average Accuracy: " + averageImprovedAccuracy*100 + "%");
        System.out.println("Improved Average Accuracy: " + averageAccuracy*100 + "%");
        /*---------------------------------------------------------------------------------------------*/
         
    }
}