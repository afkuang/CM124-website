package cs124project3;

import java.util.Random;
import java.util.Vector;
 
public class createData {
     
    public String createSNP(double maf, int snp_size) {
        //Create random number to account for maf
        //Have to use 0-10 instead of 0-1 because of integer properties of rand
        //Therefore, a maf of 0.4 is 4 in the code (4/10 = 0.4)
 
        Random prob = new Random();
        maf *= 10;
        String curSNP = "";
 
        for (int curAllele = 0; curAllele < snp_size; curAllele++) {
            //Returns the next pseudorandom, uniformly distributed int value from this random number generator's sequence.
            int prob_maf = prob.nextInt(10);
 
            if (prob_maf <= maf)
                curSNP += "1";      //Minor allele
            else
                curSNP += "0";      //Major allele
        }
 
        return curSNP;
    }
     
    public void createChildren(double maf, int snp_size, int sample_size, Vector<String> data, boolean related) {
        //Create random number to account for maf
        //Have to use 0-10 instead of 0-1 because of integer properties of rand
        //Therefore, a maf of 0.4 is 4 in the code (4/10 = 0.4)
         
        Vector<String> parentData = new Vector<String>();
 
        Random prob = new Random();
        double prob_00 = maf*maf;
        double prob_01 = 2*(1-maf)*maf;
        double prob_11 = (1-maf)*(1-maf);
        
        String parentSNP = "";
        String childSNP = "";
 
        for (int curPerson = 0; curPerson < sample_size; curPerson++) {
 
            parentSNP = "";
 
            for (int curAllele = 0; curAllele < snp_size; curAllele++) {
                //Returns the next pseudorandom, uniformly distributed int value from this random number generator's sequence.
                int int_prob_maf = prob.nextInt(100);
                double prob_maf = ((double) int_prob_maf/(double) 100);
                 
                if (prob_maf <= prob_00)
                    parentSNP += "0";       //Individual has 00
                else if (prob_00 < prob_maf && prob_maf <= prob_01)
                    parentSNP += "1";       //Individual has 01 or 10
                else
                    parentSNP += "2";       //Individual has 11
            }
             
            //Finished generating SNP, add it to data.
            //System.out.println("Parent: " + parentSNP);
            parentData.add(parentSNP);
        }
         
        String parent1;
        String parent2;
        char p1Allele;
        char p2Allele;
        for (int curChild = 0; curChild < sample_size; curChild++){
            //Reset for new child
            childSNP = "";
             
            if (!related) {
                parent1 = parentData.elementAt(curChild);
                curChild++;
                parent2 = parentData.elementAt(curChild);
            }
             
            else { //Re-use parents 0 and 1, instead of using 2 and 3.
                parent1 = parentData.elementAt(curChild % 2);
                curChild++;
                parent2 = parentData.elementAt(curChild % 2);
            }
             
            int probAllele = 0;
            for (int curAllele = 0; curAllele < snp_size; curAllele++) {
                 
                //Generate child's inherited allele from parent 1
                if (parent1.charAt(curAllele) == '0')
                    p1Allele = '0';
                else if (parent1.charAt(curAllele) == '1'){
                    probAllele = prob.nextInt(10);
                    if (probAllele < 5)
                        p1Allele = '0';
                    else
                        p1Allele = '1';
                }
                else    //parent1 is a 2 (11 homozygous)
                    p1Allele = '1';
                 
                //Generate child's inherited allele from parent 2
                if (parent2.charAt(curAllele) == '0')
                    p2Allele = '0';
                else if (parent2.charAt(curAllele) == '1'){
                    probAllele = prob.nextInt(10);
                    if (probAllele < 5)
                        p2Allele = '0';
                    else
                        p2Allele = '1';
                }
                else    //parent1 is a 2 (11 homozygous)
                    p2Allele = '1';
                 
                //Generate child's actual inherited allele from combination of p1 and p2
                if (p1Allele == '0' && p2Allele == '0')
                    childSNP += '0';
                else if (p1Allele == '1' && p2Allele == '1')
                    childSNP += '2';
                else {
                    childSNP += '1';
                }
            }
             
            //System.out.println("Child:  " + childSNP);
            data.add(childSNP);
            //System.out.println(data.size());
        }
         
        return;
    }
     
    public void breed(double maf, int snp_size, int num_child, Vector<String> data){
         
                 
        //Randomly generate one mother and father which will be used to generate related children.
        String mother = createSNP(maf, snp_size);
        String father = createSNP(maf, snp_size);
         
//      System.out.println("Mother: " + mother);
//      System.out.println("Father: " + father);
         
        Random prob = new Random();
        String child = "";
        char m_allele;
        char f_allele;
         
        //Generate (breed) each child SNP
        for (int curChild = 0; curChild < num_child; curChild++) {
             
            for (int curAllele = 0; curAllele < snp_size; curAllele++) {
                m_allele = mother.charAt(curAllele);
                f_allele = father.charAt(curAllele);
                 
                int prob_maf = prob.nextInt(10);
                if (prob_maf <= 5) 
                    child += m_allele;
                else
                    child += f_allele;
            }
         
            //Add child to data vector
            data.add(child);
            child = "";
             
        }
         
        return;
    }
}