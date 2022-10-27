/*
 * Copyright (C) 2016 weizegan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.alignment.hmm;
import edu.msu.cme.rdp.readseq.SequenceType;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
/**
 *
 * @author weizegan
 */
public class HMMPacBioParserBLAST {
    public static ArrayList readModel(File f, int minLength, int buildNum) throws IOException {
        return readModel(new FileInputStream(f), minLength, buildNum);
    }

    public static ArrayList readModel(InputStream is, int minLength, int buildNum) throws IOException {
        ArrayList hmmList = readModelInternal(is, minLength, buildNum);

        hmmListConfigureGlocal(hmmList,true);
        
        return hmmList;
    }

    public static ArrayList readUnnormalized(File f, int minLength, int buildNum) throws IOException {
        return readUnnormalized(new FileInputStream(f), minLength, buildNum);
    }

    public static ArrayList readUnnormalized(InputStream is, int minLength, int buildNum) throws IOException {
        ArrayList hmm = readModelInternal(is, minLength, buildNum);

        //hmm.configureGlocal(false);

        return hmm;
    }

    private static ArrayList readModelInternal(InputStream is, int minLength, int buildNum) throws IOException {
        ArrayList hmmList = new ArrayList();
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        

        //hmm.version = reader.readLine().trim().split("\\s+")[0];
        /*
        if (!hmm.version.equals("HMMER3/b") && !hmm.version.equals("HMMER3/f")) {
            throw new IOException("Cannot parse " + hmm.version + " version hmmer model, only supports HMMER3/b");
        }
        */
        String line;
        String baseCall=null; // sequence is
        String[] qualityValue = null,deleQV = null,deleTag,insertQV = null ,mergeQV = null,substitueQV=null,substituteTag=null;
        String alpha = null,qualityValue2;
        int order=0;  //the index of the sequence in the file
        long startTime = System.currentTimeMillis();
        while ((line = reader.readLine()) != null) {
            PacBioHMM hmm = new PacBioHMM();
            line = line.trim();
            if (line.startsWith(">")){
                baseCall      = reader.readLine().trim();                         
                qualityValue  = reader.readLine().trim().split("\\s+");  
                deleQV        = reader.readLine().trim().split("\\s+");
                deleTag       = reader.readLine().trim().split("\\s+");
                insertQV      = reader.readLine().trim().split("\\s+");
                mergeQV       = reader.readLine().trim().split("\\s+");
                substitueQV   = reader.readLine().trim().split("\\s+");
                substituteTag = reader.readLine().trim().split("\\s+");
                if (baseCall.length() < minLength){
                    continue;
                }
                order ++;
                hmm.index = order;
            
            }else{
                continue;
            }
            if (hmm.index > buildNum){
                break;
            }               
            //throw new IOException("Unknown sequence alphabet \"" + alpha + "\"");
            hmm.sequenceIs = baseCall;
            if (qualityValue.length > 1){
                hmm.m = qualityValue.length;
            
            }else{
                throw new IOException("The length of sequence is 0 at "  + order + "th" + "\"");
            
            }
            hmm.alphabet =  SequenceType.Nucleotide;
            //String[] lexemes = line.trim().split("\\s+");
            
            parseAlpha(hmm, "A  C  G  T");
            
            hmm.transitions = new double[PacBioHMM.NUM_TRANSITIONS][hmm.m]; // We disard the B->M transition
            hmm.emissions = new double[hmm.m][hmm.k][PacBioHMM.NUM_EMISSION_STATES];
            double[] probabilities;
            
            for (int index = 0; index < hmm.m; index++) {
                if( index == 1508){
                        int aa = 0;
                }
                // for insert emission
                if (index < hmm.m - 1) {
                    //probabilities = readProbabilities(reader.readLine(), 1, hmm.k);
                    //Arrays.fill(hmm.alphaMapping, 0.25);
                    probabilities = insertEmissionProbabilities(String.valueOf(baseCall.charAt(index)), deleQV[index + 1], mergeQV[index], index + 1, hmm.index);
                    for (int emission = 0; emission < hmm.k; emission++) {
                        hmm.emissions[index][emission][PacBioHMM.isc] = probabilities[emission];
                
                    }                 
                    
                }else{
                    for (int emission = 0; emission < hmm.k; emission++) {
                        hmm.emissions[index][emission][PacBioHMM.isc] = 0.0;
                
                    }
                    
                    
                }
                
                // for match emission
                //if (index > 0){
                     probabilities = matchEmissionProbabilities(String.valueOf(baseCall.charAt(index)),substituteTag[index],substitueQV[index],index + 1,hmm.index);
                     for (int emission = 0; emission < hmm.k; emission++) {
                        hmm.emissions[index][emission][PacBioHMM.msc] = probabilities[emission];
                     }
                //}
               
                // write the hmm
                //File hmmOut = new File("/work/weizegan/PacBioViterbiSequence/pacBioHMMProbility.txt");
               // hmmOut.createNewFile();
                //BufferedWriter outPut = new BufferedWriter(new FileWriter(hmmOut));
            
                // for transition
                if (index < hmm.m - 1) {
                    if( index == 1020){
                        int aa = 0;
                    }
                    probabilities = tranProbabilities(deleQV[index + 1],insertQV[index + 1],insertQV[index],mergeQV[index],String.valueOf(baseCall.charAt(index)),index + 1,hmm.index);
                    //outPut.write(index + " ");
                    
                    for (int trans = 0; trans < PacBioHMM.NUM_TRANSITIONS; trans++) {
                        hmm.transitions[trans][index] = probabilities[trans];
                        //outPut.write(String.valueOf(probabilities[trans]));
                    }
                    //outPut.write("\n");
                }else {
                    for (int trans = 0; trans < PacBioHMM.NUM_TRANSITIONS; trans++) {
                        hmm.transitions[trans][index] = Double.NEGATIVE_INFINITY;
                    }
                }
                           
            //outPut.close();     
            }
            // find the max probality of the transitions in the whole model
            
            double maxP = 0.0;
            for (int index = 0; index < hmm.m; index++){
                for (int trans = 0; trans < hmm.NUM_TRANSITIONS; trans++){
                    if (hmm.transitions[trans][index] > maxP){
                        maxP = hmm.transitions[trans][index];
                    }
                }
            }
            hmm.maxProbability = maxP;
            
            // normalize the transition probality
            //for (int index = 0; index < hmm.m; index++){
            //    for (int trans = 0; trans < hmm.NUM_TRANSITIONS; trans++){
             //       hmm.transitions[trans][index] = hmm.transitions[trans][index] + maxP;
                    
             //   }
            //}
            
            
            /*
            // write the hmm
            File hmmOut = new File("/work/weizegan/PacBioViterbiSequence/pacBioHMMProbility.txt");
            hmmOut.createNewFile();
            BufferedWriter outPut = new BufferedWriter(new FileWriter(hmmOut));
            outPut.write("Sequence: ");
            //write sequence
            for (int index = 0; index < hmm.m; index++){
                outPut.write(String.valueOf(baseCall.charAt(index)));
                if (index < hmm.m -1){
                    outPut.write(" ");
                }
            }
            outPut.write("\n");
            //write transition probalitity.
            for (int nn = 0; nn < PacBioHMM.NUM_TRANSITIONS; nn ++){
                outPut.write("MM ");
                for (int index = 0; index <= hmm.m; index++){
                
                }
            
            
            
            }
            for (int index = 0; index <= hmm.m; index++){
                outPut.write(String.valueOf(baseCall.charAt(index)));
                if (index < hmm.m -1){
                    outPut.write(" ");
                }
            }
            
            hmmList.add(hmm);
            
            if ((hmm.index + 1) % 10 == 0) {
                System.err.println("Build " + (hmm.index + 1) + " sequences in " + (System.currentTimeMillis() - startTime)/1000 + " s");
            }   */
            hmmList.add(hmm);
         
        }
             
      


        reader.close();
       

        return hmmList;
    }

    private static void parseAlpha(PacBioHMM hmm, String line) throws IOException {
        String[] lexemes = line.trim().split("\\s+");

        hmm.alphaMapping = new int[127];
        Arrays.fill(hmm.alphaMapping, -1);

        for (int index = 0; index < lexemes.length; index++) {
            if (lexemes[index].length() != 1) {
                throw new IOException("Alphabet symbol " + lexemes[index] + " too long");
            }

            char sym = Character.toUpperCase(lexemes[index].charAt(0));
            hmm.alphaMapping[(int) sym] = index;
            hmm.alphaMapping[(int) Character.toLowerCase(sym)] = index;

            if (sym == 'U' && hmm.alphabet == SequenceType.Nucleotide) {
                hmm.alphaMapping[(int) 'T'] = index;
                hmm.alphaMapping[(int) 't'] = index;
            }
        }

        hmm.k = lexemes.length;
    }
    
    private static double[] insertEmissionProbabilities(String basecallIs, String deleteQV, String mergeQV,int position, int index) throws IOException {
        if (basecallIs == null || deleteQV == null || mergeQV == null) {
            throw new IOException("Unexpected char");
        }
        double delete   =  Math.pow(10,-0.1 * Double.valueOf(deleteQV));
        double merge =  Math.pow(10,-0.1 * Double.valueOf(mergeQV));
        String nucleotide = "ACGT";
        int orderBasecall = nucleotide.indexOf(basecallIs); //the order of the basecallIs
        double[] ret = new double[4];
        ret[orderBasecall] = merge / (merge + delete);
        for (int i = 0; i < 4; i++) {
            if (i == orderBasecall) {
                continue;
            } else {
                ret[i] =  delete / ((merge + delete) * 3.0);
            }
        }
        
        // use blast score strategy
        ret[0] = ret[1] = ret[2] = ret[3] = 0.25;

        return ret;    
    }

    /**
     * Parses the probailities for a model from the provided string
     *
     * @param line
     * @param offset
     * @param expectedValues
     * @return
     */
    private static double[] matchEmissionProbabilities(String basecallIs, String subTag, String subQV,int position, int index) throws IOException {
        if (basecallIs == null || subTag == null || subQV == null) {
            throw new IOException("Unexpected char");
        }
        //String[] lexemes = line.trim().split("\\s+");
        double[] ret = new double[4];

        //if (Nucleotide.orlexemes.length < expectedValues + offset) {
           // throw new IOException("Too few probabilities on line \"" + line + "\"");
        //}
        String nucleotide = "ACGT";
        int orderBasecall = nucleotide.indexOf(basecallIs); //the order of the basecallIs
        int orderSubTag   = -1; //the order of the subTag
        switch (Integer.valueOf(subTag)){
                case 65:
                    orderSubTag = 0;
                    break;
                case 67:
                    orderSubTag = 1;
                    break;
                case 71:
                    orderSubTag = 2;
                    break;
                case 84:
                    orderSubTag = 3;
                    break;
                default:
                    throw new IOException("SubsitituteTag is wrong at +" + (position + 1) +"th " + " in the " + index + "th sequence" + "\"");     
                
        }
                        
        if (orderBasecall < 0 || orderSubTag < 0){
            throw new IOException("Unexpected char, basecallIs = " + basecallIs + ". orderBasecall = " + orderBasecall + "\"");        
        }
        
        if (orderBasecall == orderSubTag) {
            throw new IOException("SubsitituteTag = basecall !! at +" + (position + 1) +"th " + " in the " + index + "th sequence" + "\""); 
        }
        ret[orderBasecall] = 1 -  Math.pow(10,-0.1 * Double.valueOf(subQV));
        // use blast score strategy
        ret[orderBasecall] = 1.0 - 3.0 * 0.0504671584454;
        ret[orderSubTag] =  Math.pow(10,-0.1 * Double.valueOf(subQV)) * 0.5;
        // use blast score strategy
        ret[orderSubTag] =  0.0504671584454;//0.0;
        for (int i = 0; i < 4; i++) {
            if (i == orderBasecall || i == orderSubTag) {
                continue;
            } else {
                ret[i] = Math.pow(10,-0.1 * Double.valueOf(subQV)) * 0.25;
                // use blast score strategy
                ret[i] = 0.0504671584454;//0; 
            }
        }

        return ret;
    }
    
    
    private static double[] tranProbabilities(String deleQV, String insertQV, String insertQVBefore, String mergeQV, String mergeIs,int position, int index) throws IOException {
        if (deleQV == null || insertQV == null || insertQVBefore == null) {
            throw new IOException("Unexpected char at +" + (position + 1) +"th " + " in the " + index + "th sequence" + "\"");
        }
        //String[] lexemes = line.trim().split("\\s+");
        double[] ret = new double[7];

        ret[0] = 0.9865241;
        ret[1] = ret[2] = 0.00673794699;
        ret[3] = ret[5] = 0.9932620530009145;
        ret[4] = ret[6] = 0.00673794699;
        return ret;
    }
    
    
    
    
    

    public static void main(String[] args) throws IOException {
        //ArrayList hmm = HMMPacBioParser.readModel(new File("/work/weizegan/E01_1/Analysis_Results/pacBio1000_p0.1.bax.h5.hmm"));
/*
        System.out.print("double expected_compo[] = {");
        for (int index = 0; index < hmm.compo.length; index++) {
            System.out.print(hmm.compo[index]);
            if (index + 1 != hmm.compo.length) {
                System.out.print(",");
            }
        }
        System.out.println("};");

        System.out.println("double expected_msc[][20] = {");
        for (int state = 0; state < hmm.M() + 1; state++) {
            System.out.print("\t{");
            for (int index = 0; index < hmm.compo.length; index++) {
                if (hmm.msc(state, index) == Double.NEGATIVE_INFINITY) {
                    System.out.print("-std::numeric_limits<double>::infinity()");
                } else {
                    System.out.print(hmm.msc(state, index));
                }
                if (index + 1 != hmm.compo.length) {
                    System.out.print(",");
                }
            }
            System.out.print("}");
            if (state + 1 != hmm.M() + 1) {
                System.out.println(",");
            }
        }
        System.out.println("};");

        System.out.println("double expected_tsc[][7] = {");
        for (int state = 0; state < hmm.M() + 1; state++) {
            System.out.print("\t{");
            for (TSC tsc : new TSC[]{TSC.MM, TSC.MI, TSC.MD, TSC.IM, TSC.II, TSC.DM, TSC.DD}) {

                if (hmm.tsc(state, tsc) == Double.NEGATIVE_INFINITY) {
                    System.out.print("-std::numeric_limits<double>::infinity()");
                } else {
                    System.out.print(hmm.tsc(state, tsc));
                }
                if (tsc != TSC.DD) {
                    System.out.print(",");
                }
            }
            System.out.print("}");
            if (state + 1 != hmm.M() + 1) {
                System.out.println(",");
            }
        }
        System.out.println("};");

        System.out.print("double expected_max_emission[] = {");
        for(int state = 0;state <= hmm.M();state++) {
            System.out.print(hmm.getMaxMatchEmission(state));
            if(state != hmm.M()) {
                System.out.println(", ");
            }
        }
        System.out.println("};");

        System.out.println("double expected_m_hcost = " + hmm.getHCost().computeHeuristicCost('m', 0) + ";");
        System.out.println("double expected_d_hcost = " + hmm.getHCost().computeHeuristicCost('d', 0) + ";");
        System.out.println("double expected_i_hcost = " + hmm.getHCost().computeHeuristicCost('i', 0) + ";");
        */
    }

    private static void hmmListConfigureGlocal(ArrayList hmmlist, boolean b) {
        for (int i = 0; i < hmmlist.size(); i++) {
            PacBioHMM hmm = (PacBioHMM) hmmlist.get(i);
            hmm.configureGlocal(b);
        }
    }
}
