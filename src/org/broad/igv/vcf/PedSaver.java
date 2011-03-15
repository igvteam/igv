//package org.broad.igv.vcf;
//
///**
// * Created by IntelliJ IDEA.
// * User: jesse
// * Date: Jul 28, 2010
// * Time: 11:32:22 AM
// * To change this template use File | Settings | File Templates.
// */
//
//import org.broad.igv.feature.tribble.TribbleFeatureSource;
//import org.broad.tribble.Feature;
//import org.broad.tribble.vcf.VCFGenotypeEncoding;
//import org.broad.tribble.vcf.VCFGenotypeRecord;
//import org.broad.tribble.vcf.VCFHeader;
//import org.broad.tribble.vcf.VCFRecord;
//
//import java.io.*;
//import java.util.ArrayList;
//import java.util.Iterator;
//
///**
// * @author jrobinso
// */
//public class PedSaver {
//
//    public static void batch(String[] args) throws IOException {
//        createFeatureTestFile();
//    }
//
//    public static void createFeatureTestFile() throws IOException {
//        TribbleFeatureSource source = new TribbleFeatureSource("/Users/jesse/Desktop/VCF/ASW-YRI-LWK_chr10_gatk.recalibrated.vcf");
//        VCFHeader header = (VCFHeader) source.getHeader();
//        ArrayList<String> samples = new ArrayList<String>(header.getGenotypeSamples());
//        File outFile = new File("test/data/output.ped");
//        BufferedWriter output = new BufferedWriter(new FileWriter(outFile));
//        String chr = "10";
//        int start = 50450980;
//        int end = 50451551;
////        50,451,195-50,451,337
////        chr10:50,450,980-50,451,551
//
//
//        Iterator<Feature> features = source.getFeatures(chr, start, end);
//        ArrayList<Feature> list = new ArrayList<Feature>();
//        while (features.hasNext()){
//            list.add(features.next());
//        }
//
//        for (int i=0; i < samples.size(); i++) {
//            String sample = samples.get(i);
//            String fileLine = "0\t" + sample + "\t0\t0\t0\t0\t";
//
//            for (int j=0; j < list.size(); j++){
//                VCFRecord record = (VCFRecord) list.get(j);
//                VCFGenotypeRecord genotype = record.getGenotype(sample);
//                char b1;
//                char b2;
//                if (genotype.isEmptyGenotype()){
//                    b1 = '0';
//                    b2 = '0';
//                }else if(genotype.isNoCall()){
//                    b1 = '0';
//                    b2 = '0';
//                }else{
//                    String bases = genotype.getBases();
//                    if (bases.length() < 2){
//                        b1 = '0';
//                        b2 = '0';
//                    }else{
//                        b1 = genotype.getBases().charAt(0);
//                        b2 = genotype.getBases().charAt(1);
//                    }
//                }
//                char[] array = {b1, ' ', b2};
//                String bases = new String(array);
//
////                System.out.println(bases);
//                fileLine = fileLine.concat(bases + "\t");
//
//            }
//            output.write(fileLine);
//            output.newLine();
//        }
//        output.close();
//    }
//
//}
