/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.vcf;

import org.broad.igv.feature.*;
import org.broad.igv.track.tribble.TribbleFeatureSource;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;
import org.broad.tribble.iterators.CloseableTribbleIterator;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Experimental utilities, primarily for project "Sigma"
 * <p/>
 * 1) SNP
 * 2) coding / noncoding
 * 3) missense or not
 * 4) concoordinate withitn in families yes or no, add # individuals in column header
 * 5) gene name
 * 6) intron / exon / intergenic
 * 7) conservation
 * 8) TFBS
 * 9) VNTRs
 * hg18: 152,668,596 -> 154,746,948
 *
 * @author jrobinso
 * @date Jan 30, 2011
 */
public class SigmaUtils {

    static String chr = "chr1";
    static int roiStart = 152668595;
    static int roiEnd = 154746948;


    static String refSeqURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/sigma/refGene.hg18.gz";
    static String ucscURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/sigma/ucscGene.hg18.gz";
    static String rnaGeneURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/sigma/rnaGene.hg18.bed.txt";
    static String tfbsURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/sigma/tfbsConserved.hg18.bed.gz";
    static String dbSnpURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/dbSnp/dbsnp_130_hg18.bed.gz";
    static String phastconURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/conservation/phastCons44way.wig.tdf";
    static String pi12merURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/conservation/pi.12mer.wig.tdf";
    static String vntrURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/sigma/VNTR_coding_scores.bed";
    static String codingRepeatsURL = "http://iwww.broadinstitute.org/igvdata/sigma/annotations/CodingRepeats.bed";
    static String cufflinksURL = "http://iwww.broadinstitute.org/igvdata/sigma/annotations/kidney.Cufflinks.hg19ToHg18.chr1.bed";
    static String scriptureURL = "http://iwww.broadinstitute.org/igvdata/sigma/annotations/kidney.Scripture.hg19ToHg18.chr1.bed";
    static String trfURL = "http://iwww.broadinstitute.org/igvdata/sigma/annotations/simple.repeats.bed";
    /*
   kidney.Cufflinks.hg19ToHg18.bed
  kidney.Scripture.hg19ToHg18.bed
  kidneyScriptureCufflinksHG19.toHG18.bed
    */


    static String phastconsPrimateURL = "/Users/jrobinso/Sigma/phastcon44.bed";
    static String phastconsVeterbrateURL = "/Users/jrobinso/Sigma/phastconVerterbrates44.bed";
    static String snp130URL = "/Users/jrobinso/Sigma/snp130.bed";
    static String indelFilterURL = "/Users/jrobinso/Sigma/indel_filter.bed";


    TribbleFeatureSource source;

    Map<String, List<String>> samples;
    private VCFHeader header;

    private List<Feature> refSeqGenes;
    private List<Feature> tfbsFeatures;
    private List<Feature> vntrFeatures;
    private List<Feature> cufflinksFeatures;
    private List<Feature> scriptureFeatures;
    private List<Feature> trfFeatures;
    private List<Feature> ucscGenes;
    private List<Feature> phastPrimateFeatures;
    private List<Feature> phastVertFeatures;
    private List<Feature> snp130Features;
    private List<Feature> indelFilterFeatures;
    private List<Feature> codingRepeatFeatures;

    static Set<String> segregatingFamilies = new HashSet(Arrays.asList("AA", "L", "LA", "BIP", "OK", "CC", "PA"));

    static String[] allFamiliesList = {"AA", "L", "LA", "BIP", "OK", "CC", "PA"};

    String unaffected;

    static String fileBase;

    public static void main(String[] args) throws IOException {

        SigmaUtils utils = new SigmaUtils(args[0]);

        fileBase = args[1];
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(args[1] + ".bed")));
        PrintWriter excelWriter = new PrintWriter(new BufferedWriter(new FileWriter(args[1] + ".xls")));


        excelWriter.println("Type\tStart\tLocus\tID\tFamilies\tGene\tCoding\tFunctional Class\tAA change\tPosition class\t" +
                "UCSC Genes\tCufflinks\tScripture\tTRF\tCodingRepeat\tVNTR\tTFBS\tCons Primates\tCons Vertebrates\tDiscordant families\t" +
                "Disc count\tAA (3 C)\tL (2 C)\tLA (4 C)\tBIP (2 S)\tOK (2 S)\tCC (1)\tPA (1)");


        pw.println("#coords=1");
        pw.println("track name=\"Segregating snps\"");
        utils.walkVariants(pw, excelWriter);

        pw.close();
        excelWriter.close();


    }

    private void initFeatureLists() throws IOException {

        refSeqGenes = getFeatures(refSeqURL);
        tfbsFeatures = getFeatures(tfbsURL);
        vntrFeatures = getFeatures(vntrURL);
        cufflinksFeatures = getFeatures(cufflinksURL);
        scriptureFeatures = getFeatures(scriptureURL);
        trfFeatures = getFeatures(trfURL);
        ucscGenes = getFeatures(ucscURL);
        codingRepeatFeatures = getFeatures(codingRepeatsURL);
        phastPrimateFeatures = getFeatures(phastconsPrimateURL);
        phastVertFeatures = getFeatures(phastconsVeterbrateURL);
        snp130Features = getFeatures(snp130URL);
        indelFilterFeatures = getFeatures(indelFilterURL);
    }

    private List<Feature> getFeatures(String url) {

        AsciiLineReader reader = null;
        try {
            ResourceLocator rl = new ResourceLocator(url);
            reader = ParsingUtils.openAsciiReader(rl);
            List<Feature> tmp = AbstractFeatureParser.getInstanceFor(rl).loadFeatures(reader);
            List<Feature> features = new ArrayList(tmp.size());
            for (Feature f : tmp) {
                if (f.getChr().equals("chr1")) {
                    features.add(f);
                }
            }
            FeatureUtils.sortFeatureList(features);
            return features;
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    public SigmaUtils(String path) throws IOException {

        initFeatureLists();

        source = new TribbleFeatureSource(path);

        header = (VCFHeader) source.getHeader();

        Set<String> allSamples = header.getGenotypeSamples();
        samples = new LinkedHashMap();
        if (UIConstants.isSigmaProject()) {
            samples.put("AA", new ArrayList());
            samples.put("BIP", new ArrayList());
            samples.put("L", new ArrayList());
            samples.put("LA", new ArrayList());
            samples.put("OK", new ArrayList());
            //samples.put("S", new ArrayList());
            samples.put("WC", new ArrayList());
            samples.put("CC", new ArrayList());
            samples.put("WC", new ArrayList());
            samples.put("PA", new ArrayList());
            samples.put("Other", new ArrayList());

            for (String s : header.getGenotypeSamples()) {
                if (s.endsWith("-375") || s.equals("375")) {
                    unaffected = s;
                } else if (s.endsWith("-259") || s.endsWith("-265") || s.endsWith("-266") ||
                        s.equals("259") || s.equals("265") || s.equals("266")) {
                    samples.get("AA").add(s);
                } else if (s.endsWith("-701") || s.endsWith("-564") || s.equals("701") || s.equals("564")) {
                    samples.get("BIP").add(s);
                } else if (s.endsWith("-352") || s.endsWith("-414") || s.equals("352") || s.equals("414")) {
                    samples.get("L").add(s);
                } else if (s.endsWith("-4") || s.endsWith("-8") || s.endsWith("-13") || s.endsWith("-15") ||
                        s.equals("4") || s.equals("8") || s.equals("13") || s.equals("15")) {
                    samples.get("LA").add(s);
                } else if (s.endsWith("-563") || s.endsWith("-566") || s.equals("563") || s.equals("566")) {
                    samples.get("OK").add(s);
                } else if (s.endsWith("-384") || s.endsWith("-391") || s.equals("384") || s.equals("391")) {
                    // Ignore S
                    //   samples.get("S").add(s);
                } else if (s.endsWith("-384") || s.endsWith("-391") || s.equals("384") || s.equals("391")) {
                    samples.get("S").add(s);
                } else if (s.endsWith("-467") || s.equals("467")) {
                    samples.get("PA").add(s);
                } else if (s.endsWith("-469") || s.equals("469")) {
                    samples.get("CC").add(s);
                } else if (s.endsWith("-491") || s.endsWith("-497") || s.equals("491") || s.equals("497")) {
                    samples.get("WC").add(s);
                } else {
                    samples.get("Other").add(s);
                }

            }

        }

    }


    public void walkVariants(PrintWriter bedWriter, PrintWriter excelWriter) throws IOException {


        PrintWriter hetWriter = new PrintWriter(new BufferedWriter(new FileWriter(fileBase + ".hetNormal.bed")));
        hetWriter.println("#coords=1");
        hetWriter.println("trac name=Normal color=200,0,0");

        Map<String, PrintWriter> writers = new HashMap();
        writers.put("AA", new PrintWriter(new BufferedWriter(new FileWriter(fileBase + ".aa.bed"))));
        writers.put("BIP", new PrintWriter(new BufferedWriter(new FileWriter(fileBase + ".bip.bed"))));
        writers.put("LA", new PrintWriter(new BufferedWriter(new FileWriter(fileBase + ".la.bed"))));
        writers.put("L", new PrintWriter(new BufferedWriter(new FileWriter(fileBase + ".l.bed"))));
        writers.put("OK", new PrintWriter(new BufferedWriter(new FileWriter(fileBase + ".ok.bed"))));
        for (Map.Entry<String, PrintWriter> entry : writers.entrySet()) {
            PrintWriter pw = entry.getValue();
            pw.println("#coords=1");
            pw.println("track name=" + entry.getKey());
        }


        CloseableTribbleIterator<Feature> features = source.getFeatures(chr, roiStart, roiEnd);

        List<FilteredVariantRecord> disconcordantRecords = new ArrayList();
        List<FilteredVariantRecord> normalHetRecords = new ArrayList();

        for (Feature feature : features) {

            Set<String> affectedSegregatingFamilies = new LinkedHashSet();
            Set<String> affectedSingles = new LinkedHashSet();
            Set<String> discordantFamilies = new LinkedHashSet();
            boolean normalHet = false;

            VariantContext variant = (VariantContext) feature;

            Genotype genotype = variant.getGenotype(unaffected);
            if (genotype.isHomVar()) {
                continue;
            } else if (genotype.isHet()) {
                discordantFamilies.add("L-375");
                normalHet = true;
                hetWriter.println("chr1\t" + variant.getStart() + "\t" + variant.getEnd() + "\t" + "Het L-375");
            }

            if (variant.getType() == VariantContext.Type.INDEL) {
                Feature tmp = FeatureUtils.getFeatureAt(variant.getStart() - 1, 0, indelFilterFeatures);
                if (tmp == null) {
                    tmp = FeatureUtils.getFeatureAt(variant.getStart() - 1, 0, snp130Features);
                }
                if (tmp != null) {
                    continue;
                }


            }

            // Loop through the families
            for (Map.Entry<String, List<String>> entry : samples.entrySet()) {


                String family = entry.getKey();

                // Skip "S" group
                if (family.equals("S")) {
                    continue;
                }


                // Loop through samples for this familuy
                final List<String> sampleList = entry.getValue();
                if (sampleList.isEmpty()) continue;

                boolean anyHomVar = false;
                boolean allHet = true;
                boolean allRef = true;
                for (String s : sampleList) {
                    genotype = variant.getGenotype(s);
                    String gt = genotype.toString();
                    if (!genotype.isHet()) {
                        allHet = false;
                    }
                    if (!genotype.isHomRef()) {
                        allRef = false;
                    }
                    if (genotype.isHomVar()) {
                        anyHomVar = true;
                        break;
                    }


                }
                if (anyHomVar || !(allHet || allRef)) {
                    discordantFamilies.add(family);
                    //if (discordantFamilies.size() > 1) {
                    // Not this variant
                    //    break;
                    //}
                } else if (allHet) {
                    if(segregatingFamilies.contains(family)) {
                        affectedSegregatingFamilies.add(family);
                    }
                    else {
                        affectedSingles.add(family);
                    }

                }

            }

            if (affectedSegregatingFamilies.size() > 0) {

                FilteredVariantRecord variantRecord = new FilteredVariantRecord(variant, affectedSegregatingFamilies,
                        affectedSingles, discordantFamilies);
                if (discordantFamilies.size() > 0 && discordantFamilies.size() < 2) {
                    if (normalHet) {
                        normalHetRecords.add(variantRecord);
                    } else {
                        disconcordantRecords.add(variantRecord);
                    }
                } else if (discordantFamilies.size() == 0) {
                    writeVariant(bedWriter, excelWriter, variantRecord);
                }
            }

            int s = variant.getStart();
            int e = variant.getEnd();
            String type = variant.getType().toString();
            for (String f : affectedSegregatingFamilies) {
                PrintWriter pw = writers.get(f);
                if (pw != null) {
                    pw.println(variant.getChr() + "\t" + s + "\t" + e + "\t" + type + "\t1000\t.\t" +
                            s + "\t" + e + "\t0,200,0");
                }
            }
            for (String f : discordantFamilies) {
                PrintWriter pw = writers.get(f);
                if (pw != null) {
                    pw.println(variant.getChr() + "\t" + s + "\t" + e + "\t" + type + "\t1000\t.\t" +
                            s + "\t" + e + "\t200,0,0");
                }
            }

        }

        for (PrintWriter pw : writers.values()) {
            pw.close();
        }
        hetWriter.close();

        excelWriter.println();
        excelWriter.println("One discordant family");
        for (FilteredVariantRecord record : disconcordantRecords) {
            writeVariant(bedWriter, excelWriter, record);
        }

        excelWriter.println();
        excelWriter.println("Het call in normal (L-375)");
        for (FilteredVariantRecord record : normalHetRecords) {
            writeVariant(bedWriter, excelWriter, record);
        }
    }

    class FilteredVariantRecord {
        VariantContext variant;
        Set<String> affectedFamilies;
        Set<String> affectedSingles;
        Set<String> discordantFamilies;

        FilteredVariantRecord(VariantContext variant, Set<String> affectedFamilies, Set<String> affectedSingles,
                              Set<String> discordantFamilies) {
            this.affectedSingles = affectedSingles;
            this.affectedFamilies = affectedFamilies;
            this.discordantFamilies = discordantFamilies;
            this.variant = variant;
        }
    }

    private void writeVariant(PrintWriter bedWriter, PrintWriter excelWriter, FilteredVariantRecord record) {

        VariantContext variant = record.variant;
        Set<String> affectedFamilies = record.affectedFamilies;

        String chr = variant.getChr();
        int s = variant.getStart();
        int e = variant.getEnd();

        /*VCF annotations
refseq.inCodingRegion_1
refseq.positionType=intron
refseq.changesAA_1=false
refseq.functionalClass_2=silent;
refseq.positionType_1=CDS;
refseq.spliceDist_2=-41;
refseq.spliceDist_6=109;
refseq.variantAA_4
refseq.variantCodon_
refseq.name2
*/
        String geneName = "";
        String coding = "";
        String functionalClass = "";
        String positionType = "";
        String aaChange = "";


        if (variant.getType() == VariantContext.Type.INDEL) {


            IGVFeature gene = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 1, refSeqGenes);
            if (gene != null) {
                geneName = gene.getName();

                Exon exon = gene.getExonAt(s - 1);
                if (exon == null) {
                    positionType = "intron";
                } else {
                    boolean isUTR = exon.isUTR(s - 1);
                    if (isUTR) {
                        positionType = "utr";
                    } else {
                        positionType = "CDS";
                        coding = "Yes";
                    }
                }
            }

        } else {
            String tmp = variant.getAttributeAsString("refseq.name2", null);
            if (tmp != null) {
                geneName = tmp;
                boolean isCoding = variant.getAttributeAsString("refseq.inCodingRegion", "FALSE").toUpperCase().equals("TRUE");
                if (isCoding) {
                    coding = "Yes";
                }
                functionalClass = variant.getAttributeAsString("refseq.functionalClass", "");
                positionType = variant.getAttributeAsString("refseq.positionType", "");
                String vcf_refAA = variant.getAttributeAsString("refseq.referenceAA", null);
                String vcf_variantAA = variant.getAttributeAsString("refseq.variantAA", null);
                if (vcf_refAA != null && vcf_variantAA != null) {
                    aaChange = vcf_refAA + "->" + vcf_variantAA;
                }

            } else {
                tmp = variant.getAttributeAsString("refseq.name2_1", null);
                if (tmp == null) {
                    // Not in a gene
                } else {
                    geneName = tmp;
                    boolean isCoding = variant.getAttributeAsString("refseq.inCodingRegion_1", "FALSE").toUpperCase().equals("TRUE");
                    if (isCoding) {
                        coding = "Yes";
                    }
                    functionalClass = variant.getAttributeAsString("refseq.functionalClass_1", "");
                    positionType = variant.getAttributeAsString("refseq.positionType_1", "");
                    String vcf_refAA = variant.getAttributeAsString("refseq.referenceAA", null);
                    String vcf_variantAA = variant.getAttributeAsString("refseq.variantAA", null);
                    if (vcf_refAA != null && vcf_variantAA != null) {
                        aaChange = vcf_refAA + "->" + vcf_variantAA;
                    }

                    int i = 2;
                    while (i < 15) {
                        tmp = variant.getAttributeAsString("refseq.name2_" + i, null);
                        if (tmp == null) break;
                        if (geneName.indexOf(tmp) < 0) {
                            if (geneName.length() > 0) geneName += ", ";
                            geneName += tmp;
                        }

                        isCoding = variant.getAttributeAsString("refseq.inCodingRegion_" + i, "FALSE").toUpperCase().equals("TRUE");
                        if (isCoding) {
                            // If any are coding, mark yes
                            coding = "Yes";
                        }

                        tmp = variant.getAttributeAsString("refseq.functionalClass_" + i, null);
                        if (tmp != null && functionalClass.indexOf(tmp) < 0) {
                            if (functionalClass.length() > 0) geneName += ", ";
                            functionalClass += tmp;
                        }

                        tmp = variant.getAttributeAsString("refseq.positionType_" + i, null);
                        if (tmp != null && positionType.indexOf(tmp) < 0) {
                            if (positionType.length() > 0) positionType += ", ";
                            positionType += tmp;

                        }

                        vcf_refAA = variant.getAttributeAsString("refseq.referenceAA_" + i, null);
                        vcf_variantAA = variant.getAttributeAsString("refseq.variantAA_" + i, null);
                        if (vcf_refAA != null && vcf_variantAA != null) {
                            tmp = vcf_refAA + "->" + vcf_variantAA;
                            if (tmp != null && aaChange.indexOf(tmp) < 0) {
                                if (aaChange.length() > 0) aaChange += ", ";
                                aaChange += tmp;
                            }
                        }

                        i++;

                    }
                }
            }
        }


        String ucscGene = "";
        if (ucscGenes != null) {
            IGVFeature ucsc = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 5, ucscGenes);
            if (ucsc != null) ucscGene = ucsc.getName();
        }

        String tfbsName = "";
        if (tfbsFeatures != null) {

            IGVFeature tfbs = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 5, tfbsFeatures);
            if (tfbs != null) tfbsName = tfbs.getName();
        }

        String vntrName = "";
        if (vntrFeatures != null) {
            IGVFeature vntr = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 5, vntrFeatures);
            if (vntr != null) vntrName = vntr.getName();
        }

        String clName = "";
        if (cufflinksFeatures != null) {
            List<Feature> features = FeatureUtils.getAllFeaturesAt(variant.getStart(), 10000, 5, cufflinksFeatures, false);
            if (features != null && features.size() > 0) {
                clName = "intron";
                for (Feature f : features) {
                    IGVFeature cl = (IGVFeature) f;
                    Exon exon = cl.getExonAt(s - 1);
                    if (exon != null) {
                        clName = "exon";
                        break;
                    }

                }
            }
        }

        String scripture = "";
        if (scriptureFeatures != null) {
            List<Feature> features = FeatureUtils.getAllFeaturesAt(variant.getStart(), 10000, 5, scriptureFeatures, false);
            if (features != null && features.size() > 0) {
                scripture = "intron";
                for (Feature f : features) {
                    IGVFeature cl = (IGVFeature) f;
                    Exon exon = cl.getExonAt(s - 1);
                    if (exon != null) {
                        scripture = "exon";
                        break;
                    }

                }
            }
        }

        String trf = "";
        if (trfFeatures != null) {
            IGVFeature cl = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 5, trfFeatures);
            if (cl != null) trf = cl.getName();
        }

        String codingRepeat = "";
        if (codingRepeatFeatures != null) {
            IGVFeature cl = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 5, codingRepeatFeatures);
            if (cl != null) codingRepeat = cl.getName();
        }

        String phastPrimate = "";
        if (phastPrimateFeatures != null) {
            IGVFeature cl = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 2, phastPrimateFeatures);
            if (cl != null) phastPrimate = cl.getName();
        }

        String phastVert = "";
        if (phastVertFeatures != null) {

            IGVFeature cl = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 2, phastVertFeatures);
            if (cl != null) phastVert = cl.getName();
        }

        //
        // String functionalClass = variant.getAttributeAsString("refseq.functionalClass");

        // Calculated annotations
        /*
String coding = "";
String aaChange = ""; //refAA == null ? "" : refAA + "->" + altAA;
String geneName = ""; //variant.getAttributeAsString("refseq.name2", null);
String location = "";
String silent = "";


*/
        String families = "";
        for (String f : affectedFamilies) {
            if (families.length() > 0) {
                families += ",";
            }
            families += f;
        }
        for (String f : record.affectedSingles) {
            if (families.length() > 0) {
                families += ",";
            }
            families += f;
        }
        bedWriter.println(chr + "\t" + s + "\t" + e + "\t" + families);


        /*
  vcf_genename = tmp;
   vcf_coding = variant.getAttributeAsString("refseq.inCodingRegion", null);
   vcf_functionClass = variant.getAttributeAsString("refseq.functionalClass", null);
   vcf_positionType = variant.getAttributeAsString("refseq.positionType", null);
   String vcf_refAA = variant.getAttributeAsString("refseq.referenceAA", "");
   String vcf_variantAA = variant.getAttributeAsString("refseq.variantAA", "");
   vcf_aaChange = vcf_refAA + "->" + vcf_variantAA;

*/
        if (excelWriter != null) {

            String discFamiles = "";
            for (String f : record.discordantFamilies) {
                if (discFamiles.length() > 0) {
                    discFamiles += ",";
                }
                discFamiles += f;
            }

            String type = variant.getType().toString();
            String locusString = chr + ":" + s + "-" + e;
            String id = variant.getID();
            excelWriter.print(type + "\t" + variant.getStart() + "\t" + locusString + "\t" + id + "\t" + families + "\t" + geneName + "\t"
                    + coding + "\t" + functionalClass + "\t" + aaChange + "\t" + positionType +
                    "\t" + ucscGene + "\t" + clName + "\t" + scripture + "\t" +
                    trf + "\t" + codingRepeat + "\t" + vntrName + "\t" + tfbsName + "\t" + phastPrimate + "\t" + phastVert + "\t"
                    + discFamiles + "\t" + record.discordantFamilies.size());


            for (String f : allFamiliesList) {
                excelWriter.print("\t");
                if (affectedFamilies.contains(f)) {
                    excelWriter.print("Yes");
                }

            }
            excelWriter.println();

        }
    }
}
