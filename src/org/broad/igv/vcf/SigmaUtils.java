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

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.feature.AbstractFeatureParser;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.Locus;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.PackedFeatures;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
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


    TribbleFeatureSource source;

    Map<String, List<String>> samples;
    private VCFHeader header;
    private List<Feature> refSeqGenes;

    public static void main(String[] args) throws IOException {

        SigmaUtils utils = new SigmaUtils(args[0]);

        PrintWriter pw = new PrintWriter(new FileWriter("segregating_snps.bed"));
        PrintWriter excelWriter = new PrintWriter(new FileWriter("segregating_snps.xls"));

        excelWriter.println("Locus\tID\tAllele freq\tFamilies\tGene\tCoding ?");


        pw.println("#coord=1");
        pw.println("track name=\"Segregating snps\"");
        utils.walkVariants(pw, excelWriter);

        pw.close();
        excelWriter.close();


    }

    private void initFeatureLists() throws IOException {

        ResourceLocator rl = new ResourceLocator(refSeqURL);
        AsciiLineReader reader = ParsingUtils.openAsciiReader(rl);

        refSeqGenes = AbstractFeatureParser.getInstanceFor(rl).loadFeatures(reader);

        reader.close();

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
            samples.put("S", new ArrayList());
            samples.put("WC", new ArrayList());
            samples.put("CC", new ArrayList());
            samples.put("WC", new ArrayList());
            samples.put("PA", new ArrayList());
            samples.put("Other", new ArrayList());

            for (String s : header.getGenotypeSamples()) {
                if (s.endsWith("-259") || s.endsWith("-265") || s.endsWith("-266") ||
                        s.equals("259") || s.equals("265") || s.equals("266")) {
                    samples.get("AA").add(s);
                } else if (s.endsWith("-701") || s.endsWith("-564") || s.equals("701") || s.equals("564")) {
                    samples.get("BIP").add(s);
                } else if (s.endsWith("-352") || s.endsWith("-414") || s.equals("352") || s.equals("414")
                        || s.endsWith("-375") || s.equals("375")) {
                    samples.get("L").add(s);
                } else if (s.endsWith("-4") || s.endsWith("-8") || s.endsWith("-13") || s.endsWith("-15") ||
                        s.equals("4") || s.equals("8") || s.equals("13") || s.equals("15")) {
                    samples.get("LA").add(s);
                } else if (s.endsWith("-563") || s.endsWith("-566") || s.equals("563") || s.equals("566")) {
                    samples.get("OK").add(s);
                } else if (s.endsWith("-384") || s.endsWith("-391") || s.equals("384") || s.equals("391")) {
                    samples.get("S").add(s);
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


        CloseableTribbleIterator<Feature> features = source.getFeatures(chr, roiStart, roiEnd);

        for (Feature feature : features) {

            VariantContext variant = (VariantContext) feature;
            //char ref = getReference(variant, windowStart, reference);

            // 1 -> 0 based coordinates
            int start = variant.getStart();
            int end = variant.getEnd();

            boolean inNormal = false;
            boolean segregates = false;
            String families = "";
            int familyCount = 0;
            int individualCount = 0;
            for (Map.Entry<String, List<String>> entry : samples.entrySet()) {

                String group = entry.getKey();

                //if (group.equals("S")) {
                //    continue;
                //}

                boolean isHet = true;

                final List<String> sampleList = entry.getValue();
                if (sampleList.size() > 1) {
                    for (String s : sampleList) {
                        Genotype genotype = variant.getGenotype(s);


                        // Check the one un-affected, if it is not homRef this does not segregate
                        if (group.equals("L") && s.endsWith("375")) {
                            if (!genotype.isHomRef()) {
                                inNormal = true;
                                break;
                            }

                        } else {
                            if (!genotype.isHet()) {
                                isHet = false;
                                break;
                            }
                        }
                    }
                    if (inNormal) {
                        break;
                    } else if (isHet) {
                        segregates = true;
                        if (families.length() > 0) {
                            families += ",";
                        }
                        families += group;
                        familyCount++;
                    }
                }
            }

            if (segregates) {
                String chr = variant.getChr();
                int s = variant.getStart();
                int e = variant.getEnd();
                bedWriter.println(chr + "\t" + s + "\t" + e + "\t" + families);

                IGVFeature gene = (IGVFeature) FeatureUtils.getFeatureAt(variant.getStart(), 1, refSeqGenes);
                String geneName = "";
                String coding = "";
                if (gene != null) {
                    geneName = gene.getName();
                    if (gene.getExonAt(s - 1) != null) {
                        coding = "Yes";
                    }
                }

                if (excelWriter != null) {
                    String locusString = chr + ":" + s + "-" + e;
                    String id = variant.getID();
                    String af = (id == null || id.length() == 0 || id.equals(".")) ? "" : variant.getAttributeAsString("AF");
                    excelWriter.println(locusString + "\t" + id + "\t" + af + "\t" + families + "\t" + geneName + "\t" + coding);

                }


            }

        }

    }


}
