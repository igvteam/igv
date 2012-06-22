/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * User: jacob
 * Date: 2012-Jun-22
 */
public class GFFFeatureSource extends TribbleFeatureSource{

    private static Logger log = Logger.getLogger(GFFFeatureSource.class);
    public static final String PHASE_STRING = "XXPHASE_STRINGXX";

    public GFFFeatureSource(String path, Genome genome) throws IOException{
        super(path, genome);
        this.isVCF = false;
    }

    @Override
    public CloseableTribbleIterator<Feature> getFeatures(String chr, int start, int end) throws IOException {
        CloseableTribbleIterator<Feature> rawIter = super.getFeatures(chr, start, end);

        return new WrappedIterator((new GFFCombiner()).combineFeatures(rawIter).iterator());
    }

    public static class GFFCombiner{

        Map<String, GFF3Transcript> transcriptCache = new HashMap(50000);
        Map<String, BasicFeature> geneCache = new HashMap(50000);


        private void clearCombinedFeatures() {
            transcriptCache.clear();
            geneCache.clear();
        }

        public List<Feature> combineFeatures(Iterator<Feature> rawIter) {
            List<org.broad.tribble.Feature> features = new ArrayList();

            while (rawIter.hasNext()) {
                BasicFeature bf = (BasicFeature) rawIter.next();
                String featureType = bf.getType();

                if (featureType.equalsIgnoreCase("CDS_parts") || featureType.equalsIgnoreCase("intron")) {
                    String[] parentIds = bf.getParentIds();
                    for (String pid : parentIds) {
                        getGFF3Transcript(pid).addCDSParts(bf.getChr(), bf.getStart(), bf.getEnd());
                    }
                } else if (GFFCodec.exonTerms.contains(featureType)) {
                    incorporateExon(bf);
                } else {
                    Feature f = incorporateFeature(bf);
                    if (f != null) {
                        features.add(f);
                    }
                }
            }

            // Create and add IGV genes
            for (GFF3Transcript transcript : transcriptCache.values()) {
                Feature igvTranscript = transcript.createTranscript();
                if (igvTranscript != null) {
                    features.add(igvTranscript);
                }
            }

            FeatureUtils.sortFeatureList(features);

            return features;
        }

        private void incorporateExon(BasicFeature bf) {
            String featureType = bf.getType();
            String[] parentIds = bf.getParentIds();
            // Make a copy of the exon record for each parent
            for (String pid : parentIds) {

                Exon exon = new Exon(bf.getChr(), bf.getStart(), bf.getEnd(), bf.getStrand());

                String sPhase = bf.getAttributes().remove(PHASE_STRING);
                exon.setAttributes(bf.getAttributes());

                exon.setUTR(GFFCodec.utrTerms.contains(featureType));


                int phase = parsePhase(sPhase);
                if (phase >= 0) {
                    exon.setPhase(phase);

                }
                exon.setName(bf.getName());

                if (featureType.equalsIgnoreCase("exon")) {
                    getGFF3Transcript(pid).addExon(exon);
                } else if (featureType.equalsIgnoreCase("CDS")) {
                    getGFF3Transcript(pid).addCDS(exon);
                } else if (featureType.equalsIgnoreCase("five_prime_UTR") || featureType.equalsIgnoreCase("5'-UTR")) {
                    getGFF3Transcript(pid).setFivePrimeUTR(exon);
                } else if (featureType.equalsIgnoreCase("three_prime_UTR") || featureType.equalsIgnoreCase("3'-UTR")) {
                    getGFF3Transcript(pid).setThreePrimeUTR(exon);
                }
            }
        }

        private BasicFeature incorporateFeature(BasicFeature bf) {

            String featureType = bf.getType();
            String id = bf.getIdentifier();

            if (featureType.equalsIgnoreCase("gene")) {
                geneCache.put(id, bf);
                return bf;
            } else if (featureType.equalsIgnoreCase("mRNA") || featureType.equalsIgnoreCase("transcript")) {
                String pid = null;
                String[] parentIds = bf.getParentIds();
                if (parentIds != null && parentIds.length > 0) {
                    pid = parentIds[0];
                }
                //Transcripts get turned into features at end
                getGFF3Transcript(id).transcript(bf, pid);
            } else {
                return bf;
            }
            return null;
        }

        private int parsePhase(String phaseString) {
            int phase = -1;
            if (!phaseString.equals(".")) {
                try {
                    phase = Integer.parseInt(phaseString);
                } catch (NumberFormatException numberFormatException) {
                    // Just skip setting the phase
                    log.error("GFF3 error: non numeric phase: " + phaseString);
                }
            }
            return phase;
        }


        private GFF3Transcript getGFF3Transcript(String id) {
            GFF3Transcript transcript = transcriptCache.get(id);
            if (transcript == null) {
                transcript = new GFF3Transcript(id, geneCache);
                transcriptCache.put(id, transcript);
            }
            return transcript;
        }
    }


    static class GFF3Transcript {

        private String id;
        private Set<Exon> exons = new HashSet();
        private List<Exon> cdss = new ArrayList();
        private Exon fivePrimeUTR;
        private Exon threePrimeUTR;
        private BasicFeature transcript;
        private String parentId;
        String chr = null;
        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        String description;
        Map<String, String> attributes;
        Map<String, BasicFeature> geneCache;

        GFF3Transcript(String id, Map<String, BasicFeature> geneCache) {
            this.id = id;
            this.geneCache = geneCache;
        }

        void transcript(BasicFeature mRNA, String parent) {
            this.transcript = mRNA;
            this.parentId = parent;
            if (mRNA.getName() == null) {
                mRNA.setName(mRNA.getIdentifier());
            }
            int prefixIndex = mRNA.getName().indexOf(":");
            if (prefixIndex > 0) {
                mRNA.setName(mRNA.getName().substring(prefixIndex + 1));
            }

        }

        void setFivePrimeUTR(Exon exon) {
            fivePrimeUTR = exon;
            this.start = Math.min(exon.getStart(), start);
            this.end = Math.max(exon.getEnd(), end);
        }

        void setThreePrimeUTR(Exon exon) {
            threePrimeUTR = exon;
            this.start = Math.min(exon.getStart(), start);
            this.end = Math.max(exon.getEnd(), end);
        }

        void addExon(Exon exon) {
            exons.add(exon);
            this.start = Math.min(exon.getStart(), start);
            this.end = Math.max(exon.getEnd(), end);
        }

        void addCDS(Exon cds) {
            cdss.add(cds);
            this.start = Math.min(cds.getStart(), start);
            this.end = Math.max(cds.getEnd(), end);
        }

        void addCDSParts(String chr, int start, int end) {
            this.chr = chr;
            this.start = Math.min(this.start, start);
            this.end = Math.max(this.end, end);
        }

        /**
         * Create a transcript from its constituent parts.
         *
         * @return
         */
        Feature createTranscript() {

            Strand strand = Strand.NONE;
            String name = null;

            // Combine CDS and exons
            while (!cdss.isEmpty()) {
                Exon cds = cdss.get(0);
                Exon exon = findMatchingExon(cds);
                if (exon == null) {
                    cds.setCodingStart(cds.getStart());
                    cds.setCodingEnd(cds.getEnd());
                    exons.add(cds);
                } else {
                    exon.setCodingStart(cds.getStart());
                    exon.setCodingEnd(cds.getEnd());
                    exon.setReadingFrame(cds.getReadingShift());
                }
                cdss.remove(0);
            }
            for (Exon exon : exons) {
                chr = exon.getChr();
                strand = exon.getStrand();
                start = Math.min(exon.getStart(), start);
                end = Math.max(exon.getEnd(), end);
                name = exon.getName();
            }


            if (transcript == null) {
                transcript = new BasicFeature(chr, start, end, strand);
                transcript.setIdentifier(id);
                transcript.setName(name == null ? id : name);
                transcript.setDescription(description);
                transcript.setAttributes(attributes);
            }

            if ((parentId != null) && geneCache.containsKey(parentId)) {
                BasicFeature gene = geneCache.get(parentId);
                geneCache.remove(parentId);
                if (transcript.getName() == null && gene.getName() != null) {
                    transcript.setName(gene.getName());
                }
                transcript.setDescription("Transcript<br>" + transcript.getDescription() + "<br>--------<br>Gene<br>" + gene.getDescription());
            }

            for (Exon exon : exons) {
                transcript.addExon(exon);
            }


            transcript.sortExons();

            // If 5'UTR is represented by an exon, adjust its start, else add an exon to represent 5'utr
            if (fivePrimeUTR != null) {
                fivePrimeUTR.setUTR(true);
                transcript.addExon(fivePrimeUTR);
                Exon exon = findMatchingExon(fivePrimeUTR);
                if (exon != null) {
                    if (exon.getStrand() == Strand.POSITIVE) {
                        exon.setStart(fivePrimeUTR.getEnd());
                    } else {
                        exon.setEnd(fivePrimeUTR.getStart());
                    }
                }
            }

            if (threePrimeUTR != null) {
                threePrimeUTR.setUTR(true);
                transcript.addExon(threePrimeUTR);
                Exon exon = findMatchingExon(threePrimeUTR);
                if (exon != null) {
                    if (exon.getStrand() == Strand.POSITIVE) {
                        exon.setEnd(threePrimeUTR.getStart());
                    } else {
                        exon.setStart(threePrimeUTR.getEnd());
                    }
                }
            }

            return transcript;
        }

        Exon findMatchingExon(IGVFeature cds) {
            for (Exon exon : exons) {
                if (exon.contains(cds)) {
                    return exon;
                }
            }
            return null;
        }
    }
}
