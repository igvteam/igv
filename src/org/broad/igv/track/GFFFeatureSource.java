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
import org.broad.igv.util.collections.MultiMap;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * User: jacob
 * Date: 2012-Jun-22
 */
public class GFFFeatureSource extends TribbleFeatureSource {

    private static Logger log = Logger.getLogger(GFFFeatureSource.class);
    public static final String PHASE_STRING = "XXPHASE_STRINGXX";


    public static boolean isGFF(String path) {
        String lowpath = path.toLowerCase();
        if (lowpath.endsWith(".gz")) {
            int idx = lowpath.length() - 3;
            lowpath = lowpath.substring(0, idx);
        }
        if (lowpath.endsWith(".txt")) {
            int idx = lowpath.length() - 4;
            lowpath = lowpath.substring(0, idx);
        }
        return lowpath.endsWith("gff3") || lowpath.endsWith("gvf") || lowpath.endsWith("gff") || lowpath.endsWith("gtf");
    }

    public GFFFeatureSource(String path, Genome genome) throws IOException {
        super(path, genome);
        this.isVCF = false;
    }

    @Override
    public CloseableTribbleIterator<Feature> getFeatures(String chr, int start, int end) throws IOException {

        CloseableTribbleIterator<Feature> rawIter = super.getFeatures(chr, start, end);
        GFFCombiner combiner = new GFFCombiner();
        combiner.addFeatures(rawIter);

        return new WrappedIterator(combiner.combineFeatures().iterator());
    }

    /**
     * The GFF spec is available at http://www.sequenceontology.org/gff3.shtml
     * <p/>
     * GFF Combiner is needed because certain features belong to parent features
     * See GFFCodec for some details there. We want to combine these, and make sure
     * that features are modified properly. mRNAs should contain exons, and be modified
     * by 5'/3'-UTR, and coding start sites, etc.
     * <p/>
     * The GFF spec seems to indicate that multiple lines of type "CDS" should be treated
     * as the same feature iff they have the same ID attribute. Example from spec:
     * <p/>
     * <p/>
     * <br/>
     * chr6	.	CDS	3301	3902	.	+	0	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
     * chr6	.	CDS	5000	5500	.	+	1	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
     * chr6	.	CDS	7000	7600	.	+	1	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
     * chr6	.	CDS	3391	3902	.	+	0	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
     * chr6	.	CDS	5000	5500	.	+	1	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
     * chr6	.	CDS	7000	7600	.	+	1	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
     * <br/>
     * <p/>
     * This should represent 2 separate protein coding sequences for mRNA00003: cds00003 and cds00004
     * This would be fine, except data we see in the wild often assigns a unique ID to each CDS
     * feature. So the above data might look like:
     * <br/>
     * chr6	.	CDS	3301	3902	.	+	0	ID=cds00003.1;Parent=mRNA00003;Name=edenprotein.3
     * chr6	.	CDS	5000	5500	.	+	1	ID=cds00003.2;Parent=mRNA00003;Name=edenprotein.3
     * chr6	.	CDS	7000	7600	.	+	1	ID=cds00003.3;Parent=mRNA00003;Name=edenprotein.3
     * chr6	.	CDS	3391	3902	.	+	0	ID=cds00004.1;Parent=mRNA00003;Name=edenprotein.4
     * chr6	.	CDS	5000	5500	.	+	1	ID=cds00004.2;Parent=mRNA00003;Name=edenprotein.4
     * chr6	.	CDS	7000	7600	.	+	1	ID=cds00004.3;Parent=mRNA00003;Name=edenprotein.4
     * <br/>
     * Can't say as I blame the users here, having unique IDs is usually a good thing.
     * Which leaves us in a bit of a conundrum. How do we tell the difference between two protein coding
     * sequences, or members of the same coding sequence?
     * <p/>
     * We can say for sure that multiple CDS features with the same ID value are part of the same
     * coding sequence. That's the case where the file (gasp!) conforms to the spec
     * <p/>
     * With multiple ID values, we could look at spatial information (start/end coords) and perform
     * some type of packing procedure. Since we can't assume any naming convention on the IDs that's about
     * the best we could do. But that's likely to show misleading data in subtle ways
     * <p/>
     * So here's the plan: On a per parent basis, we look at the CDS elements. There should never be more
     * coding sequences created than unique ID values.
     * <br/>
     * <p/>
     * 0. If each one has a unique ID value, we assign
     * them all to the same coding sequence ({@link org.broad.tribble.Feature}).
     * <p/>
     * <br/>
     * 1. If some sequences have the same ID, assume that features have an ID
     * per coding sequence. Create a {@link org.broad.tribble.Feature} for each one.
     * <p/>
     * This may not always be correct if people mix and match unique IDs with grouping IDs.
     */
    public static class GFFCombiner {

        Map<String, GFF3Transcript> transcriptCache;
        Map<String, BasicFeature> geneCache;
        Queue<BasicFeature> subFeatures;

        public GFFCombiner() {
            int numElements = 50000;
            transcriptCache = new HashMap(numElements / 2);
            geneCache = new HashMap(numElements / 10);
            subFeatures = new ArrayDeque<BasicFeature>(numElements);
        }

        /**
         * First pass, create transcripts so everything
         * has a parent (if it exists in the iterator)
         *
         * @param rawIter
         */
        public void addFeatures(Iterator<Feature> rawIter) {
            while (rawIter.hasNext()) {
                addFeature((BasicFeature) rawIter.next());
            }
        }


        public void addFeature(BasicFeature bf) {
            String featureType = bf.getType();
            if (featureType.equalsIgnoreCase("CDS_parts") || featureType.equalsIgnoreCase("intron")
                    || GFFCodec.exonTerms.contains(featureType)) {
                //Save for later, add parents first
                subFeatures.add(bf);
            } else {
                //Can't guarantee it's a parent, but better safe than sorry
                //Each of these gets its own transcript
                createTranscript(bf);
            }
        }

        public List<Feature> combineFeatures() {

            BasicFeature bf;
            while ((bf = subFeatures.poll()) != null) {
                String featureType = bf.getType();

                if (featureType.equalsIgnoreCase("CDS_parts") || featureType.equalsIgnoreCase("intron")) {
                    String[] parentIds = bf.getParentIds();
                    for (String pid : parentIds) {
                        getGFF3Transcript(pid).addCDSParts(bf.getChr(), bf.getStart(), bf.getEnd());
                    }
                } else if (GFFCodec.exonTerms.contains(featureType)) {
                    incorporateExon(bf);
                }
            }

            // Create and add IGV genes
            List<Feature> finalFeatures = new ArrayList<Feature>(transcriptCache.size());
            Iterator<GFF3Transcript> iter = transcriptCache.values().iterator();
            while (iter.hasNext()) {
                List<? extends Feature> igvTranscript = iter.next().createTranscript();
                //To save memory, remove transcript from cache once built
                if (igvTranscript != null) {
                    finalFeatures.addAll(igvTranscript);
                }
                //To save memeroy, remove GFF3Transcript directly after builder
                iter.remove();
            }

            FeatureUtils.sortFeatureList(finalFeatures);

            return finalFeatures;
        }

        /**
         * Combine all features in the iterator. Match elements to their parents,
         * including sorting out exons/CDS/UTR business.
         *
         * @return
         */
//        public List<Feature> combineFeatures_old(Iterator<Feature> rawIter) {
//            List<Feature> features = new ArrayList();
//
//            while(rawIter.hasNext()) {
//                BasicFeature bf = (BasicFeature) rawIter.next();
//                String featureType = bf.getType();
//
//                if (featureType.equalsIgnoreCase("CDS_parts") || featureType.equalsIgnoreCase("intron")) {
//                    String[] parentIds = bf.getParentIds();
//                    for (String pid : parentIds) {
//                        getGFF3Transcript(pid).addCDSParts(bf.getChr(), bf.getStart(), bf.getEnd());
//                    }
//                } else if (GFFCodec.exonTerms.contains(featureType)) {
//                    incorporateExon(bf);
//                } else {
//                    createTranscript(bf);
//
//                }
//            }
//
//            // Create and add IGV genes
//            for (GFF3Transcript transcript : transcriptCache.values()) {
//                List<? extends Feature> igvTranscript = transcript.createTranscript();
//                if (igvTranscript != null) {
//                    features.addAll(igvTranscript);
//                }
//            }
//
//            FeatureUtils.sortFeatureList(features);
//
//            return features;
//        }
        private void incorporateExon(BasicFeature bf) {
            String featureType = bf.getType();
            String[] parentIds = bf.getParentIds();

            //If the exon has no parent, we effectively create a new feature
            //and treat it as it's own parent.
            if (!isValidParentIds(parentIds) && bf.getIdentifier() != null) {
                parentIds = new String[]{bf.getIdentifier()};
                bf.setParentIds(parentIds);
                createGFF3Transcript(bf.getIdentifier(), bf, parentIds[0]);
            }

            // Make a copy of the exon record for each parent
            for (String pid : parentIds) {

                Exon exon = new Exon(bf.getChr(), bf.getStart(), bf.getEnd(), bf.getStrand());
                exon.setPhase(bf.getReadingFrame());

                if(bf.getColor() != null) exon.setColor(bf.getColor());

                //String sPhase = bf.getAttributes().remove(PHASE_STRING);
                exon.setAttributes(bf.getAttributes());
                exon.setUTR(GFFCodec.utrTerms.contains(featureType));


//                if (sPhase != null) {
//                    int phase = parsePhase(sPhase);
//                    if (phase >= 0) {
//                        exon.setPhase(phase);
//
//                    }
//                }

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

        private void createTranscript(BasicFeature bf) {

            String featureType = bf.getType();
            String id = bf.getIdentifier();

            if (featureType.equalsIgnoreCase("gene")) {
                geneCache.put(id, bf);
            }

            String pid = null;
            String[] parentIds = bf.getParentIds();
            if (isValidParentIds(parentIds)) {
                pid = parentIds[0];
            }
            //Transcripts get turned into features at end
            createGFF3Transcript(id, bf, pid);
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
//            return transcriptCache.get(id);
            GFF3Transcript transcript = transcriptCache.get(id);
            if (transcript == null) {
                transcript = new GFF3Transcript(id, geneCache);
                transcriptCache.put(id, transcript);
            }
            return transcript;
        }

        private GFF3Transcript createGFF3Transcript(String id, BasicFeature mRNA, String parentId) {
            GFF3Transcript transcript = new GFF3Transcript(id, geneCache);
            transcript.setBasicTranscript(mRNA, parentId);
            transcriptCache.put(id, transcript);
            return transcript;
        }

        private boolean isValidParentIds(String[] parentIds) {
            return parentIds != null && parentIds.length > 0 && parentIds[0] != null &&
                    parentIds[0].trim().length() > 0 && !parentIds[0].equals(".");
        }

    }

    static class GFF3Transcript {


        /**
         * We store certain objects under a map, by their ID.
         * This is the special key used when no ID field is present.
         * Normally, all objects which have the same ID represent the same feature.
         * However, that is not true
         */
        private static final String ID_ABSENT = "SPECIALIDFOROBJECTSWITHNOID";

        private String id;
        private Set<Exon> exons = new HashSet<Exon>();
        private List<Exon> cdss = new ArrayList<Exon>();

        private Exon fivePrimeUTR;
        private Exon threePrimeUTR;
        private BasicFeature basicTranscript;
        private String parentId;
        String chr = null;
        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        String description;
        MultiMap<String, String> attributes;
        Map<String, BasicFeature> geneCache;

        GFF3Transcript(String id, Map<String, BasicFeature> geneCache) {
            this.id = id;
            this.geneCache = geneCache;
        }

        void setBasicTranscript(BasicFeature basicTranscript, String parent) {
            this.basicTranscript = basicTranscript;
            this.parentId = parent;
            if (basicTranscript.getName() == null) {
                basicTranscript.setName(basicTranscript.getIdentifier());
            }

            if (basicTranscript.getName() == null) return;

            int prefixIndex = basicTranscript.getName().indexOf(":");
            if (prefixIndex > 0) {
                basicTranscript.setName(basicTranscript.getName().substring(prefixIndex + 1));
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
         * Due to CDS configuration, this can be multiple features.
         * One for each protein coding sequence.
         *
         * @return
         */
        List<? extends Feature> createTranscript() {

            /**
             * We create a separate set of exons for each CDS id
             * This is to handle alternate transcripts.
             */
            Map<String, List<Exon>> cdsByID = new HashMap<String, List<Exon>>();
            for (Exon cds : cdss) {
                String cdsID = cds.getAttributes().get("ID");
                List<Exon> cdsArr = cdsByID.get(cdsID);
                if (cdsArr == null) {
                    cdsArr = new ArrayList<Exon>();
                    cdsByID.put(cdsID, cdsArr);
                }
                cdsArr.add(cds);
            }

            if (cdsByID.size() == cdss.size()) {
                //Fairly easy case, all CDSs have unique IDs
                //We have only one protein coding sequence,
                //made up of this sequence

                cdsByID = new HashMap<String, List<Exon>>(1);
                cdsByID.put(this.id, cdss);
            } else {
                //Assume each unique key belongs to its own CDS
                //Which is what we built already
            }

            //At this point, cdsByID should be a map from transcript id -> list of coding sequences


            //We modify some exons live, so we keep track of which ones
            //and make a copy if necessary
            Set<Exon> foundExons = new HashSet<Exon>();

            /**
             * Loop through each CDS, and modify exons as necessary to include coding start/stop
             * information.
             */
            List<BasicFeature> transcriptList = new ArrayList<BasicFeature>(cdsByID.size());
            for (String cdsID : cdsByID.keySet()) {

                Set<Exon> exonSet = new HashSet<Exon>();
                for (Exon cds : cdsByID.get(cdsID)) {

                    Exon exon = findMatchingExon(exons, cds);
                    if (exon == null) {
                        cds.setCodingStart(cds.getStart());
                        cds.setCodingEnd(cds.getEnd());
                        exonSet.add(cds);
                    } else {
                        if (foundExons.contains(exon)) {
                            exon = exon.copy();
                        } else {
                            foundExons.add(exon);
                        }
                        exon.setCodingStart(cds.getStart());
                        exon.setCodingEnd(cds.getEnd());
                        exon.setReadingFrame(cds.getReadingShift());
                        exonSet.add(exon);
                    }
                }

                /**
                 * If there are any exon features in the file which don't have
                 * a corresponding CDS, we simply add them wholesale. The most likely case
                 * is that they are 5'/3' UTR and will be modified later accordingly.
                 */
                for(Exon exon: exons){
                    if(!foundExons.contains(exon)){
                        exonSet.add(exon);
                    }
                }


                Strand strand = Strand.NONE;
                String name = null;
                String curChr = chr;
                int curStart = start;
                int curEnd = end;
                for (Exon exon : exonSet) {
                    curChr = exon.getChr();
                    strand = exon.getStrand();
                    curStart = Math.min(exon.getStart(), curStart);
                    curEnd = Math.max(exon.getEnd(), curEnd);
                    name = exon.getName();
                }


                BasicFeature curTranscript = basicTranscript;
                if (curTranscript == null) {
                    curTranscript = new BasicFeature(curChr, curStart, curEnd, strand);
                    curTranscript.setIdentifier(id);
                    curTranscript.setName(name == null ? id : name);
                    curTranscript.setDescription(description);
                    curTranscript.setAttributes(attributes);
                } else if (cdsByID.size() >= 2) {
                    curTranscript = basicTranscript.copy();
                }

                BasicFeature gene = null;
                if ((parentId != null) && geneCache.containsKey(parentId) && gene == null) {
                    gene = geneCache.get(parentId);
                    geneCache.remove(parentId);
                }

                if (gene != null) {
                    if (curTranscript.getName() == null && gene.getName() != null) {
                        curTranscript.setName(gene.getName());
                    }
                    curTranscript.setDescription("Transcript<br>" + curTranscript.getDescription() + "<br>--------<br>Gene<br>" + gene.getDescription());
                }

                for (Exon exon : exonSet) {
                    curTranscript.addExon(exon);
                }

                curTranscript.sortExons();

                // If 5'UTR is represented by an exon, adjust its start
                if (fivePrimeUTR != null) {
                    adjustBoundariesByUTR(curTranscript, fivePrimeUTR, false);
                }

                if (threePrimeUTR != null) {
                    adjustBoundariesByUTR(curTranscript, threePrimeUTR, true);
                }

                transcriptList.add(curTranscript);
            }

            return transcriptList;
        }

        private void adjustBoundariesByUTR(BasicFeature transcript, Exon UTR, boolean threePrime) {
            UTR.setUTR(true);
            transcript.addExon(UTR);
            Exon exon = findMatchingExon(transcript.getExons(), UTR);
            if (exon != null) {
                //This UTR either truncates the end or the beginning.
                //Being on the negative strand or 3'-UTR would truncate the end
                //both is a double negative
                boolean truncStart = (exon.getStrand() == Strand.POSITIVE) ^ threePrime;
                if (truncStart) {
                    exon.setStart(UTR.getEnd());
                } else {
                    exon.setEnd(UTR.getStart());
                }
            }
        }

        Exon findMatchingExon(Iterable<Exon> exons, IGVFeature cds) {
            for (Exon exon : exons) {
                if (exon.contains(cds)) {
                    return exon;
                }
            }
            return null;
        }
    }
}
