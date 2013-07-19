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
import org.broad.igv.feature.tribble.IGVFeatureReader;
import org.broad.igv.util.collections.MultiMap;
import org.broad.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * User: jacob
 * Date: 2012-Jun-22
 */
public class GFFFeatureSource extends TribbleFeatureSource {

    private static Logger log = Logger.getLogger(GFFFeatureSource.class);

    static Set<String> fivePrimeUTRTypes = new HashSet<String>();
    static Set<String> threePrimeUTRTypes = new HashSet<String>();
    static Set<String> utrTypes = new HashSet<String>();
    static Set<String> cdsTypes = new HashSet<String>();
    static Set<String> transcriptTypes = new HashSet<String>();
    static Set<String> transcriptPartTypes = new HashSet<String>();

    static {
        fivePrimeUTRTypes.add("five_prime_UTR");
        fivePrimeUTRTypes.add("5'-UTR");
        fivePrimeUTRTypes.add("5'-utr");
        fivePrimeUTRTypes.add("5UTR");
        threePrimeUTRTypes.add("three_prime_UTR");
        threePrimeUTRTypes.add("3'-utr");
        threePrimeUTRTypes.add("3'-UTR");
        threePrimeUTRTypes.add("3UTR");
        utrTypes.addAll(fivePrimeUTRTypes);
        utrTypes.addAll(threePrimeUTRTypes);
        cdsTypes.add("CDS");
        cdsTypes.add("cds");
        transcriptPartTypes.addAll(utrTypes);
        transcriptPartTypes.addAll(cdsTypes);

        transcriptTypes.add("transcript");
        transcriptTypes.add("mature_transcript");
        transcriptTypes.add("processed_transcript");
        transcriptTypes.add("mrna");
        transcriptTypes.add("mRNA");
    }


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
    public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {

        Iterator<Feature> rawIter = super.getFeatures(chr, start, end);
        GFFCombiner combiner = (new GFFCombiner()).addFeatures(rawIter);


        return new WrappedIterator(combiner.combineFeatures().iterator());
    }

    /**
     * The GFF spec is available at http://www.sequenceontology.org/gff3.shtml
     * <p/>
     * GFF Combiner is needed because IGV represents a transcript (e.g. an mRNAs) as a single feature with,
     * optionally, a collection of child exons.  The feature can have a "thick" start and end, corresponding
     * to coding start and end.  This representation comes from a literal transcript of certain UCSC formats.
     * The same feature in a GFF file, on the other hand, is represented as a graph of sub features with
     * defined by parent-child relations.  GFF3 formalizes the feature types and their relationships to some
     * degree in the Sequence Ontology,
     */
    public static class GFFCombiner {


        Map<String, GFFTranscript> gffTranscripts;
        List<Feature> igvFeatures;

        public GFFCombiner() {
            int numElements = 50000;
            gffTranscripts = new HashMap(numElements / 2);
            igvFeatures = new ArrayList(numElements);
        }

        /**
         * First pass, create transcripts so everything
         * has a parent (if it exists in the iterator)
         *
         * @param rawIter
         * @return this, for chaining
         */
        public GFFCombiner addFeatures(Iterator<Feature> rawIter) {
            while (rawIter.hasNext()) {
                addFeature((BasicFeature) rawIter.next());
            }
            return this;
        }


        public void addFeature(BasicFeature bf) {

            String featureType = bf.getType();

            if (transcriptTypes.contains(featureType)) {
                GFFTranscript gffTranscript = getGFFTranscript(bf.getIdentifier());
                gffTranscript.setTranscript(bf);
            } else if (transcriptPartTypes.contains(featureType)){
                for (String parentID : bf.getParentIds()) {
                    GFFTranscript gffTranscript = getGFFTranscript(parentID);
                    gffTranscript.addSubpart(bf);
                }
            }
            else {
                igvFeatures.add(bf);
            }
        }


        private GFFTranscript getGFFTranscript(String id) {
//            return transcriptCache.get(id);
            GFFTranscript transcript = gffTranscripts.get(id);
            if (transcript == null) {
                transcript = new GFFTranscript(id);
                gffTranscripts.put(id, transcript);
            }
            return transcript;
        }


        public List<Feature> combineFeatures() {
            for (GFFTranscript gffTranscript : gffTranscripts.values()) {
                igvFeatures.add(gffTranscript.createIGVFeature());
            }
            FeatureUtils.sortFeatureList(igvFeatures);
            return igvFeatures;
        }
    }

    private static class GFFTranscript {



        static String[] nameFields = {"Name", "name", "Alias", "gene", "primary_name", "locus", "alias", "systematic_id", "ID", "transcript_id"};


        /**
         * We store certain objects under a map, by their ID.
         * This is the special key used when no ID field is present.
         * Normally, all objects which have the same ID represent the same feature.
         * However, that is not true
         */
        private static final String ID_ABSENT = "SPECIALIDFOROBJECTSWITHNOID";

        private String id;
        private String name;
        private Set<BasicFeature> soExons = new HashSet<BasicFeature>();
        private List<BasicFeature> subFeatures = new ArrayList<BasicFeature>();

        String chr = null;
        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        Strand strand;
        String description;
        MultiMap<String, String> attributes;

        GFFTranscript(String id) {
            this.id = id;
        }

        void setTranscript(BasicFeature transcript) {

            chr = transcript.getChr();
            start = transcript.getStart();
            end = transcript.getEnd();
            strand = transcript.getStrand();
            final MultiMap<String, String> transcriptAttributes = transcript.getAttributes();
            if (attributes == null) {
                attributes = transcriptAttributes;
            } else {
                for (String key : transcriptAttributes.keys()) {
                    attributes.put(key, transcriptAttributes.get(key));
                }
            }

            if (transcript.getName() != null) {
                name = transcript.getName();
            } else {
                for (String nm : nameFields) {
                    if (attributes.containsKey(nm)) {
                        name = attributes.get(nm);
                        break;
                    }
                }
            }

        }


        public void addSubpart(BasicFeature bf) {
            this.chr = bf.getChr();
            this.start = Math.min(bf.getStart(), start);
            this.end = Math.max(bf.getEnd(), end);
            this.strand = bf.getStrand();
            String type = bf.getType();

            if (type.equals("exon")) {
                soExons.add(bf);
            } else if (type.equals("intron")) {
                // IGV does not explicitly represent introns
            } else {
                subFeatures.add(bf);
            }


        }


        BasicFeature createIGVFeature() {

            // First add explicit 'exon' types
            List<Exon> exons = new ArrayList<Exon>();
            for (BasicFeature bf : soExons) {
                Exon exon = new Exon(bf.getChr(), bf.getStart(), bf.getEnd(), bf.getStrand());
                exon.setAttributes(bf.getAttributes());
                exons.add(exon);
            }

            // Now create exons for UTRs and CDS types, or modify existing exons
            for (BasicFeature bf : subFeatures) {
                Exon exon = findMatchingExon(exons, bf);
                if (exon == null) {
                    exon = new Exon(bf.getChr(), bf.getStart(), bf.getEnd(), bf.getStrand());
                    exon.setAttributes(bf.getAttributes());
                } else {
                    exon.getAttributes().addAll(bf.getAttributes());
                }
                adjustCodingPositions(exon, bf);
                exons.add(exon);
            }

            FeatureUtils.sortFeatureList(exons);

            BasicFeature bf = new BasicFeature(chr, start, end, strand);
            bf.setName(name);
            bf.setIdentifier(id);
            bf.setAttributes(attributes);
            for (Exon exon : exons) {
                bf.addExon(exon);
            }
            return bf;
        }


        private void adjustCodingPositions(Exon exon, BasicFeature sf) {

            final String type = sf.getType();
            if (cdsTypes.contains(type)) {
                // Don't do this for "CDS_parts" for now.  I'm not sure what this type really means
                if (!type.equals("CDS_parts")) {
                    exon.setCodingStart(sf.getStart());
                    exon.setCodingEnd(sf.getEnd());
                }
            } else if (utrTypes.contains(type)) {
                exon.setUTR(true);
                boolean rhs =
                        (fivePrimeUTRTypes.contains(type) && sf.getStrand() == Strand.POSITIVE) ||
                                (threePrimeUTRTypes.contains(type) && sf.getStrand() == Strand.NEGATIVE);
                if (rhs) {
                    exon.setCodingStart(sf.getEnd());
                } else {
                    exon.setCodingEnd(sf.getStart());
                }
            }
        }

        Exon findMatchingExon(List<Exon> exons, BasicFeature bf) {
            for (Exon exon : exons) {
                if (exon.contains(bf)) {
                    return exon;
                }
            }
            return null;
        }

    }
}
