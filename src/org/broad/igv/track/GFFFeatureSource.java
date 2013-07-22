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
import org.broad.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * User: jacob
 * Date: 2012-Jun-22
 */
public class GFFFeatureSource extends TribbleFeatureSource {

    private static Logger log = Logger.getLogger(GFFFeatureSource.class);


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

        List<Feature> igvFeatures;
        Map<String, BasicFeature> gffFeatures;
        List<BasicFeature> gffExons;
        Map<String, GFFCdsCltn> gffCdss;
        List<BasicFeature> gffUtrs;

        public GFFCombiner() {
            int numElements = 10000;
            igvFeatures = new ArrayList<Feature>(numElements);
            gffFeatures = new HashMap<String, BasicFeature>(numElements);
            gffExons = new ArrayList<BasicFeature>(numElements);
            gffCdss = new LinkedHashMap<String, GFFCdsCltn>(numElements);
            gffUtrs = new ArrayList<BasicFeature>(numElements);

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
            String[] parentIDs = bf.getParentIds();
            String id = bf.getIdentifier();
            if (SequenceOntology.exonTypes.contains(featureType) && parentIDs != null) {
                gffExons.add(bf);
            } else if (SequenceOntology.utrTypes.contains(featureType) && parentIDs != null) {
                gffUtrs.add(bf);
            } else if (SequenceOntology.cdsTypes.contains(featureType) && parentIDs != null) {
                for (String pid : parentIDs) {
                    GFFCdsCltn cds = gffCdss.get(pid);
                    if (cds == null) {
                        cds = new GFFCdsCltn(pid);
                        gffCdss.put(pid, cds);
                    }
                    cds.addPart(bf);
                }
            } else if (id != null) {
                gffFeatures.put(id, bf);
            } else {
                igvFeatures.add(bf); // Just use this feature  as is.
            }
        }


        public List<Feature> combineFeatures() {


            for (BasicFeature gffExon : gffExons) {
                final String[] parentIds = gffExon.getParentIds();
                for (String parentId : parentIds) {
                    BasicFeature parent = gffFeatures.get(parentId);
                    if (parent == null) {
                        igvFeatures.add(gffExon);
                    } else {
                        final Exon exon = new Exon(gffExon);
                        exon.setNonCoding(!SequenceOntology.isCoding(gffExon.getType()));
                        parent.addExon(exon);
                    }
                }
            }

            // Now process utrs.  Modify exon if its already defined, create a new one if not.
            for (BasicFeature utr : gffUtrs) {
                for (String parentId : utr.getParentIds()) {
                    BasicFeature parent = gffFeatures.get(parentId);
                    if (parent == null) {
                        igvFeatures.add(utr);
                    } else {
                        parent.addUTRorCDS(utr);
                    }
                }
            }

            // Finally overlay cdss
            for (GFFCdsCltn gffCdsCltn : gffCdss.values()) {

                // Get the parent.
                String parentId = gffCdsCltn.getParentId();
                BasicFeature parent = gffFeatures.get(parentId);
                if (parent == null) {
                    // Create a "dummy" transcript for the orphaned cds records
                    parent = new BasicFeature(gffCdsCltn.chr, gffCdsCltn.start, gffCdsCltn.end, gffCdsCltn.strand);
                    igvFeatures.add(parent);
                }

                // Now add the cds objects.  There are 2 conventions in use for describing the coding section of mRNAs
                // (1) All cds records for the same isoform get the same id.  CDS objects with different ids then
                // imply different isoforms.  In IGV we need to create a parent object for each.  (2) each CDS has
                // a unique ID, and all cds records with the same parent id belong to the same isoform.

                if (gffCdsCltn.isUniqueIds()) {
                    for (BasicFeature cdsPart : gffCdsCltn.getParts()) {
                        parent.addUTRorCDS(cdsPart);
                    }

                } else {
                    Map<String, List<BasicFeature>> cdsPartsMap = gffCdsCltn.getPartsById();
                    boolean first = true;
                    for (Map.Entry<String, List<BasicFeature>> entry : cdsPartsMap.entrySet()) {
                        String cdsId = entry.getKey();
                        List<BasicFeature> cdsParts = entry.getValue();

                        BasicFeature isoform;
                        if(first) {
                            isoform = parent;
                        }
                        else {
                            isoform = copyForCDS(parent);
                            igvFeatures.add(isoform);
                        }
                        for (BasicFeature cds : cdsParts) {
                            isoform.addUTRorCDS(cds);
                        }
                        first = false;
                    }
                }
            }

            igvFeatures.addAll(gffFeatures.values());

            for (Feature bf : igvFeatures) {
                ((BasicFeature) bf).sortExons();
            }

            FeatureUtils.sortFeatureList(igvFeatures);
            return igvFeatures;

        }


        private static BasicFeature copyForCDS(BasicFeature bf) {

            BasicFeature copy = new BasicFeature(bf.getChr(), bf.getStart(), bf.getEnd(), bf.getStrand());
            copy.setName(bf.getName());
            copy.setColor(bf.getColor());
            copy.setIdentifier(bf.getIdentifier());
            copy.setURL(bf.getURL());
            copy.setType(bf.getType());
            for (Exon ex : bf.getExons()) {
                Exon newExon = new Exon(ex);
                newExon.setNonCoding(true);
                copy.addExon(newExon);

            }

            return copy;
        }

        /**
         * Container to hold all the cds records for a given parent (usually an mRNA).
         */
        public static class GFFCdsCltn {
            String parentId;
            List<BasicFeature> cdsParts;
            String chr;
            int start = Integer.MAX_VALUE;
            int end = Integer.MIN_VALUE;
            Strand strand;

            public GFFCdsCltn(String parentId) {
                this.parentId = parentId;
                cdsParts = new ArrayList(5);
            }

            public void addPart(BasicFeature bf) {
                cdsParts.add(bf);
                this.chr = bf.getChr();
                this.start = Math.min(this.start, bf.getStart());
                this.end = Math.max(this.end, bf.getEnd());
                this.strand = bf.getStrand();
            }

            public String getParentId() {
                return parentId;
            }

            public List<BasicFeature> getParts() {
                return cdsParts;
            }

            public Map<String, List<BasicFeature>> getPartsById() {

                Map<String, List<BasicFeature>> map = new HashMap<String, List<BasicFeature>>();
                for (BasicFeature bf : cdsParts) {
                    String id = bf.getIdentifier();
                    List<BasicFeature> parts = map.get(id);
                    if (parts == null) {
                        parts = new ArrayList<BasicFeature>();
                        map.put(id, parts);
                    }
                    parts.add(bf);

                }
                return map;
            }

            public boolean isUniqueIds() {
                if (cdsParts.isEmpty()) return true;
                Iterator<BasicFeature> iter = cdsParts.iterator();
                String firstId = iter.next().getIdentifier();
                while (iter.hasNext()) {
                    BasicFeature bf = iter.next();
                    if (firstId.equals(bf.getIdentifier())) return false;
                }
                return true;
            }
        }

    }
}
