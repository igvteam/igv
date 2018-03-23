/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.plugin.mongovariant;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broad.igv.annotations.ForTesting;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.Track;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.*;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Dec-14
 */
public class VariantReviewSource implements FeatureSource<VCFVariant> {

    private static Logger log = Logger.getLogger(VariantReviewSource.class);

    private int featureWindowSize = 1000000;
    private NA12878DBArgumentCollection args;
    private NA12878KnowledgeBase kb;
    private GenomeLocParser parser;

    /**
     * Whether to return consensus sites only, as opposed to all sites
     */
    @ForTesting
    boolean consensusOnly = true;

    public VariantReviewSource(ResourceLocator locator) {
        this.args = new NA12878DBArgumentCollection(locator.getPath());
        parser = createGenomeLocParser();
    }

    private void initKB() {
        kb = new NA12878KnowledgeBase(parser, this.args);
        log.info("Connected to " + this.args.getLocator());
    }

    private void closeKB() {
        if (kb != null) {
            kb.close();
            kb = null;
        }
    }

    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        this.closeKB();
    }

    @Override
    public Iterator<VCFVariant> getFeatures(String chr, int start, int end) throws IOException {
        if (kb == null) {
            initKB();
        }

        //Convert from 0-based to 1-based
        GenomeLoc loc = this.parser.createGenomeLoc(chromoNameToStandard(chr), start + 1, end);
        GenomeLocSortedSet glss = new GenomeLocSortedSet(this.parser, loc);
        SiteManager manager = new SiteManager(parser, glss, parser.getContigs());

        SiteIterator<MongoVariantContext> iterator;
        if (consensusOnly) {
            iterator = kb.getConsensusSites(manager);
        } else {
            iterator = kb.getCalls(manager);
        }
        List<VCFVariant> variants = new ArrayList<VCFVariant>();
        while (iterator.hasNext()) {
            MongoVariantContext mvc = iterator.next();
            MongoVCFVariant vcf = new MongoVCFVariant(mvc, mvc.getChr());
            variants.add(vcf);
        }
        return variants.iterator();
    }


    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    @Override
    public int getFeatureWindowSize() {
        return featureWindowSize;
    }

    @Override
    public void setFeatureWindowSize(int size) {
        this.featureWindowSize = size;
    }

    /**
     * TODO This chromosome replacement is a hack
     * The db uses digits (or X/Y, note the case), IGV uses "chr#"
     *
     * @param chromoName
     * @return
     */
    static String chromoNameToStandard(String chromoName) {
        if(chromoName.toLowerCase().startsWith("chr")){
            return chromoName.substring(3);
        }else{
            return chromoName;
        }
    }

    private GenomeLocParser createGenomeLocParser() {
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        for (Chromosome chr : GenomeManager.getInstance().getCurrentGenome().getChromosomes()) {
            dict.addSequence(new SAMSequenceRecord(chromoNameToStandard(chr.getName()), chr.getLength()));
        }
        return new GenomeLocParser(dict);
    }

    static MongoVariantContext createMVC(int allele0, int allele1, String callsetName, VariantContext variantContext, TruthStatus truthStatus, boolean isComplexEvent) {
        List<Allele> alleleList = variantContext.getAlleles();

        MongoGenotype mgt = new MongoGenotype(allele0, allele1);
        Genotype gt = mgt.toGenotype(alleleList);
        MongoVariantContext mvc = MongoVariantContext.create(callsetName, variantContext, truthStatus, gt);
        mvc.setReviewed(true);
        mvc.setChr(chromoNameToStandard(mvc.getChr()));
        mvc.setIsComplexEvent(isComplexEvent);
        mvc.setConfidence(NA12878KnowledgeBase.InputCallsetConfidence.REVIEW.confidence);
        return mvc;
    }

    public static VariantTrack loadVariantReview(ResourceLocator locator, List<Track> newTracks) {
        //TODO Figure out how to name the samples properly
        List<String> allSamples = Collections.emptyList();
        VariantReviewSource source = new VariantReviewSource(locator);
        VariantTrack track = new VariantTrack(locator, source, allSamples, false);
        track.setRenderer(new VariantReviewRenderer(track));
        newTracks.add(track);
        track.setMargin(0);
        return track;
    }

}
