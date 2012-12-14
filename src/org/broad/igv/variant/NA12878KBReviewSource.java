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

package org.broad.igv.variant;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.*;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Dec-14
 */
public class NA12878KBReviewSource implements FeatureSource<VCFVariant> {

    private int featureWindowSize = 1000000;
    private NA12878DBArgumentCollection args;
    private NA12878KnowledgeBase kb;
    private GenomeLocParser parser;

    public NA12878KBReviewSource(String dbSpecPath){
        this.args = new NA12878DBArgumentCollection(dbSpecPath);
        parser = createGenomeLocParser();

    }

    private void initKB() {
        kb = new NA12878KnowledgeBase(parser, this.args);
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
        if(kb == null){
            initKB();
        }

        SiteSelector criteria = new SiteSelector(parser);
        //Convert from 0-based to 1-based
        criteria.addInterval(chr, start + 1, end).onlyReviewed();
        SiteIterator<MongoVariantContext> iterator  = kb.getCalls(criteria);
        List<VCFVariant> variants = new ArrayList<VCFVariant>();
        while(iterator.hasNext()){
            MongoVariantContext mvc = iterator.next();
            VCFVariant vcf = new VCFVariant(mvc.getVariantContext(), mvc.getChr());
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

    private GenomeLocParser createGenomeLocParser(){
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        for(Chromosome chr: GenomeManager.getInstance().getCurrentGenome().getChromosomes()){
            dict.addSequence(new SAMSequenceRecord(chr.getName(), chr.getLength()));
        }
        return new GenomeLocParser(dict);
    }

    static MongoVariantContext createMVC(int allele0, int allele1, String callsetName, VariantContext variantContext, TruthStatus truthStatus){
        List<Allele> alleleList = variantContext.getAlleles();

        MongoGenotype mgt = new MongoGenotype(allele0, allele1);
        Genotype gt = mgt.toGenotype(alleleList);
        MongoVariantContext mvc = MongoVariantContext.create(callsetName, variantContext, truthStatus, gt);
        mvc.setReviewed(true);
        return mvc;
    }

}
