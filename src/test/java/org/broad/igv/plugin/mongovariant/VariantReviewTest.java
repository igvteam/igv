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

import com.mongodb.WriteResult;
import junit.framework.Assert;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.variant.Allele;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.NA12878KnowledgeBase;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.TruthStatus;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.junit.Assume;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Test retrieving variants from mongo database.
 * This test requires that the relevant mongo db be running,
 * so we ignore it because we can't be sure about foreign conditions.
 * <p/>
 * User: jacob
 * Date: 2012-Dec-14
 */
@Ignore("Test only intended to be run manually")
public class VariantReviewTest extends AbstractHeadlessTest {

    String dbSpecPath = "resources/NA12878kb_local.json";

    private int allele0 = 0;
    private int allele1 = 1;
    private String callsetName = "test_callset";
    private TruthStatus truthStatus = TruthStatus.SUSPECT;

    private String chr = VariantReviewSource.chromoNameToStandard("chr10");
    //0-based coords
    private int start = (int) 1e6;
    private int end = start + 1;
    private List<String> alleles = Arrays.asList("A", "G");
    private MongoVariantContext mvc;

    private VariantReviewSource source;

    @Before
    public void setUp() throws Exception {
        super.setUp();
        VariantContextBuilder builder = new VariantContextBuilder();
        //Convert from exclusive end to inclusive end
        builder.start(start + 1).stop(end).chr(chr).alleles(alleles);
        VariantContext vc = builder.make();
        mvc = VariantReviewSource.createMVC(allele0, allele1, callsetName, vc, truthStatus, false);

        int errorsResetting = 0;
        try {
            errorsResetting = resetDB();
        } catch (Exception e) {
            System.out.println(e);
            Assume.assumeNoException(e);
        }
        Assume.assumeTrue(errorsResetting == 0);

        source = new VariantReviewSource(new ResourceLocator(dbSpecPath));
    }

    @Override
    public void tearDown() throws Exception {
        super.tearDown();
        try {
            resetDB();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private boolean checkFeatureNotPresent() throws Exception {
        Iterator<VCFVariant> result = getFeature();
        return !result.hasNext();
    }

    /*
     * Basically here to check the remove call
     */
    @Test
    public void testFeatureNotPresent() throws Exception {
        boolean featNotPresent = checkFeatureNotPresent();
        if (!featNotPresent) {
            Iterator<VCFVariant> result = getFeature();
            while (result.hasNext()) {
                System.out.println(result.next());
            }
        }
        assertTrue(featNotPresent);
    }

    //TODO Separate into add/get methods, but that requires prepopulation of data
    @Test
    public void testAddGetFeature() throws Exception {
        source.consensusOnly = false;
        Assume.assumeTrue(checkFeatureNotPresent());

        String errorMessage = VariantReviewDialog.addCall(dbSpecPath, mvc);
        if (errorMessage != null) System.out.println(errorMessage);
        assertNull(errorMessage);

        Iterator<VCFVariant> result = getFeature();
        Assert.assertTrue("Empty result after adding call", result.hasNext());
        VCFVariant variant = result.next();
        Assert.assertFalse("Empty result after adding call", result.hasNext());

        assertEquals(chr, variant.getChr());
        assertEquals(start, variant.getStart());
        assertEquals(end, variant.getEnd());
        assertEquals(alleles.get(0), variant.getReference().replace("*", ""));
        int index = 1;
        for (Allele al : variant.getAlternateAlleles()) {
            assertEquals(alleles.get(index++), al.toString());
        }
        assertTrue(variant.getSource().equals(callsetName));

    }

    private Iterator<VCFVariant> getFeature() throws IOException {
        return source.getFeatures(chr, start, end);
    }

    private int resetDB() throws Exception {
        NA12878DBArgumentCollection args = new NA12878DBArgumentCollection(dbSpecPath);
        NA12878KnowledgeBase kb = new NA12878KnowledgeBase(null, args);

        List<WriteResult> writeResults = kb.removeCall(mvc);
        kb.close();
        int errCount = 0;
        for (WriteResult wr : writeResults) {
            if (wr.getError() != null) {
                System.out.println(wr.getError());
                errCount++;
            }
        }
        return errCount;
    }
}
