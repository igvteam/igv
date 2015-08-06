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

package org.broad.igv.ga4gh;

import com.google.gson.JsonObject;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

public class Ga4ghAPIHelperTest {

    static Ga4ghProvider provider = Ga4ghAPIHelper.GA4GH_GOOGLE_PROVIDER;


    @Test
    public void testReadGroupsetsSearch() throws Exception {

        String datasetId = "10473108253681171589";

        List<Ga4ghReadset> readGroupsets = Ga4ghAPIHelper.searchReadGroupsets(provider, datasetId, 1000);

        assertNotNull(readGroupsets);

    }

    @Test
    public void testReferenceSearch() throws Exception {

        String referenceSetId = "EOSt9JOVhp3jkwE";

        List<JsonObject> references = Ga4ghAPIHelper.searchReferences(provider, referenceSetId, 1000);

        assertNotNull(references);

    }


    @Test
    public void testReads() throws Exception {

        String readGroupsetId = "CMvnhpKTFhCwvIWYw9eikzQ";
    }
}