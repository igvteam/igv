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

import org.apache.log4j.Logger;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.TruthStatus;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Wrapper for VCFVariant so we can include additional fields from MongoDB
 * User: jacob
 * Date: 2013-Jan-28
 */
public class MongoVCFVariant extends VCFVariant {

    private static Logger log = Logger.getLogger(MongoVCFVariant.class);

    private final MongoVariantContext mongoVariantContext;

    @Override
    public String getAttributeAsString(String key) {
        if (key.equalsIgnoreCase("date")) {
            try {
                long date = Long.parseLong((String) super.getAttributes().get(key));
                SimpleDateFormat formatter = new SimpleDateFormat("MMM dd, yyyy. HH:mm");
                return formatter.format(new Date(date));
            } catch (Exception e) {
                //Don't much care, but if we drop the ball and the API changes
                //we don't want to crash everything
                log.error(e.getMessage(), e);
            }
        }
        return super.getAttributeAsString(key);
    }

    public MongoVCFVariant(MongoVariantContext mongoVariantContext, String chr) {
        super(mongoVariantContext.getVariantContext(), chr);
        this.mongoVariantContext = mongoVariantContext;
    }

    public TruthStatus getTruthStatus() {
        return mongoVariantContext.getType();
    }

    public boolean isReviewed() {
        return mongoVariantContext.isReviewed();
    }
}
