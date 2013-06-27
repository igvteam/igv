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

package org.broad.igv.plugin.mongovariant;

import org.apache.log4j.Logger;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.TruthStatus;

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
                long date = (Long) super.getAttributes().get(key);
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
