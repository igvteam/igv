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

package org.broad.igv.dev.api;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Annotation for a method which indicates it is part of IGVs
 * public API. This method should be documented (correctly) and
 * only removed with appropriate notice (see deprecation policy)
 *
 * TODO Write deprecation policy
 * TODO Note: API doesn't actually exist yet
 * User: jacob
 * Date: 2012-Dec-17
 */
@Retention(RetentionPolicy.SOURCE)
public @interface api{

    /**
     * The version in which this method was created
     * @return
     */
    String since() default "alpha";

    /**
     * The version in which this method will be removed.
     * For example, a value of "2.2" would mean this
     * method will be available in IGV 2.2, but not 2.3
     * @return
     */
    String until() default "";
}
