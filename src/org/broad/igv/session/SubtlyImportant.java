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

package org.broad.igv.session;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Annotates something which seem like it's not used anywhere,
 * but is. An example might be private getter/setter/constructors
 * used in tracks which JAXB requires. Since these are only
 * called reflectively, and IDE might report them as
 * never being used.
 *
 * User: jacob
 * Date: 2013-Jan-08
 */
@Retention(RetentionPolicy.SOURCE)
public @interface SubtlyImportant {
    String description() default "";

    /**
     * Comma-separated list of methods which
     * rely on this
     * @return
     */
    String whereUsed() default "";
}
