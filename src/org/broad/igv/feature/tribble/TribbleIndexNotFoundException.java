/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature.tribble;

/**
 * Thrown when an index is not found, but required.  This class extends Exception, as opposed to RuntimeException, so
 * clients will be forced to deal with it.
 *
 * @author jrobinso
 *         Date: 11/16/13
 *         Time: 10:10 PM
 */
public class TribbleIndexNotFoundException extends Exception{

    public TribbleIndexNotFoundException(String message) {
        super(message);
    }
}
