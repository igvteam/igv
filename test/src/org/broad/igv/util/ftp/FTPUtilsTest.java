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

package org.broad.igv.util.ftp;

import net.sf.samtools.util.ftp.FTPUtils;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.net.URL;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

/**
 * @author jrobinso
 * @date Feb 28, 2011
 */
public class FTPUtilsTest {

    @Test
    public void testResourceAvailable() throws Exception {
        assertTrue(FTPUtils.resourceAvailable(new URL(TestUtils.AVAILABLE_FTP_URL)));
        assertFalse(FTPUtils.resourceAvailable(new URL(TestUtils.UNAVAILABLE_FTP_URL)));
    }


}
