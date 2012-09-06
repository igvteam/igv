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

package org.broad.igv.feature.genome;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.FilenameFilter;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Sep-06
 */
public class GenomeImporterTest extends AbstractHeadlessTest {

    @Test
    public void testCreateGenomeArchiveFromDir() throws Exception {

        File genomeFile = new File(TestUtils.DATA_DIR, "out/testSetGenome.genome");
        String genomeId = "testSet";
        String genomeDisplayName = genomeId;
        String fastaPath = TestUtils.DATA_DIR + "fasta/set";

        deleteFaiFiles(fastaPath);

        File genomeArchive = (new GenomeImporter()).createGenomeArchive(
                genomeFile, genomeId, genomeDisplayName, fastaPath,
                null, null, null);

        assertNotNull(genomeArchive);
        assertTrue(genomeArchive.exists());

        deleteFaiFiles(fastaPath);
    }

    /**
     * Deletes all fasta index files in the provided path
     *
     * @param dir
     */
    private void deleteFaiFiles(String dir) {
        //Delete index files, if they exist
        File fastaDir = new File(dir);
        File[] idxFiles = fastaDir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".fai");
            }
        });

        for (File idxFile : idxFiles) {
            idxFile.delete();
        }
    }
}
