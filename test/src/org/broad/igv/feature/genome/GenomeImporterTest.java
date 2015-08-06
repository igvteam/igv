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

package org.broad.igv.feature.genome;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.FilenameFilter;
import java.util.zip.ZipFile;

import static junit.framework.Assert.*;

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

    @Test
    public void testCreateGenomeArchiveFromFiles() throws Exception {

        File genomeFile = new File(TestUtils.DATA_DIR, "out/testSetGenome.genome");
        String genomeId = "testSet";
        String genomeDisplayName = genomeId;
        String fastaPath = TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta";
        File fastaFile = new File(fastaPath);

        String parentDir = TestUtils.DATA_DIR + "genomes/genome_raw_files/hg18.unittest";
        File cytobandFile = new File(parentDir, "hg18_cytoBand.txt");
        File geneAnnotFile = new File(parentDir, "hg18_refGene_head1k.txt");
        File chrAliasFile = new File(parentDir, "pointless_alias.tab");

        deleteFaiFiles(fastaPath);

        File genomeArchive = (new GenomeImporter()).createGenomeArchive(
                genomeFile, genomeId, genomeDisplayName, fastaPath,
                geneAnnotFile, cytobandFile, chrAliasFile);

        assertNotNull(genomeArchive);
        assertTrue(genomeArchive.exists());

        GenomeDescriptor descriptor = GenomeManager.getInstance().parseGenomeArchiveFile(genomeArchive);

        assertEquals(cytobandFile.getName(), descriptor.cytoBandFileName);
        assertEquals(geneAnnotFile.getName(), descriptor.getGeneFileName());
        assertEquals(chrAliasFile.getName(), descriptor.chrAliasFileName);
        assertEquals(fastaFile.getAbsolutePath(), descriptor.getSequenceLocation());

        assertTrue(descriptor.hasCytobands());
        assertTrue(descriptor.isChromosomesAreOrdered());
        assertTrue(descriptor.isFasta());

        //Check that files seem to be accurate. We just look at sizes
        ZipFile zipFile = new ZipFile(genomeArchive);
        File[] files = new File[]{cytobandFile, geneAnnotFile, chrAliasFile};
        for(File file: files){
            assertEquals("File sizes unequal for " + file.getName(), file.length(), zipFile.getEntry(file.getName()).getSize());
        }

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

        if(idxFiles == null) return;

        for (File idxFile : idxFiles) {
            idxFile.delete();
        }
    }
}
