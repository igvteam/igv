package org.broad.igv.sam.reader;

import org.broad.igv.sam.Alignment;
import org.broad.igv.tools.IGVToolsTest;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/03/01
 */
public class AlignmentReaderFactoryTest {
    @Before
    public void setUp() throws Exception {
        TestUtils.setUpHeadless();

    }

    @After
    public void tearDown() throws Exception {

    }


    @Test
    public void testGetMergedReader() throws Exception {
        String relfiname = "HG00171.hg18.bam";
        String path = TestUtils.DATA_DIR + "/out/temprel.bam.list";
        String[] relnames = IGVToolsTest.generateRepLargebamsList(path, relfiname, 3, false);
        for(String s: relnames){
            assertFalse((new File(s)).isAbsolute());
        }
        
        String apath = TestUtils.DATA_DIR + "/out/tempabs.bam.list";
        String[] absnames = IGVToolsTest.generateRepLargebamsList(apath, relfiname, 3, true);
        for(String s: absnames){
            assertTrue((new File(s)).isAbsolute());
        }

        AlignmentReader relreader = AlignmentReaderFactory.getReader(path, false);
        AlignmentReader absreader = AlignmentReaderFactory.getReader(apath, false);
        Iterator<Alignment> iter = relreader.iterator();
        Iterator<Alignment> absiter = absreader.iterator();

        while(iter.hasNext()){
            Alignment relA = iter.next();
            Alignment absA = absiter.next();
            //Only check 1 field for speed, these are long strings so unlikely to be equal by accident
            assertEquals(absA.getCigarString(), relA.getCigarString());
        }
    }

    public static void main(String[] args){
        String path = "http://www.example.com";
        File fi = new File(path);
        System.out.println(fi.isAbsolute());
    }


}
