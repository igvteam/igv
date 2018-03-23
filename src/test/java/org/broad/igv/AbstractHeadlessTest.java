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


package org.broad.igv;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.util.TestUtils;
import org.junit.*;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;
import org.junit.runner.RunWith;
import util.IGVTestRunner;

import java.io.IOException;
import java.io.PrintStream;


/**
 * General setup/takedown/rules for headless tests
 * User: jacob
 * Date: 2012/05/17
 */
@Ignore
@RunWith(IGVTestRunner.class)
public class AbstractHeadlessTest {

    protected static Genome genome;

    protected PrintStream oldOut = System.out;

    @Rule
    public TestRule testTimeout = new Timeout((int) 30e9);

    @BeforeClass
    public static void setUpClass() throws Exception {
        setUpHeadless();
        GenomeManager.getInstance().setCurrentGenome(null);
        genome = TestUtils.loadGenome();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        TestUtils.clearOutputDir();
        GenomeManager.getInstance().setCurrentGenome(null);
    }

    @Before
    public void setUp() throws Exception {
        oldOut = System.out;
        TestUtils.resetPrefsFile();
        TestUtils.resetTestUserDefinedGenomes();
    }

    @After
    public void tearDown() throws Exception {
        TestUtils.resetPrefsFile();
        TestUtils.resetTestUserDefinedGenomes();
        TestUtils.clearOutputDir();
        System.setOut(oldOut);
    }


    private static void setUpHeadless() throws IOException {
        Globals.setHeadless(true);
        TestUtils.setUpTestEnvironment();
    }
}
