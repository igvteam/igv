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
        GenomeManager.getInstance().setCurrentGenomeForTest(null);
        genome = TestUtils.loadGenome();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        TestUtils.clearOutputDir();
        GenomeManager.getInstance().setCurrentGenomeForTest(null);
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
