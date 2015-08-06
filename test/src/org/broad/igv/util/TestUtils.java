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

package org.broad.igv.util;

import com.google.common.base.Function;
import com.google.common.base.Supplier;
import htsjdk.samtools.util.ftp.FTPClient;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.AsciiLineReader;
import org.junit.Assert;
import org.junit.Ignore;
import org.w3c.dom.Document;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.namespace.QName;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.awt.*;
import java.io.*;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.lang.management.ThreadInfo;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.*;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

/**
 * @author jrobinso
 * @date Jul 28, 2010
 */
@Ignore
public class TestUtils {
    public static final String DATA_DIR = "test/data/";
    public static final String TMP_OUTPUT_DIR = DATA_DIR + "out/";
    public static final String defaultGenome = DATA_DIR + "/genomes/hg18.unittest.genome";
    public static final String AVAILABLE_FTP_URL = "ftp://ftp.broadinstitute.org/pub/igv/TEST/test.txt";
    public static final String UNAVAILABLE_FTP_URL = "ftp://www.example.com/file.txt";
    //This is so ant can set the large data directory
    private static String LARGE_DATA_DIR_KEY = "LARGE_DATA_DIR";
    public static String LARGE_DATA_DIR = "test/largedata/";

    static {
        LARGE_DATA_DIR = System.getProperty(LARGE_DATA_DIR_KEY, LARGE_DATA_DIR);
    }

    public static void setUpTestEnvironment() throws IOException {
        Globals.setTesting(true);
        //Globals.setBatch(true);
        resetPrefsFile();
        resetTestUserDefinedGenomes();
        Globals.READ_TIMEOUT = 60 * 1000;
        Globals.CONNECT_TIMEOUT = 60 * 1000;
        FTPClient.READ_TIMEOUT = 60 * 1000;

        //Create output directory if it doesn't exist
        File outDir = new File(TestUtils.TMP_OUTPUT_DIR);
        if (!outDir.exists()) {
            outDir.mkdir();
        }
        clearOutputDir();
    }

    public static void resetPrefsFile(){
        File prefsFile = new File("testprefs.properties");
        prefsFile.delete();
        prefsFile.deleteOnExit();
        PreferenceManager.getInstance().setPrefsFile(prefsFile.getAbsolutePath());
    }

    /**
     * Loads the session into IGV. This blocks until the session
     * is loaded.
     *
     * @param igv
     * @param sessionPath
     * @throws InterruptedException
     */
    public static void loadSession(IGV igv, String sessionPath) throws InterruptedException {
        igv.doRestoreSession(sessionPath, null, false);
    }

    /**
     * See {@link #createIndex(String, int, int)}
     *
     * @param file
     * @throws IOException
     */
    public static void createIndex(String file) throws IOException {
        createIndex(file, IgvTools.LINEAR_INDEX, IgvTools.LINEAR_BIN_SIZE);
    }

    /**
     * Destroys index file if it exists, and creates new one under
     * the specified parameters
     *
     * @param file
     * @param indexType
     * @param binSize
     * @throws IOException
     */
    public static void createIndex(String file, int indexType, int binSize) throws IOException {
        String indexPath = (new IgvTools()).doIndex(file, null, indexType, binSize);
        File indexFile = new File(indexPath);
        indexFile.deleteOnExit();
    }

    /**
     * Load a test genome
     *
     * @return
     * @throws IOException
     */
    public static Genome loadGenome() throws IOException {
        final String genomeFile = defaultGenome;
        return IgvTools.loadGenome(genomeFile);
    }

    public static void clearOutputDir() throws IOException {
        File outputDir = new File(TMP_OUTPUT_DIR);
        if (outputDir.isDirectory()) {
            File[] listFiles = outputDir.listFiles();
            for (File fi : listFiles) {
                //Keep hidden files and directories
                if (!fi.isHidden()) {
                    if(fi.isFile()){
                        fi.delete();
                    }else if(fi.isDirectory()){
                        FileUtils.deleteDir(fi);
                    }
                }
            }
        }
    }

    /**
     * Returns either 1 or 2, representing the number of
     * bytes used to end a line. Reads only from first line of a file
     *
     * @param filePath
     * @return
     * @throws IOException
     */
    public static int getBytesAtEnd(String filePath) throws IOException {
        InputStream is = new FileInputStream(filePath);
        AsciiLineReader reader = new AsciiLineReader(is);
        String line = reader.readLine();
        int bytesThisLine = (int) reader.getPosition();

        reader.close();

        return bytesThisLine - line.length();

    }

    /**
     * Check that the features are all the same. Checks size of list, chr, start, and end of
     * each feature
     *
     * @param expected
     * @param actIter
     * @throws AssertionError
     */
    public static void assertFeatureListsEqual(Iterator<? extends Feature> expected, Iterator<? extends Feature> actIter) throws AssertionError {
        while (expected.hasNext()) {
            Feature exp = expected.next();
            assertTrue(actIter.hasNext());
            Feature act = actIter.next();
            assertFeaturesEqual(exp, act);
        }
        assertFalse(actIter.hasNext());
    }


    public static void assertFeaturesEqual(Feature exp, Feature act) {
        assertEquals(exp.getChr(), act.getChr());
        assertEquals(exp.getStart(), act.getStart());
        assertEquals(exp.getEnd(), act.getEnd());
    }

    /**
     * Name matching is case-insensitive
     * @param exp
     * @param act
     */
    public static void assertNamedFeaturesEqual(NamedFeature exp, NamedFeature act) {
        assertFeaturesEqual(exp, act);
        assertEquals(exp.getName().toUpperCase(), act.getName().toUpperCase());
    }

    public static void assertTrackLoaded(IGV igv, String trackName){
        boolean found = false;
        for(Track t: igv.getAllTracks()){
            if(t.getName().equals(trackName)){
                found = true;
                break;
            }
        }
        Assert.assertTrue("Track not loaded", found);
    }

    /**
     *
     * @param featureIterator
     * @return Number of features in the iterator
     * @throws Exception
     */
    public static int assertFeatureIteratorSorted(Iterator<? extends Feature> featureIterator){
        int lastStart = -1;
        int count = 0;
        while (featureIterator.hasNext()) {
            Feature f0 = featureIterator.next();
            assertTrue(f0.getStart() >= lastStart);
            lastStart = f0.getStart();
            count++;
        }
        return count;
    }

    /*
     * The FEST library finds components by name
     * We set the name property to equal the variable name.
     * Intended for testing ONLY
     * @param parent Element in which to set child names. Each field of this Object which is a Component
     *               will have it's name set
     * @param recursive Whether to set names on grand-child components. Note that this is all-or-none,
     *                  if this is true it goes recursively all the way down. Recursion is breadth-first,
     *                  stops at each level when an object has only non-Component fields
     */
    public static void setAllNames(Object parent, boolean recursive) {
        Field[] fields = parent.getClass().getDeclaredFields();
        java.util.List<Component> childComponents = new ArrayList<Component>(fields.length);
        try {
            for (Field f : fields) {
                Component c = null;
                f.setAccessible(true);
                try {
                    c = (Component) f.get(parent);
                } catch (ClassCastException e) {
                    continue;
                }
                //Null valued fields don't throw a CCE
                //We don't overwrite names, this should also prevent
                //infinite recursion
                if (c == null || c.getName() != null) {
                    continue;
                }
                //At this point, we've established
                //that the field in question is a Component
                c.setName(f.getName());
                childComponents.add(c);

                f.setAccessible(false);
            }

            if (recursive) {
                for (Component c : childComponents) {
                    setAllNames(c, recursive);
                }
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }


    public static void checkDeadlockedThreads() throws IllegalThreadStateException {
        ThreadMXBean tmx = ManagementFactory.getThreadMXBean();
        long[] ids = tmx.findDeadlockedThreads();
        if (ids != null) {
            ThreadInfo[] infos = tmx.getThreadInfo(ids, true, true);
            String log = "The following threads are deadlocked:\n";
            for (ThreadInfo ti : infos) {
                System.out.println(ti);
                log += ti.toString() + "\n";
            }

            throw new IllegalThreadStateException(log);
        }
    }

    public static Timer startDeadlockChecker(final int period) {

        TimerTask checker = new TimerTask() {
            @Override
            public void run() {
                //System.out.println("deadlock checker thread:" + Thread.currentThread().getName());
                checkDeadlockedThreads();
            }
        };
        Timer timer = new Timer();
        timer.scheduleAtFixedRate(checker, 0, period);
        return timer;
    }

    /**
     * Executes {@code predicate} on those objects supplied by {@code supplier},
     * timing each iteration. Goes {@code nTrials} times
     *
     * @param supplier
     * @param predicate
     * @param nTrials
     * @param <T>
     * @return The runtime, in nanoseconds, of each call of predicate with input of supplier.
     *         Array is sorted ascending
     */
    public static <T> long[] timeMethod(Supplier<T> supplier, Function<T, Void> predicate, int nTrials) {
        long total = 0;
        long[] times = new long[nTrials];
        System.gc();

        for (int tri = 0; tri < nTrials; tri++) {

            T input = supplier.get();

            long startTime = System.nanoTime();

            predicate.apply(input);

            long endTime = System.nanoTime();
            long elapsed = endTime - startTime;
            times[tri] = elapsed;
            total += elapsed;
        }

        Arrays.sort(times);
        long minTime = times[0];
        long maxTime = times[nTrials - 1];
        long median = times[times.length / 2];
        double average = (total * 1.0 / nTrials);
        double stdev = -1;
        try {
            stdev = stdev(times, (long) average);
        } catch (ArithmeticException e) {
            //pass
        }

        System.out.println(String.format("Avg. time: %2.2e sec. Median: %2.2e sec", average * 1.0 / 1e9, median * 1.0 / 1e9));
        System.out.println(String.format("Best: %2.2e sec. Worst: %2.2e sec", minTime * 1.0 / 1e9, maxTime * 1.0 / 1e9));
        System.out.println(String.format("Standard deviation: %2.2e sec", stdev / 1e9));
        System.out.println(String.format("Total time: %2.2e sec", total * 1.0 / 1e9));

        return times;
    }

    public static double average(long[] vals) {
        long sum = 0;
        double n = (double) vals.length;
        for (Long i : vals)
            sum += ((double) i);
        return sum / n;
    }

    public static double stdev(long[] vals, long average) {
        long sum = 0;
        for (Long i : vals)
            sum += Math.pow((i - average), 2);
        return Math.sqrt(sum / (vals.length - 1)); // sample
    }

    private static long benchmarkTime = -1;

    public static long getBenchmarkTime() {
        if (benchmarkTime < 0) {
            //Generate some numbers to average
            Random r = new Random();
            int numNumbers = 1000000;
            long[] vals = new long[numNumbers];
            for (int rr = 0; rr < numNumbers; rr++) {
                vals[rr] = r.nextInt();
            }
            System.gc();

            long startTime = System.nanoTime();

            double avg = average(vals);

            long endTime = System.nanoTime();
            benchmarkTime = endTime - startTime;
            System.out.println("Benchmark Time (s): " + ((double) benchmarkTime)/1e9d);
        }
        return benchmarkTime;
    }



    /**
     * Marshalls {@code inObj} and unmarshalls the result, returning the
     * unmarshalled version
     *
     * @param inObj
     * @return
     * @throws Exception
     */
    public static <T> T marshallUnmarshall(T inObj) throws Exception{

        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document doc = builder.newDocument();

        JAXBContext jc = JAXBContext.newInstance(inObj.getClass());
        Marshaller m = jc.createMarshaller();
        m.setProperty(Marshaller.JAXB_FRAGMENT, true);

        //This JAXBElement business is necessary because we don't know if we have @XmlRootElement on inObj
        JAXBElement inel = new JAXBElement(new QName("", "obj"), inObj.getClass(), inObj);
        //m.marshal(inel, System.out);
        m.marshal(inel, doc);

        Unmarshaller u = jc.createUnmarshaller();
        JAXBElement el = (JAXBElement) u.unmarshal(doc, inObj.getClass());
        return (T) el.getValue();
    }

    private static Object getField(Object object, Class clazz, String fieldName) throws Exception{
        if(clazz == null) throw new NoSuchFieldException(fieldName + " not found all the way up");
        Field field;
        try{
            field = object.getClass().getDeclaredField(fieldName);
        }catch (NoSuchFieldException e){
            return getField(object, clazz.getSuperclass(), fieldName);
        }
        field.setAccessible(true);
        return field.get(object);
    }

    /**
     * Get the specified field, ignoring access restrictions
     * @param object
     * @param fieldName
     * @return
     */
    public static Object getField(Object object, String fieldName) throws Exception{
        return getField(object, object.getClass(), fieldName);
    }

    private static Object runMethod(Object object, Class clazz, String methodName, Object... args) throws Exception{
        if(clazz == null) throw new NoSuchFieldException(methodName + " not found all the way up");
        Method method;
        try{
            method = object.getClass().getDeclaredMethod(methodName);
        }catch (NoSuchMethodException e){
            return runMethod(object, clazz.getSuperclass(), methodName, args);
        }
        method.setAccessible(true);
        return method.invoke(object, args);
    }

    public static Object runMethod(Object object, String methodName, Object... args) throws Exception{
        return runMethod(object, object.getClass(), methodName, args);
    }


    private static Map<String, String> replaceMap = new HashMap<String, String>(2);
    static{
        replaceMap.put("${DATA_DIR}", TestUtils.DATA_DIR);
        replaceMap.put("${LARGE_DATA_DIR}", TestUtils.LARGE_DATA_DIR);
    }
    /**
     * Mainly for session files, we want use some test paths stored in otherwise hardcoded
     * test data files. This re-writes the input text file
     * @param inputPath
     * @return File of the output location
     */
    public static File replaceTestPaths(File inputPath) throws Exception{

        BufferedReader reader = new BufferedReader(new FileReader(inputPath));

        File outputFile = new File("tmpsession.xml");
        outputFile.delete();
        outputFile.deleteOnExit();
        PrintWriter writer = new PrintWriter(new FileWriter(outputFile));

        String line;
        while((line = reader.readLine()) != null){
            for(Map.Entry<String, String> entry: replaceMap.entrySet()){
                line = line.replace(entry.getKey(), entry.getValue());
            }
            writer.println(line);
        }
        writer.flush();
        writer.close();
        return outputFile;
    }

    /**
     * Note: THIS METHOD IS VERY FRAGILE, NOT GUARANTEED TO WORK ON ALL JVMs
     * @return
     */
    static long getProcessID(){
        RuntimeMXBean bean = ManagementFactory.getRuntimeMXBean();

        // The name representing the running Java virtual machine.
        // It returns something like 12345@blah. The value
        // before the @ symbol is the PID.
        String jvmName = bean.getName();
        return Long.valueOf(jvmName.split("@")[0]);
    }

    /**
     * Uses lsof, only works on *nix systems
     *
     * @return
     */
    public static int getNumberOpenFileHandles() throws IOException{
        long pid = getProcessID();
        String result = RuntimeUtils.executeShellCommand(new String[]{"lsof", "-p", "" + pid}, null, null, false);
        return result.split("\n").length;
    }

    public static void resetTestUserDefinedGenomes() throws IOException{
        File userDefinedGenomeListFile = new File(DirectoryManager.getGenomeCacheDirectory(), GenomeManager.TEST_USER_DEFINED_GENOME_LIST_FILE);
        userDefinedGenomeListFile.delete();
        userDefinedGenomeListFile.deleteOnExit();

        Collection<GenomeListItem> userDefined = GenomeManager.getInstance().getUserDefinedGenomeArchiveList();
        userDefined.clear();
        GenomeManager.getInstance().buildGenomeItemList();
    }
}
