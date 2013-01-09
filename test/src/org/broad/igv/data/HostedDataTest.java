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

package org.broad.igv.data;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.reader.BAMHttpReader;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.ui.ResourceTree;
import org.broad.igv.ui.action.LoadFromServerAction;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;
import org.w3c.dom.Document;

import javax.swing.tree.DefaultMutableTreeNode;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import static org.junit.Assert.*;

/**
 * Test loading the data we host. Sometimes it gets corrupted,
 * or we make backwards-incompatible changes
 * <p/>
 * User: jacob
 * Date: 2012-Jul-27
 */
public class HostedDataTest extends AbstractHeadlessTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 1200e4);

    private PrintStream errorWriter = System.out;

    /**
     * Test loading all the data hosted for each genome
     *
     * @throws Exception
     */

    @Test
    public void testLoadServerData() throws Exception {
        DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd-HH-mm-ss");
        Date date = new Date();

        String outPath = TestUtils.DATA_DIR + "failed_loaded_files_" + dateFormat.format(date) + ".txt";
        errorWriter = new PrintStream(outPath);

        List<GenomeListItem> serverSideGenomeList = getServerGenomes();


        Map<ResourceLocator, Exception> failedFiles = new LinkedHashMap<ResourceLocator, Exception>(10);
        Set<ResourceLocator> loadedResources = new HashSet<ResourceLocator>(1000);
        LinkedHashSet<String> nodeURLs;

        int counter = 0;
        //Large BAM files have have index files ~10 Mb, don't want to use
        //too much space on disk at once
        int clearInterval = 50;

        DefaultMutableTreeNode maxNode = null;
        int maxChildCount = -1;

        for (GenomeListItem genomeItem : serverSideGenomeList) {
            //Do this within the loop, both to make sure we get a fresh genome
            //and not use too much disk space
            GenomeManager.getInstance().clearGenomeCache();

            String genomeURL = LoadFromServerAction.getGenomeDataURL(genomeItem.getId());

            TrackLoader loader = new TrackLoader();
            Genome curGenome = GenomeManager.getInstance().loadGenome(genomeItem.getLocation(), null);

            errorWriter.println("Genome: " + curGenome.getId());

            try {
                nodeURLs = LoadFromServerAction.getNodeURLs(genomeURL);
                if (nodeURLs == null) {
                    errorWriter.println("Warning: No Data found for " + genomeURL);
                    continue;
                }
            } catch (Exception e) {
                recordError(genomeURL, e, failedFiles);
                continue;
            }

            for (String nodeURL : nodeURLs) {

                errorWriter.println("NodeURL: " + nodeURL);

                try {
                    Document xmlDocument = LoadFromServerAction.createMasterDocument(Arrays.asList(nodeURL));
                    DefaultMutableTreeNode treeNode = new DefaultMutableTreeNode("HostedDataTest");
                    ResourceTree.buildLocatorTree(treeNode, xmlDocument.getDocumentElement(),
                            Collections.<ResourceLocator>emptySet(), null);

                    Enumeration enumeration = treeNode.depthFirstEnumeration();

                    while (enumeration.hasMoreElements()) {
                        Object nextEl = enumeration.nextElement();
                        DefaultMutableTreeNode node = (DefaultMutableTreeNode) nextEl;

                        Object userObject = node.getUserObject();
                        //Get resource locator from tree
                        //don't load resources we've already tried (same file can be listed multiple times)
                        ResourceTree.CheckableResource checkableResource;
                        ResourceLocator locator;
                        if (userObject instanceof ResourceTree.CheckableResource) {
                            checkableResource = (ResourceTree.CheckableResource) userObject;
                            locator = checkableResource.getResourceLocator();
                            if (locator.getPath() == null || loadedResources.contains(locator)) {
                                continue;
                            } else {
                                loadedResources.add(locator);
                            }
                        } else {
                            continue;
                        }

//                        int childCount = node.getChildCount();
//                        if(childCount > 0){
//                            System.out.println(node.getUserObject()  + " Children: " + childCount);
//                            if(childCount > maxChildCount){
//                                maxChildCount = childCount;
//                                maxNode = node;
//                            }
//                        }

                        FeatureDB.clearFeatures();

                        try {
//                            if(locator.getServerURL() != null){
//                                //System.out.println("server url " + locator.getServerURL());
//                                //System.out.println("path " + locator.getPath());
//                            }else{
//                                continue;
//                            }
//                            errorWriter.println("Loading " + locator);
                            loader.load(locator, curGenome);
                        } catch (Exception e) {
                            recordError(locator, e, failedFiles);
                        }

                        counter = (counter + 1) % clearInterval;
                        if (counter == 0) {
                            BAMHttpReader.cleanTempDir(1l);
                        }
                    }

                } catch (Exception e) {
                    recordError(nodeURL, e, failedFiles);
                }
            }

        }

        //System.out.println("Max Node: " + maxNode + ". Children: " + maxChildCount);

        for (Map.Entry<ResourceLocator, Exception> entry : failedFiles.entrySet()) {
            ResourceLocator item = entry.getKey();
            errorWriter.println(formatLocator(item) + "\terror: " + entry.getValue().getMessage());
        }

        errorWriter.flush();
        errorWriter.close();
        assertEquals(0, failedFiles.size());

    }

    private String formatLocator(ResourceLocator locator) {
        return String.format("Name: %s\tPath: %s\t serverURL: %s",
                locator.getName(), locator.getPath(), locator.getServerURL());
    }

    private void recordError(String path, Exception e, Map<ResourceLocator, Exception> failures) {
        ResourceLocator locator = new ResourceLocator(path);
        recordError(locator, e, failures);
    }

    private void recordError(ResourceLocator locator, Exception e, Map<ResourceLocator, Exception> failures) {
        failures.put(locator, e);
        errorWriter.println(formatLocator(locator) + "\terror: " + e.getMessage());
        errorWriter.println("StackTrace: ");
        for (StackTraceElement el : e.getStackTrace()) {
            errorWriter.println(el);
        }
        errorWriter.flush();
    }

    private List<GenomeListItem> getServerGenomes() throws IOException {
        String genomeListPath = PreferenceManager.DEFAULT_GENOME_URL;
        PreferenceManager.getInstance().overrideGenomeServerURL(genomeListPath);
        List<GenomeListItem> serverSideItemList = GenomeManager.getInstance().getServerGenomeArchiveList(null);
        assertNotNull("Could not retrieve genome list from server", serverSideItemList);
        assertTrue("Genome list empty", serverSideItemList.size() > 0);
        return serverSideItemList;
    }


    @Ignore
    @Test
    public void testLoadServerGenomes() throws Exception {

        List<GenomeListItem> serverSideItemList = getServerGenomes();

        Map<GenomeListItem, Exception> failedGenomes = new LinkedHashMap<GenomeListItem, Exception>(10);

        int count = 0;
        for (GenomeListItem genome : serverSideItemList) {
            try {
                count++;
                tstLoadGenome(genome.getLocation());
                Runtime.getRuntime().gc();
            } catch (Exception e) {
                failedGenomes.put(genome, e);
            }
        }
        errorWriter.println("Attempted to load " + count + " genomes");
        errorWriter.println(failedGenomes.size() + " of them failed");
        for (Map.Entry<GenomeListItem, Exception> entry : failedGenomes.entrySet()) {
            GenomeListItem item = entry.getKey();
            System.out.println(String.format("Exception loading (%s\t%s\t%s): %s", item.getDisplayableName(),
                    item.getLocation(), item.getId(), entry.getValue()));
        }

        assertEquals(0, failedGenomes.size());
    }

    public void tstLoadGenome(String path) throws Exception {
        FeatureDB.clearFeatures();
        Genome genome = GenomeManager.getInstance().loadGenome(path, null);
        assertTrue(genome.getAllChromosomeNames().size() > 0);
    }

}
