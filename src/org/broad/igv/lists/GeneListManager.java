/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.lists;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 * @date Sep 26, 2010
 */
public class GeneListManager {

    private static Logger log = Logger.getLogger(GeneListManager.class);

    public static final String[] DEFAULT_GENE_LISTS = {
            "examples.gmt", /*"biocarta_cancer_cp.gmt",*/  "reactome_cp.gmt", "kegg_cancer_cp.gmt"};


    public static LinkedHashSet<String> groups = new LinkedHashSet();

    private static LinkedHashMap<String, GeneList> geneLists = new LinkedHashMap();


    /**
     * Static initializer, called once when class is loaded.
     */
    static {
        loadDefaultLists();
        loadUserLists();
    }


    public static GeneList getGeneList(String listID) {
        return geneLists.get(listID);
    }

    public static LinkedHashMap<String, GeneList> getGeneLists() {
        return geneLists;
    }
    // Gene lists -- these don't belong here obviously

    /**
     * Import a .gmt file
     *
     * @param file
     */
    public static void importGMT(File file) {

    }


    public static void addGeneList(GeneList genes) {
        geneLists.put(genes.getName(), genes);
        groups.add(genes.getGroup());
    }


    private static void loadDefaultLists() {

        for (String geneListFile : DEFAULT_GENE_LISTS) {
            InputStream is = GeneListManager.class.getResourceAsStream(geneListFile);
            if (is == null) {
                log.info("Could not find gene list resource: " + geneListFile);
                return;
            }
            BufferedReader reader = null;

            try {
                reader = new BufferedReader(new InputStreamReader(is));
                new BufferedReader(new InputStreamReader(is));
                List<GeneList> lists = loadGMT(geneListFile, reader);
                for (GeneList gl : lists) {
                    gl.setEditable(false);
                    addGeneList(gl);

                }
            } catch (IOException e) {
                log.error("Error loading default gene lists", e);
            } finally {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
        }


    }

    private static void loadUserLists() {
        File dir = Globals.getGeneListDirectory();
        if (dir.exists()) {
            for (File f : dir.listFiles()) {
                GeneList geneList = loadGRPFile(f);
                geneList.setGroup("My lists");
                if (geneList != null) {
                    addGeneList(geneList);
                }
            }

        }
    }

    static GeneList loadGRPFile(File f) {
        String name = f.getName();
        String group = "My lists";
        String description = null;
        List<String> genes = new ArrayList();
        BufferedReader reader = null;

        try {
            reader = new BufferedReader(new FileReader(f));
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("#")) {

                    if (nextLine.startsWith("#name")) {
                        String[] tokens = nextLine.split("=");
                        if (tokens.length > 1) {
                            name = tokens[1];
                        }
                    } else if (nextLine.startsWith("#group")) {

                    } else if (nextLine.startsWith("#description")) {
                        String[] tokens = nextLine.split("=");
                        if (tokens.length > 1) {
                            group = tokens[1];
                        }
                        description = tokens[1];
                    }
                } else {
                    String[] tokens = nextLine.split("\\s+");
                    for (String s : tokens) {
                        genes.add(s);
                    }
                }
            }
            if (genes.size() > 0) {
                if (name == null) {
                    name = f.getName();
                }
                return new GeneList(name, description, group, genes);

            }
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
        return null;
    }

    static List<GeneList> loadGMTFile(File f) {

        String group = f.getName();
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(f));
            return loadGMT(group, reader);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
        return null;
    }

    private static List<GeneList> loadGMT(String group, BufferedReader reader) throws IOException {
        String nextLine;
        List<GeneList> lists = new ArrayList();
        while ((nextLine = reader.readLine()) != null) {
            if (nextLine.startsWith("#")) {
                if (nextLine.startsWith("#group") || nextLine.startsWith("#name")) {
                    String[] tokens = nextLine.split("=");
                    if (tokens.length > 1) {
                        group = tokens[1];
                    }
                }
            } else {
                String[] tokens = nextLine.split("\t");
                String name = tokens[0];
                String description = tokens[1].replaceFirst(">", "");
                // TODO -- description = tokes[1];
                List<String> genes = new ArrayList();
                for (int i = 2; i < tokens.length; i++) {
                    genes.add(tokens[i]);
                }
                lists.add(new GeneList(name, description, group, genes));

            }
        }
        return lists;
    }
}
