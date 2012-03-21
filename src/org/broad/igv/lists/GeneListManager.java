/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
import org.broad.igv.DirectoryManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.net.URLEncoder;
import java.util.*;

/**
 * @author jrobinso
 * @date Sep 26, 2010
 */
public class GeneListManager {

    private static Logger log = Logger.getLogger(GeneListManager.class);

    public static final List<String> DEFAULT_GENE_LISTS = Arrays.asList(
            /*"biocarta_cancer_cp.gmt",*/  "examples.gmt", "reactome_cp.gmt", "kegg_cancer_cp.gmt");

    public static final String USER_GROUP = "My lists";

    private LinkedHashSet<String> groups = new LinkedHashSet();

    private HashMap<String, File> importedFiles = new HashMap();

    private LinkedHashMap<String, GeneList> geneLists = new LinkedHashMap();

    static GeneListManager theInstance;

    public static GeneListManager getInstance() {
        if (theInstance == null) {
            theInstance = new GeneListManager();
        }
        return theInstance;
    }

    private GeneListManager() {
        loadDefaultLists();
        loadUserLists();
    }


    public GeneList getGeneList(String listID) {
        return geneLists.get(listID);
    }

    public LinkedHashMap<String, GeneList> getGeneLists() {
        return geneLists;
    }
    // Gene lists -- these don't belong here obviously


    public void addGeneList(GeneList genes) {
        geneLists.put(genes.getName(), genes);
        groups.add(genes.getGroup());
    }


    private void loadDefaultLists() {

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
                MessageUtils.showMessage("<html>Error encountered loading gene lists (" + e.toString() + ")" +
                        "<br/>See log for more details");
            } finally {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
        }
    }

    private void loadUserLists() {
        File dir = DirectoryManager.getGeneListDirectory();
        if (dir.exists()) {
            for (File f : dir.listFiles()) {
                try {
                    if (f.getName().toLowerCase().endsWith(".gmt")) {
                        importGMTFile(f);
                    } else {
                        GeneList geneList = loadGRPFile(f);
                        geneList.setEditable(true);
                        geneList.setGroup(USER_GROUP);
                        if (geneList != null) {
                            addGeneList(geneList);
                        }
                    }
                } catch (IOException e) {
                    log.error("Error loading user gene lists: ", e);
                    MessageUtils.showMessage("<html>Error encountered loading user gene lists (" + e.toString() + ")" +
                            "<br/>See log for more details");
                }
            }

        }

        // Add empty group if there are no lists
        if (!groups.contains(USER_GROUP)) {
            groups.add(USER_GROUP);
        }
    }

    GeneList loadGRPFile(File grpFile) throws IOException {

        // First copy file to gene list directory
        File f = grpFile;
        File dir = DirectoryManager.getGeneListDirectory();
        if (!dir.equals(grpFile.getParentFile())) {
            f = new File(dir, grpFile.getName());
            FileUtils.copyFile(grpFile, f);
        }

        String name = f.getName();
        String group = USER_GROUP;
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
                    } else if (nextLine.startsWith("#description")) {
                        String[] tokens = nextLine.split("=");
                        if (tokens.length > 1) {
                            description = tokens[1];
                        }

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
                importedFiles.put(name, f);
                return new GeneList(name, description, group, genes);
            }
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
        }
        return null;
    }

    // TODO -- why are there 2 of these methods?
    public void importGMTFile(File gmtFile) throws IOException {

        File f = gmtFile;
        File dir = DirectoryManager.getGeneListDirectory();
        if (!dir.equals(gmtFile.getParentFile())) {
            f = new File(dir, gmtFile.getName());
            FileUtils.copyFile(gmtFile, f);
        }

        String group = f.getName().replace(".gmt", "");
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(f));
            importedFiles.put(group, f);
            List<GeneList> lists = loadGMT(group, reader);
            for (GeneList gl : lists) {
                gl.setEditable(true);
                addGeneList(gl);
            }
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
        }
    }


    // TODO -- why are there 2 of these methods?
    public List<GeneList> importGMTFile(String path) throws IOException {

        BufferedReader reader = null;
        try {
            String group = new File(path).getName();
            reader = ParsingUtils.openBufferedReader(path);
            List<GeneList> lists = loadGMT(group, reader);
            for (GeneList gl : lists) {
                gl.setEditable(false);
                addGeneList(gl);
            }
            return lists;

        } finally {
            if (reader != null) reader.close();
        }

    }

    private List<GeneList> loadGMT(String group, BufferedReader reader) throws IOException {
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
                if (tokens.length > 2) {
                    String name = tokens[0];
                    String description = tokens[1].replaceFirst(">", "");
                    List<String> genes = new ArrayList();
                    for (int i = 2; i < tokens.length; i++) {
                        genes.add(tokens[i]);
                    }
                    lists.add(new GeneList(name, description, group, genes));
                }
            }
        }
        return lists;
    }

    /**
     * #name=Example gene lists
     * Proneural dev genes	Proneural dev genes	SOX1    SOX2	SOX3	SOX21	DCX	DLL3	ASCL1	TCF4
     *
     * @param group
     * @param outputFile
     */
    void exportGMT(String group, File outputFile) throws IOException {

        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            List<GeneList> lists = getListsForGroup(group);
            if (lists.isEmpty()) {
                MessageUtils.showMessage("Nothing to export.");
                return;
            }

            pw.println("#name=" + group);
            for (GeneList gl : lists) {
                pw.print(gl.getName());
                for (String gene : gl.getLoci()) {
                    pw.print("\t");
                    pw.print(gene);
                }
                pw.println();
            }
        } finally {
            if (pw != null) pw.close();
        }

    }

    // TODO -- this is really ineffecient, redesign
    private List<GeneList> getListsForGroup(String group) {
        List<GeneList> list = new ArrayList<GeneList>();
        for (GeneList gl : geneLists.values()) {
            if (gl.getGroup().equals(group)) {
                list.add(gl);
            }
        }
        return list;
    }

    public void saveGeneList(GeneList geneList) {

        File file = null;
        PrintWriter pw = null;
        try {
            final String listName = geneList.getName();
            String description = geneList.getDescription();
            List<String> genes = geneList.getLoci();

            if (listName != null && genes != null) {
                file = new File(DirectoryManager.getGeneListDirectory(), getLegalFilename(listName) + ".grp");
                pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
                pw.println("#name=" + listName);
                if (description != null) pw.println("#description=" + description);
                for (String s : genes) {
                    pw.println(s);
                }
                pw.close();
                importedFiles.put(listName, file);
            }

        } catch (IOException e) {
            if (file != null) {
                MessageUtils.showMessage("Error writing gene list file: " + file.getAbsolutePath() + " " + e.getMessage());
            }
            log.error("Error saving gene list", e);
        } finally {
            if (pw != null) {
                pw.close();
            }
        }
    }


    /**
     * Return a legal filename derived from the input string.
     * todo Move this to a utility class
     *
     * @param s
     * @return
     */
    private static String getLegalFilename(String s) {
        try {
            return URLEncoder.encode(s, "UTF-8");
        } catch (UnsupportedEncodingException e) {
            return s;
        }
    }


    /**
     * Test to see of group was imported by the user  Needed to determine if group can be removed.
     */
    public boolean isImported(String groupName) {
        return importedFiles.containsKey(groupName);
    }

    public LinkedHashSet<String> getGroups() {
        return groups;
    }

    public void deleteGroup(String selectedGroup) {
        File f = importedFiles.get(selectedGroup);
        if (f.exists()) {
            f.delete();
        }
        groups.remove(selectedGroup);
        importedFiles.remove(selectedGroup);

        Collection<GeneList> tmp = new ArrayList(geneLists.values());
        for (GeneList gl : tmp) {
            if (gl.getGroup().equals(selectedGroup)) {
                geneLists.remove(gl.getName());
            }
        }
    }


    /**
     * Return true if a group was also removed.
     *
     * @param listName
     */
    public boolean deleteList(String listName) {


        File f = importedFiles.get(listName);
        if (f.exists()) {
            f.delete();
        }
        importedFiles.remove(listName);

        if (geneLists.containsKey(listName)) {
            String group = geneLists.get(listName).getGroup();
            geneLists.remove(listName);

            // If the group is empty remove it as well, except for user group
            if (!group.equals(USER_GROUP)) {
                for (GeneList gl : geneLists.values()) {
                    if (gl.getGroup().equals(group)) {
                        return false;
                    }
                }
                groups.remove(group);
                return true;
            }
        }
        return false;
    }
}
