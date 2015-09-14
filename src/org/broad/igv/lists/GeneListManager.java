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

package org.broad.igv.lists;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.track.TrackProperties;
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

    private static final HashSet<String> fileTypes = new HashSet(Arrays.asList("bed", "gmt", "grp"));

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
                List<GeneList> lists = loadGMTFile(reader);
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
                    if (fileTypes.contains(getFileType(f.getPath()))) {
                        importFile(f);
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

    private String getFileType(String path) {
        String tmp = path.toLowerCase();
        if (tmp.endsWith(".gz")) {
            tmp = tmp.substring(0, tmp.length() - 3);
        }
        int idx = path.lastIndexOf(".");
        return idx < 0 ? path : path.substring(idx + 1);
    }

    public List<GeneList> importFile(File f) throws IOException {

        String path = f.getPath();
        String name = f.getName();

        File dir = DirectoryManager.getGeneListDirectory();
        if (!dir.equals(f.getParentFile())) {
            File copy = new File(dir, f.getName());
            FileUtils.copyFile(f, copy);
        }

        List<GeneList> loadedLists = loadFile(path);

        if (loadedLists.size() > 0) {
            importedFiles.put(name, f);
        }

        return loadedLists;

    }

    private List<GeneList> loadFile(String path) throws IOException {

        String type = getFileType(path);
        if (type.equals("bed")) {
            return loadBEDFile(path);
        } else if (type.equals("gmt")) {
            return loadGMTFile(path);
        } else if (type.equals("grp")) {
            return loadGRPFile(path);
        } else {
            throw new RuntimeException("Unrecognized file extension: " + path);
        }
    }


    private List<GeneList> loadBEDFile(String path) throws IOException {

        String name = (new File(path)).getName();

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(path);

            try {
                List<String> loci = new ArrayList<String>(1000);
                String nextLine;
                while ((nextLine = reader.readLine()) != null) {

                    if (nextLine.startsWith("#") || nextLine.startsWith("browser"))
                        continue;
                    if (nextLine.startsWith("track")) {
                        TrackProperties tp = new TrackProperties();
                        ParsingUtils.parseTrackLine(nextLine, tp);
                        String tmp = tp.getName();
                        if (tmp != null) name = tmp;
                        continue;
                    }
                    String[] tokens = Globals.whitespacePattern.split(nextLine);
                    if (tokens.length > 2) {
                        loci.add(tokens[0] + ":" + tokens[1] + "-" + tokens[2]);
                    }
                }
                GeneList geneList = new GeneList(name, loci);
                geneList.setGroup(USER_GROUP);
                geneList.setEditable(false);
                addGeneList(geneList);

                return Arrays.asList(geneList);

            } finally {
                if (reader != null) reader.close();
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

    public List<GeneList> loadGMTFile(String path) throws IOException {

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(path);
            return loadGMTFile(reader);
        } finally {
            if (reader != null) reader.close();
        }
    }

    private List<GeneList> loadGMTFile(BufferedReader reader) throws IOException {

        String group = USER_GROUP;

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

        for (GeneList gl : lists) {
            gl.setEditable(false);
            addGeneList(gl);
        }
        return lists;
    }

    List<GeneList> loadGRPFile(String path) throws IOException {

        String name = (new File(path)).getName();
        String group = USER_GROUP;
        String description = null;
        //
        BufferedReader reader = null;

        try {
            List<String> genes = new ArrayList();
            reader = ParsingUtils.openBufferedReader(path);
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
                GeneList geneList = new GeneList(name, description, group, genes);
                geneList.setEditable(false);
                addGeneList(geneList);
                return Arrays.asList(geneList);
            } else {
                return Collections.EMPTY_LIST;
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
        if (f != null && f.exists()) {
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
