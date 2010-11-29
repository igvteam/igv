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

import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 26, 2010
 */
public class GeneListManager {

    private static LinkedHashMap<String, GeneList> geneLists = new LinkedHashMap();


    public static GeneList getGeneList(String listID) {
        return geneLists.get(listID);
    }

    public static LinkedHashMap<String, GeneList> getGeneLists() {
        return geneLists;
    }
    // Gene lists -- these don't belong here obviously

    static {

        geneLists.put("None", null);

        //chr6:63,599,345-171,829,806
        addNewGeneList(new GeneList("Proneural dev genes", Arrays.asList("SOX1", "SOX2", "SOX3", "SOX21", "DCX",
                "DLL3", "ASCL1", "TCF4")));
        addNewGeneList(new GeneList("Microglia markers", Arrays.asList("CD68", "PTPRC", "TNF")));
        addNewGeneList(new GeneList("MMR genes", Arrays.asList("MGMT", "MLH1", "MSH2", "MSH6", "PMS2")));
        addNewGeneList(new GeneList("PI(3)K complex", Arrays.asList("PIK3CA", "PIK3R1")));
        addNewGeneList(new GeneList("Oligodendrocytic dev genes", Arrays.asList("PDGFRA", "NKX2-2", "NG2", "OLIG2",
                "CDKN1A")));
        addNewGeneList(new GeneList("P53 signaling", Arrays.asList("CDKN2A", "TP53", "MDM2", "MDM4")));
        addNewGeneList(new GeneList("RB signaling", Arrays.asList("CDKN2A", "CDKN2B", "CDKN2C", "CDK4", "CDK6",
                "CCND1", "CCND2", "RB1")));
        addNewGeneList(new GeneList("AKT signaling", Arrays.asList("PIK3CA", "PTEN", "PIP3", "AKT1", "AKT2", "AKT3")));
        addNewGeneList(new GeneList("RTK/RAS signaling", Arrays.asList("PDGFRA", "PDGFRB", "ERBB3", "EGFR", "ERBB2",
                "FGFR1", "FGFR2", "MET", "GRB2", "NF1", "NRAS", "KRAS", "HRAS", "ARAF", "BRAF", "RAF1", "SPRY2")));

        loadGeneLists();


    }


    public static void addNewGeneList(GeneList genes) {
        geneLists.put(genes.getName(), genes);
    }

    

    private static void loadGeneLists() {
         File dir = Globals.getGeneListDirectory();
         if (dir.exists()) {
             for (File f : dir.listFiles()) {
                 String name = null;
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
                             }
                         } else {
                             String[] tokens = nextLine.split("\\s+");
                             for (String s : tokens) {
                                 genes.add(s);
                             }
                         }
                         if(genes.size() > 0) {
                             if(name == null) {
                                 name = f.getName();
                             }
                             geneLists.put(name, new GeneList(name, genes));
                         }
                     }
                 } catch (IOException e) {
                     e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                 } finally {
                     if(reader != null) {
                         try {
                             reader.close();
                         } catch (IOException e) {
                             e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                         }
                     }
                 }
             }

         }
     }


}
