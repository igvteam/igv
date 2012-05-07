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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.maf;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 */
public class MAFManager {

    private static Logger log = Logger.getLogger(MAFManager.class);
    public static String[] species = {
            "panTro2", "gorGor1", "ponAbe2", "rheMac2", "calJac1",
            "tarSyr1", "micMur1", "otoGar1", "tupBel1", "mm9", "rn4", "dipOrd1",
            "cavPor3", "speTri1", "oryCun1", "ochPri2", "vicPac1", "turTru1",
            "bosTau4", "equCab2", "felCat3", "canFam2", "myoLuc1", "pteVam1",
            "eriEur1", "sorAra1", "loxAfr2", "proCap1", "echTel1", "dasNov2",
            "choHof1", "monDom4", "ornAna1", "galGal3", "taeGut1", "anoCar1",
            "xenTro2", "tetNig1", "fr2", "gasAcu1", "oryLat28", "danRer5", "petMar1"
    };
    private List<String> selectedSpecies;
    private List<String> allSpecies;
    private Map<String, String> speciesNames;

    private int tileSize = 500;
    LRUCache<String, MAFTile> tileCache;

    MAFReader reader;
    String refId;
    Genome genome;


    public MAFManager(final ResourceLocator locator, Genome genome) throws IOException {
        tileCache = new LRUCache(this, 10);
        loadSpecies(locator.getPath());
        speciesNames.put(genome.getId(), genome.getDisplayName());
        overrides();

        if (locator.getPath().endsWith(".maf.dict")) {
            reader = new MAFListReader(locator.getPath());
        } else {
            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        reader = new MAFLocalReader(locator.getPath());
                        IGV.getInstance().repaintDataAndHeaderPanels();
                    } catch (Exception e) {
                        log.error("Error loading MAF reader (" + locator.getPath() + "):  ", e);
                        MessageUtils.showMessage("Error loading MAF file: " + e.getMessage());
                    }
                }
            };
            LongRunningTask.submit(runnable);
        }
        //.asList(new ArrayList<String>(), species);
    }

    private void overrides() {
        String[] humanIds = {"hg18", "1kg_ref", "hg19", "b37", "1kg_v37"};
        for (String id : humanIds) {
            speciesNames.put(id, "Human");
        }
    }

    public List<String> getSelectedSpecies() {
        return selectedSpecies;
    }

    /**
     * @param selectedSpecies the selectedSpecies to set
     */
    public void setSelectedSpecies(List<String> selectedSpecies) {
        this.selectedSpecies = selectedSpecies;
        this.allSpecies = new ArrayList(selectedSpecies);
        allSpecies.add(0, refId);
    }

    /**
     * Load a file containing species IDs and names.
     *
     * @param path
     */
    private void loadSpecies(String path) {

        InputStream is = null;
        speciesNames = new HashMap<String, String>();
        List<String> species = new ArrayList<String>();

        try {
            String speciesPath = path + ".species";
            if (FileUtils.resourceExists(speciesPath)) {
                is = ParsingUtils.openInputStream(speciesPath);
            } else {
                is = MAFUtils.class.getResourceAsStream("species.properties");
            }

            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("#ref")) {
                    String[] tokens = Globals.equalPattern.split(nextLine);
                    refId = tokens[1];
                } else {
                    String[] tokens = Globals.equalPattern.split(nextLine);
                    if (tokens.length == 2) {
                        String id = tokens[0];
                        String name = tokens[1];
                        speciesNames.put(id, name);
                        species.add(id);
                    } else {
                        //log.info("Skipping line: " + nextLine);
                    }
                }
            }
            setSelectedSpecies(species);

        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {

                }
            }
        }
    }

    public String getSpeciesName(String speciesId) {
        String name = speciesNames.get(speciesId);
        return name == null ? speciesId : name;
    }

    public List<String> getChrNames() {
        return reader.getChrNames();
    }

    public MAFTile[] getTiles(String chr, int start, int end) {

        int startTile = start / tileSize;
        int endTile = end / tileSize;
        MAFTile[] tiles = new MAFTile[endTile - startTile + 1];
        for (int i = startTile; i <= endTile; i++) {
            tiles[i - startTile] = getTile(chr, i);
        }
        return tiles;
    }

    private MAFTile getTile(String chr, int tileNo) {
        if (reader == null) {
            // This will be true while the index is loaded (asynchronously)
            return null;
        }

        String key = getKey(chr, tileNo);
        MAFTile tile = tileCache.get(key);
        if (tile == null) {
            int start = tileNo * tileSize;
            int end = start + tileSize + 1;
            tile = reader.loadTile(chr, start, end, allSpecies);
            tileCache.put(key, tile);
        }
        return tile;
    }

    static String getKey(String chr, int tileNo) {
        return chr + tileNo;
    }

}
