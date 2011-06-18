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
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.util.LRUCache;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

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
    public static Properties speciesNames;
    private int tileSize = 500;
    LRUCache<String, MAFTile> tileCache;
    private List<String> selectedSpecies;
    private List<String> allSpecies;
    MAFReader reader;
    static MAFManager instance;
    String refId;

    public MAFManager(ResourceLocator locator) {
        tileCache = new LRUCache(this, 10);
        if (speciesNames == null) {
            loadNames();
        }

        if (locator.isLocal()) {
            reader = new MAFLocalReader(locator.getPath());
            allSpecies = ((MAFLocalReader) reader).getSequenceIds();
            selectedSpecies = allSpecies.subList(1, allSpecies.size());
            refId = allSpecies.get(0);
        } else {
            reader = new MAFRemoteReader(locator); // MAFLocalReader();
            this.setSelectedSpecies(PreferenceManager.getInstance().getMafSpecies());
            refId = "hg18";
        }
        //.asList(new ArrayList<String>(), species);
    }

    private void loadNames() {
        try {
            InputStream is = MAFUtils.class.getResourceAsStream("species.properties");
            speciesNames = new Properties();
            speciesNames.load(is);
            is.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    public String getSpeciesName(String speciesId) {
        return speciesNames.getProperty(speciesId);
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
        String key = getKey(chr, tileNo);
        MAFTile tile = tileCache.get(key);
        if (tile == null) {
            WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
            try {

                int start = tileNo * tileSize;
                int end = start + tileSize + 1;
                tile = reader.loadTile(chr, start, end, allSpecies);
                tileCache.put(key, tile);
            } finally {
                WaitCursorManager.removeWaitCursor(token);
            }

        }
        return tile;
    }

    static String getKey(String chr, int tileNo) {
        return chr + tileNo;
    }

    /**
     * @return the selectedSpecies
     */
    public List<String> getSelectedSpecies() {
        return selectedSpecies;
    }

    /**
     * @param selectedSpecies the selectedSpecies to set
     */
    public void setSelectedSpecies(List<String> selectedSpecies) {
        this.selectedSpecies = selectedSpecies;
        this.allSpecies = new ArrayList(selectedSpecies);
        allSpecies.add(0, "hg18");
    }
}
