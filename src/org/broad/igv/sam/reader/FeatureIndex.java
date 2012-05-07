/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
package org.broad.igv.sam.reader;


import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 */
public class FeatureIndex {

    private int tileWidth;
    private LinkedHashMap<String, ChromosomeIndex> chrIndeces;
    private Logger log = Logger.getLogger(FeatureIndex.class);


    /**
     * Constructs ...
     *
     * @param tileWidth
     */
    public FeatureIndex(int tileWidth) {
        this.tileWidth = tileWidth;
        chrIndeces = new LinkedHashMap();

    }

    /**
     */
    public FeatureIndex(File f) {
        this(f.getAbsolutePath());
    }

    /**
     */
    public FeatureIndex(String path) {

        InputStream is = null;
        try {
            is = ParsingUtils.openInputStreamGZ(new ResourceLocator(path));
            chrIndeces = new LinkedHashMap();
            read(is);
        } catch (IOException ex) {
            log.error("Error reading index", ex);
            throw new DataLoadException("Error reading index: " + ex.getMessage(), path);
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }

    }

    public boolean containsChromosome(String chr) {
        return chrIndeces.containsKey(chr);
    }

    public Set<String> getIndexedChromosomes() {
        return chrIndeces.keySet();
    }

    /**
     * @param chr
     * @param idx
     * @param count
     */
    public void add(String chr, long idx, int count, int longestFeature) {

        ChromosomeIndex chrIndex = chrIndeces.get(chr);
        if (chrIndex == null) {
            chrIndex = new ChromosomeIndex(longestFeature);
            chrIndeces.put(chr, chrIndex);
        }
        chrIndex.addTile(new TileDef(idx, count));

    }

    /**
     * @param chr
     * @param tile
     * @return
     */
    public TileDef getTileDef(String chr, int tile) {

        ChromosomeIndex chrIdx = chrIndeces.get(chr);
        if (chrIdx == null) {

            // Todo -- throw an execption ?
            return null;
        } else {
            return chrIdx.getTileDefinition(tile);
        }
    }

    /**
     * Store a SamIndex to a stream.
     * <p/>
     * It is the responsibility of the caller  to close the stream.
     *
     * @param f
     * @throws IOException
     */
    public void store(File f) throws IOException {

        DataOutputStream dos = null;

        try {
            dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f)));

            dos.writeInt(getTileWidth());

            for (Map.Entry<String, ChromosomeIndex> entry : chrIndeces.entrySet()) {

                ChromosomeIndex chrIdx = entry.getValue();
                List<TileDef> tmp = chrIdx.getTileDefinitions();

                if (entry.getKey() != null) {
                    dos.writeUTF(entry.getKey());
                    dos.writeInt(tmp.size());
                    dos.writeInt(chrIdx.getLongestFeature());

                    for (int i = 0; i < tmp.size(); i++) {
                        final FeatureIndex.TileDef tileDef = tmp.get(i);
                        dos.writeLong(tileDef.getStartPosition());
                        dos.writeInt(tileDef.getCount());
                    }
                }
            }
        } finally {
            dos.close();
        }

    }

    /**
     * Read the index.  Translate the chromosome names to the current genome aliases, if any.
     *
     * @throws IOException
     */
    private void read(InputStream is) throws IOException {

        DataInputStream dis = new DataInputStream(new BufferedInputStream(is));

        tileWidth = dis.readInt();
        try {
            while (true) {
                String chr = dis.readUTF();
                int nTiles = dis.readInt();
                int longestFeature = dis.readInt();

                List<TileDef> tileDefs = new ArrayList(nTiles);
                int tileNumber = 0;
                while (tileNumber < nTiles) {
                    long pos = dis.readLong();
                    int count = dis.readInt();
                    tileDefs.add(new TileDef(pos, count));
                    tileNumber++;
                }

                chrIndeces.put(chr, new ChromosomeIndex(longestFeature, tileDefs));
            }
        } catch (EOFException e) {
            // This is normal.  Unfortuantely we don't have a better way to test EOF for this stream
        } catch (UTFDataFormatException e) {
            log.error("Error reading chromosome name. ", e);
        }

    }

    public int getTileWidth() {
        return tileWidth;
    }

    public int getLongestFeature(String chr) {
        ChromosomeIndex tmp = this.chrIndeces.get(chr);
        if (tmp == null) {
            return 1000;
        } else {
            return tmp.getLongestFeature();
        }
    }

    static class ChromosomeIndex {

        private int longestFeature;
        private List<TileDef> tileDefinitions;

        ChromosomeIndex(int longestFeature) {
            this(longestFeature, new ArrayList<TileDef>());
        }

        public ChromosomeIndex(int longestFeature, List<TileDef> tileDefinitions) {
            this.longestFeature = longestFeature;
            this.tileDefinitions = tileDefinitions;
        }

        void addTile(TileDef tileDef) {
            getTileDefinitions().add(tileDef);
        }

        TileDef getTileDefinition(int i) {
            if (getTileDefinitions().isEmpty()) {
                // TODO -- throw exception ?
                return null;
            }
            int tileNumber = Math.min(i, getTileDefinitions().size() - 1);
            return getTileDefinitions().get(tileNumber);
        }

        /**
         * @return the longestFeature
         */
        public int getLongestFeature() {
            return longestFeature;
        }

        /**
         * @return the tileDefinitions
         */
        public List<TileDef> getTileDefinitions() {
            return tileDefinitions;
        }
    }

    public static class TileDef {

        private long startPosition;
        private int count;

        /**
         * Constructs ...
         *
         * @param startPosition
         * @param count
         */
        public TileDef(long startPosition, int count) {
            this.startPosition = startPosition;
            this.count = count;
        }

        /**
         * @return the startPosition
         */
        public long getStartPosition() {
            return startPosition;
        }

        /**
         * @return the count
         */
        public int getCount() {
            return count;
        }
    }
}
