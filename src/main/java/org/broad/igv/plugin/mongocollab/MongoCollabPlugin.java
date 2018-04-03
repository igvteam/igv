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

package org.broad.igv.plugin.mongocollab;

import com.mongodb.*;
import org.apache.log4j.Logger;
import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.dev.api.LoadHandler;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.bson.BSON;
import org.bson.Transformer;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.UnknownHostException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author jacob
 * @date 2013-Sep-06
 */
public class MongoCollabPlugin implements IGVPlugin {

    private static Logger log = Logger.getLogger(MongoCollabPlugin.class);

    static{
        BSON.addEncodingHook(Color.class, new Transformer() {
            @Override
            public Object transform(Object o) {
                if(o instanceof Color){
                    return ColorUtilities.colorToString((Color) o);
                }else{
                    return o;
                }
            }
        });
    }

    @Override
    public void init() {
        TrackLoader.registerHandler("db.spec", new TrackLoadHandler());
    }

    static void insertFeaturesFromFile(DBCollection collection){

        File featFile = FileDialogUtils.chooseFile("Select feature file");
        if (featFile != null) {
            String path = featFile.getAbsolutePath();
            insertFeaturesFromFile(collection, path);
        }
    }

    /**
     *
     * @param collection
     * @param path
     * @return The number of features inserted
     */
    static int insertFeaturesFromFile(DBCollection collection, String path) {
        FeatureCodec codec = CodecFactory.getCodec(path, GenomeManager.getInstance().getCurrentGenome());

        if (codec != null) {
            AbstractFeatureReader<Feature, ?> bfs = AbstractFeatureReader.getFeatureReader(path, codec, false);

            Iterable<Feature> iter = null;
            try {
                iter = bfs.iterator();
            } catch (IOException ex) {
                log.error(ex.getMessage(), ex);
                throw new RuntimeException("Error reading file: " + path, ex);
            }

            int count = 0;
            for (Feature feat : iter) {
                String err = saveFeature(collection, DBFeature.create(feat));
                if(err == null){
                    count += 1;
                }else{
                    log.error("Error inserting feature: " + err);
                }
            }
            return count;
        } else {
            throw new RuntimeException("Cannot load features from file of this type");
        }
    }


    /**
     * Save the specified FeatDBObject to the specified collection
     * Does either an insert or update
     *
     * @param collection
     * @param dbFeat
     * @return
     */
    static String saveFeature(DBCollection collection, DBFeature dbFeat) {

        String errorMessage = "";
        try {

            if(log.isDebugEnabled()){
                log.debug("Saving feature " + Locus.getFormattedLocusString(dbFeat.getChr(), dbFeat.getStart(), dbFeat.getEnd()));
            }
            WriteResult wr = collection.save(dbFeat);
            errorMessage = wr.getError();

        } catch (Exception ex) {
            errorMessage = ex.getMessage();
            log.error(errorMessage, ex);
            if (errorMessage == null) errorMessage = "" + ex;
        }
        return errorMessage;
    }

    public static void removeFeature(DBCollection collection, DBFeature featDBObject) {
        collection.remove(featDBObject);
    }

    private static Map<String, Mongo> connections = new HashMap<String, Mongo>();

    static Mongo getMongo(String host, int port){
        String key = createConnString(host, port);
        Mongo connection = connections.get(key);
        if(connection == null){
            try {
                log.info("Connecting to MongoDB host=" + host + " port=" + port);
                connection = new MongoClient(host, port);
            } catch (UnknownHostException e) {
                log.error(e.getMessage(), e);
                throw new RuntimeException(e.getMessage(), e);
            }
            connections.put(key, connection);
        }
        return connection;
    }


    static void closeMongo(String host, int port) {
        String key = createConnString(host, port);
        Mongo connection = connections.get(key);
        if(connection != null){
            log.info("Closing connection to MongoDB host=" + host + " port=" + port);
            connection.close();
            connections.remove(key);
        }
    }

    private static String createConnString(String host, int port){
        return host + ":" + port;
    }

    static DBCollection getCollection(Locator locator) {
        Mongo mongo = getMongo(locator.host, locator.port);
        DB mongoDB = mongo.getDB(locator.dbName);
        return mongoDB.getCollection(locator.collectionName);
    }

    public static class Locator {
        public final String host;
        public final int port;
        public final String dbName;
        public final String collectionName;
        public final boolean buildLocusIndex;

        public Locator(String path) throws IOException{
            this(ParsingUtils.openInputStream(path));
        }

        public Locator(InputStream is){
            Map<String, String> fields = ParsingUtils.loadMap(is);
            this.host = fields.get("host");
            this.port = Integer.parseInt(fields.get("port"));
            this.dbName = fields.get("dbName");
            this.collectionName = fields.get("collectionName");
            boolean tmpBuildIndex = false;
            if(fields.containsKey("buildLocusIndex")){
                tmpBuildIndex = Boolean.parseBoolean(fields.get("buildLocusIndex"));
            }
            this.buildLocusIndex = tmpBuildIndex;

        }

    }

    public static final int DB_EXISTS = 0x1;
    public static final int COLLECTION_EXISTS = 0x2;

    /**
     * Check whether the database/collection specified by the given {@code locator}
     * exists.
     * @param locator
     * @return An integer consisting of 2 flags:
     * DB_EXISTS
     * COLLECTION_EXISTS
     *
     * with each set if the DB/COLLECTION exists. !DB_EXISTS && COLLECTION_EXISTS should never happen
     */
    public static int checkDestinationExists(Locator locator){
        Mongo mongo = getMongo(locator.host, locator.port);
        List<String> dbNames = mongo.getDatabaseNames();
        boolean dbExists = dbNames.indexOf(locator.dbName) >= 0;

        if(!dbExists){
            return 0;
        }

        int result = DB_EXISTS;
        DB db = mongo.getDB(locator.dbName);
        Set<String> collections = db.getCollectionNames();
        result |= collections.contains(locator.collectionName) ? COLLECTION_EXISTS : 0;

        return result;
    }

    public static class TrackLoadHandler implements LoadHandler {

        @Override
        public void load(String path, List<Track> newTracks) throws IOException{
            Locator locator = new Locator(path);

            int destExists = checkDestinationExists(locator);
            String toAsk = null;
            boolean doLoadTrack = false;

            if( (destExists & COLLECTION_EXISTS) > 0){
                //Both DB and collection exist
                assert (destExists & DB_EXISTS) > 0;
                doLoadTrack = true;
            }else if(destExists == 0){
                //Neither collection nor track exist
                toAsk = String.format("Host '%s' does not contain database '%s'. Do you wish to create it?\n", locator.host, locator.dbName);
                toAsk += String.format("If you select yes, collection '%s' will be created as well.", locator.collectionName);
            }else{
                //Collection doesn't exist but DB does
                toAsk = String.format("Host '%s', database '%s', does not contain collection '%s'.\nDo you wish to create it?",
                        locator.host, locator.dbName, locator.collectionName);
            }

            if(toAsk != null){
                doLoadTrack = MessageUtils.confirm(toAsk);
            }

            if(doLoadTrack) MongoFeatureSource.loadFeatureTrack(locator, newTracks);
        }
    }

}
