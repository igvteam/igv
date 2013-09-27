/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.plugin.mongocollab;

import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.Mongo;
import com.mongodb.WriteResult;
import org.apache.log4j.Logger;
import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.dev.api.LoadHandler;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.UnknownHostException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jacob
 * @date 2013-Sep-06
 */
public class MongoCollabPlugin implements IGVPlugin {

    private static Logger log = Logger.getLogger(MongoCollabPlugin.class);

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
                if(err == null) count += 1;
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
                connection = new Mongo(host, port);
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
        public final boolean buildIndex;

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
            if(fields.containsKey("buildIndex")){
                tmpBuildIndex = Boolean.parseBoolean(fields.get("buildIndex"));
            }
            this.buildIndex = tmpBuildIndex;

        }

    }

    public static class TrackLoadHandler implements LoadHandler {

        @Override
        public void load(String path, List<Track> newTracks) throws IOException{
            Locator locator = new Locator(path);
            MongoFeatureSource.loadFeatureTrack(locator, newTracks);
        }
    }
}
