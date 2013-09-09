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

import com.mongodb.*;
import org.apache.log4j.Logger;
import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.feature.AbstractFeature;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.PanelName;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jacob
 * @date 2013-Sep-06
 */
public class MongoCollabPlugin implements IGVPlugin {

    private static Logger log = Logger.getLogger(MongoCollabPlugin.class);

    //TODO Make into config file
    private static String host = "localhost";
    private static int port = 27017;
    private static String dbname = "mylocaldb";
    private static String collectionName = "annotations";

    @Override
    public void init() {

        JMenuItem insertFeattoDBItem = new JMenuItem("Load features into DB");
        insertFeattoDBItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                File featFile = FileDialogUtils.chooseFile("Select feature file");
                if(featFile != null){
                    String path = featFile.getAbsolutePath();
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
                        DBCollection collection = getCollection(host, port, dbname, collectionName);
                        //TODO Make this more flexible
                        collection.setObjectClass(FeatDBObject.class);
                        for(Feature feat: iter){
                            DBObject featdbobj = createFeatDBObject(feat);
                            WriteResult result = collection.insert(featdbobj);
                        }

                    }else{
                        throw new RuntimeException("Cannot load features from file of this type");
                    }
                }
            }
        });

        JMenuItem loadAnnotTrack = new JMenuItem("Load Annotation Track From DB");
        loadAnnotTrack.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                //TODO Specify db from config file
                DBCollection collection = getCollection(host, port, dbname, collectionName);
                List<Track> newTracks = new ArrayList<Track>(1);
                MongoFeatureSource.loadFeatureTrack(collection, newTracks);
                IGV.getInstance().addTracks(newTracks, PanelName.DATA_PANEL);
            }
        });



        JMenu collabPluginItem = new JMenu("Annotate via DB");
        collabPluginItem.add(insertFeattoDBItem);
        collabPluginItem.add(loadAnnotTrack);
        IGV.getInstance().addOtherToolMenu(collabPluginItem);

    }

    private static Map<String, Mongo> connections = new HashMap<String, Mongo>();

    private static Mongo getMongo(String host, int port){
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

    private static String createConnString(String host, int port){
        return host + ":" + port;
    }

    private DBCollection getCollection(String host, int port, String dbname, String collectionName) {
        Mongo mongo = getMongo(host, port);
        DB mongoDB = mongo.getDB(dbname);
        return mongoDB.getCollection(collectionName);
    }

    private DBObject createFeatDBObject(Feature feat) {
        return FeatDBObject.create(feat);
//        BasicDBObject obj = new BasicDBObject();
//        obj.put("chr", feat.getChr());
//        obj.put("start", feat.getStart());
//        obj.put("end", feat.getEnd());
//        if (feat instanceof IGVFeature) {
//            obj.put("description", ((IGVFeature) feat).getDescription());
//        }
//        return obj;
    }

    /**
     * Object mapping to Mongo database
     * ReflectionDBObject works with getters/setters, and
     * doesn't use the Java Beans case convention.
     * So (get/set)Chr maps to a field named "Chr", not "chr"
     * as we might prefer
     *
     * TODO Use existing feature interfaces/classes, which are long past
     * overdue for refactoring
     */
    public static class FeatDBObject extends ReflectionDBObject implements Feature {

        private String chr;
        private int start;
        private int end;
        private String description;
        private double score;

        @SubtlyImportant
        public FeatDBObject(){}

        public FeatDBObject(String chr, int start, int end, String description, double score){
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.description = description;
            this.score = score;
        }

        static FeatDBObject create(Feature feature){
            if(feature instanceof AbstractFeature){
                return create((AbstractFeature) feature);
            }
            return new FeatDBObject(feature.getChr(), feature.getStart(), feature.getEnd(), null, 0);
        }

        static FeatDBObject create(AbstractFeature feature){
            return new FeatDBObject(feature.getChr(), feature.getStart(), feature.getEnd(), feature.getDescription(), feature.getScore());
        }

        public String getChr() {
            return chr;
        }

        public void setChr(String chr) {
            this.chr = chr;
        }

        public String getDescription() {
            return description;
        }

        public void setDescription(String description) {
            this.description = description;
        }

        public int getEnd() {
            return end;
        }

        public void setEnd(int end) {
            this.end = end;
        }

        public double getScore() {
            return score;
        }

        public void setScore(double score) {
            this.score = score;
        }

        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }

        public BasicFeature createBasicFeature(){
            BasicFeature bf = new BasicFeature(chr, start, end);
            bf.setDescription(this.description);
            //TODO Shouldn't just cast from double to float
            bf.setScore((float) this.score);
            return bf;
        }


    }
}
