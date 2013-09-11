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
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuItemBuilder;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.PanelName;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.UnknownHostException;
import java.util.*;

/**
 * @author jacob
 * @date 2013-Sep-06
 */
public class MongoCollabPlugin implements IGVPlugin {

    private static Logger log = Logger.getLogger(MongoCollabPlugin.class);

    @Override
    public void init() {

        JMenuItem loadAnnotTrack = new JMenuItem("Load Annotation Track From DB");
        loadAnnotTrack.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                File locatorFile = FileDialogUtils.chooseFile("Select config file");
                if(locatorFile == null) return;

                Locator locator = null;
                try {
                    locator = new Locator(locatorFile);
                } catch (FileNotFoundException ex) {
                    //This shouldn't happen, user just picked a file
                    log.error(ex.getMessage(), ex);
                    return;
                }
                DBCollection collection = getCollection(locator);
                //TODO Make this more flexible
                collection.setObjectClass(FeatDBObject.class);
                List<Track> newTracks = new ArrayList<Track>(1);
                MongoFeatureSource.loadFeatureTrack(collection, newTracks);
                IGV.getInstance().addTracks(newTracks, PanelName.FEATURE_PANEL);

                //Add context menu entry
                TrackMenuUtils.addTrackMenuItemBuilder(createBuilderForLocator(locator));
            }
        });

        JMenu collabPluginItem = new JMenu("Annotate via DB");
        collabPluginItem.add(loadAnnotTrack);
        //Add "Tools" menu option for loading annotation DB track
        IGV.getInstance().addOtherToolMenu(collabPluginItem);
    }

    private TrackMenuItemBuilder createBuilderForLocator(final Locator locator){
        TrackMenuItemBuilder builder = new TrackMenuItemBuilder() {
            @Override
            public JMenuItem build(Collection<Track> selectedTracks, TrackClickEvent te) {
                ReferenceFrame frame = te.getFrame();
                boolean hasFrame = frame != null;
                if(!hasFrame) return null;

                final String chr = frame.getChrName();
                final int start = (int) te.getChromosomePosition();
                final int end = (int) Math.ceil(frame.getChromosomePosition(te.getMouseEvent().getX() + 1));

                //TODO create dialog so user can fill in score/description
                JMenuItem item = new JMenuItem("Annotate in " + locator.collectionName);
                item.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        DBCollection collection = getCollection(locator);
                        Feature feat = new FeatDBObject(chr, start, end, null, 0.0);
                        insertFeature(collection, feat);
                        log.info("Adding feature " + Locus.getFormattedLocusString(chr, start, end));
                    }
                });
                return item;
            }
        };
        return builder;
    }
    
    private static void insertFeaturesFromFile(DBCollection collection){

        File featFile = FileDialogUtils.chooseFile("Select feature file");
        if (featFile != null) {
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

                for (Feature feat : iter) {
                    insertFeature(collection, feat);
                }

            } else {
                throw new RuntimeException("Cannot load features from file of this type");
            }
        }
    }

    private static WriteResult insertFeature(DBCollection collection, Feature feat) {
        DBObject featdbobj = createFeatDBObject(feat);
        return collection.insert(featdbobj);
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

    private static DBCollection getCollection(Locator locator) {
        Mongo mongo = getMongo(locator.host, locator.port);
        DB mongoDB = mongo.getDB(locator.dbName);
        return mongoDB.getCollection(locator.collectionName);
    }

    private static DBObject createFeatDBObject(Feature feat) {
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
     * TODO overdue for refactoring
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

    public static class Locator {
        public final String host;
        public final int port;
        public final String dbName;
        public final String collectionName;

        public Locator(File file) throws FileNotFoundException{
            this(new FileInputStream(file));
        }

        public Locator(InputStream is){
            Map<String, String> fields = ParsingUtils.loadMap(is);
            this.host = fields.get("host");
            this.port = Integer.parseInt(fields.get("port"));
            this.dbName = fields.get("dbName");
            this.collectionName = fields.get("collectionName");
        }

    }
}
