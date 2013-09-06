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
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.ui.IGV;
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
import java.util.HashMap;
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


        JMenu collabPluginItem = new JMenu("Annotate via DB");
        collabPluginItem.add(insertFeattoDBItem);
        //collabPluginItem.add(loadAnnotTrack);
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
        BasicDBObject obj = new BasicDBObject();
        obj.put("chr", feat.getChr());
        obj.put("start", feat.getStart());
        obj.put("end", feat.getEnd());
        if (feat instanceof IGVFeature) {
            obj.put("description", ((IGVFeature) feat).getDescription());
        }
        return obj;
    }


    private class FeatDBObject extends ReflectionDBObject{

        public String chr;
        public int start;
        public int end;
        public String description;

        FeatDBObject(Feature feat){
            this.chr = feat.getChr();
            this.start = feat.getStart();
            this.end = feat.getEnd();
            if(feat instanceof IGVFeature){
                this.description = ((IGVFeature) feat).getDescription();
            }
        }

    }
}
