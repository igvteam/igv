/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.plugin.mongocollab;

import com.mongodb.BasicDBObject;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Dec-14
 */
public class MongoFeatureSource implements FeatureSource<DBFeature.IGVFeat> {

    private int featureWindowSize = 1000000;

    private DBCollection collection;
    private boolean hasLocusIndex = false;


    public MongoFeatureSource(DBCollection collection, boolean buildIndex) {
        this.collection = collection;
        this.collection.setObjectClass(DBFeature.class);
        checkForLocusIndex(buildIndex);
    }

    boolean hasLocusIndex(){
        return this.hasLocusIndex;
    }

    /**
     * Check to see if we have an index useful for queries
     * @param buildIndex Whether to build locus index if not found
     */
    private void checkForLocusIndex(boolean buildIndex){
        if(buildIndex){
            ensureLocusIndex(collection);
        }
        //Check to see if we have the index we want
        List<DBObject> indexes = collection.getIndexInfo();
        DBObject neededFields = getLocusIndexKeys();
        for(DBObject index: indexes){

            boolean isMatchingIndex = true;
            DBObject indexKey = (DBObject) index.get("key");
            for(String key: neededFields.keySet()){
                boolean hasKey = indexKey.containsField(key);
                if(!hasKey){
                    isMatchingIndex = false;
                    break;
                }
                Object value = indexKey.get(key);
                boolean equals = neededFields.get(key).equals(value);
                isMatchingIndex &= equals;
            }

            if(isMatchingIndex) {
                this.hasLocusIndex = true;
                break;
            }
        }
    }


    private DBObject createQueryObject(String chr, int start, int end){
        BasicDBObject query = new BasicDBObject("Chr", chr);

        //Only query over given interval
        //See http://docs.mongodb.org/manual/tutorial/query-documents/
        query.append("Start", new BasicDBObject("$lte", end));
        query.append("End", new BasicDBObject("$gte", start));

        return query;
    }

    /**
     * Return the key/value pairs for indexing
     * Store these as doubles for easier comparison, MongoDB stores everything
     * as a double in the DB
     * @return
     */
    private DBObject getLocusIndexKeys(){
        BasicDBObject indexKeys = new BasicDBObject("Chr", 1.0d);
        indexKeys.append("Start", 1.0d);
        indexKeys.append("End", 1.0d);
        return indexKeys;
    }

    /**
     * Ensures there is an index on Chr, Start, and End
     * To speed queries
     * See http://docs.mongodb.org/manual/reference/method/db.collection.ensureIndex/#db.collection.ensureIndex
     * @param collection
     */
    private void ensureLocusIndex(DBCollection collection){
        collection.ensureIndex(getLocusIndexKeys());
    }

    @Override
    public Iterator<DBFeature.IGVFeat> getFeatures(String chr, int start, int end) throws IOException {
         return getFeatures(createQueryObject(chr, start, end));
    }

    public Iterator<DBFeature.IGVFeat> getFeatures(String name) throws IOException {
        BasicDBObject dbObj = new BasicDBObject("UpperName", name.toUpperCase());
        return getFeatures(dbObj);
    }

    private Iterator<DBFeature.IGVFeat> getFeatures(DBObject queryObject) throws IOException{
        DBCursor cursor = this.collection.find(queryObject);

        //Sort by increasing start value
        //Only do this if we have an index, otherwise might be too memory intensive
        if(hasLocusIndex){
            cursor.sort(new BasicDBObject("Start", 1));
        }
        boolean isSorted = true;
        int lastStart = -1;

        List<DBFeature.IGVFeat> features = new ArrayList<DBFeature.IGVFeat>();
        while (cursor.hasNext()) {
            DBObject obj = cursor.next();
            DBFeature feat = (DBFeature) obj;
            features.add(feat.createIGVFeature());
            isSorted &= feat.getStart() >= lastStart;
            lastStart = feat.getStart();
        }

        if(!isSorted){
            FeatureUtils.sortFeatureList(features);
        }

        return features.iterator();
    }


    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    @Override
    public int getFeatureWindowSize() {
        return featureWindowSize;
    }

    @Override
    public void setFeatureWindowSize(int size) {
        this.featureWindowSize = size;
    }

    DBCollection getCollection(){
        return this.collection;
    }

    public static FeatureTrack loadFeatureTrack(MongoCollabPlugin.Locator locator, List<Track> newTracks) {

        DBCollection collection = MongoCollabPlugin.getCollection(locator);
        //TODO Make this more flexible
        collection.setObjectClass(DBFeature.class);
        MongoFeatureSource source = new MongoFeatureSource(collection, locator.buildIndex);
        FeatureTrack track = new MongoFeatureTrack(collection.getFullName(), collection.getName(), source);
        newTracks.add(track);
        track.setMargin(0);
        return track;
    }
}
