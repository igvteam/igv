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

import com.mongodb.BasicDBObject;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import org.apache.log4j.Logger;
import org.broad.igv.dev.api.NamedFeatureSearcher;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.action.SearchCommand;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Dec-14
 */
public class MongoFeatureSource implements FeatureSource<DBFeature.IGVFeat>, NamedFeatureSearcher {

    private int featureWindowSize = 1000000;

    private DBCollection collection;
    private boolean hasLocusIndex = false;

    private static Logger log = Logger.getLogger(MongoCollabPlugin.class);

    public MongoFeatureSource(DBCollection collection, boolean buildLocusIndex) {
        this.collection = collection;
        this.collection.setObjectClass(DBFeature.class);
        checkForLocusIndex(buildLocusIndex);
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
         return getFeatures(createQueryObject(chr, start, end), 0).iterator();
    }

    @Override
    public Collection<? extends NamedFeature> search(String name, int limit) {
        BasicDBObject dbObj = new BasicDBObject("UpperName", name.toUpperCase());
        try {
            return getFeatures(dbObj, limit);
        } catch (IOException e) {
            log.error(e.getMessage(), e);
            return null;
        }
    }

    /**
     *
     * @param queryObject
     * @param limit Limitation on the number of results returned. Setting to 0 is equivalent to unlimited
     * @return
     * @throws IOException
     */
    private Collection<DBFeature.IGVFeat> getFeatures(DBObject queryObject, int limit) throws IOException{
        DBCursor cursor = this.collection.find(queryObject);
        cursor.limit(limit >= 0 ? limit : 0);

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

        return features;
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
        MongoFeatureSource source = new MongoFeatureSource(collection, locator.buildLocusIndex);
        FeatureTrack track = new MongoFeatureTrack(collection.getFullName(), collection.getName(), source);
        SearchCommand.registerNamedFeatureSearcher(source);
        newTracks.add(track);
        track.setMargin(0);
        return track;
    }
}
