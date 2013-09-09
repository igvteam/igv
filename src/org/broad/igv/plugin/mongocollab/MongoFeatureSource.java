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
import org.broad.igv.feature.IGVFeature;
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
public class MongoFeatureSource implements FeatureSource {

    private int featureWindowSize = 1000000;

    private DBCollection collection;

    public MongoFeatureSource(DBCollection collection) {
        this.collection = collection;
    }

    private DBObject createQueryObject(String chr, int start, int end){
        BasicDBObject query = new BasicDBObject("Chr", chr);
        //TODO figure out how to do OR so we can query range properly
        //query.append("Start", (new BasicDBObject("$gt", start)).append("$lt", end));

        return query;
    }

    @Override
    public Iterator<IGVFeature> getFeatures(String chr, int start, int end) throws IOException {
        this.collection.setObjectClass(MongoCollabPlugin.FeatDBObject.class);
        DBCursor cursor = this.collection.find(createQueryObject(chr, start, end));

        List<IGVFeature> features = new ArrayList<IGVFeature>();


        while (cursor.hasNext()) {
            DBObject obj = cursor.next();
            MongoCollabPlugin.FeatDBObject feat = (MongoCollabPlugin.FeatDBObject) obj;
            features.add(feat.createBasicFeature());
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

    public static FeatureTrack loadFeatureTrack(DBCollection collection, List<Track> newTracks) {
        MongoFeatureSource source = new MongoFeatureSource(collection);
        FeatureTrack track = new FeatureTrack(collection.getName(), collection.getFullName(), source);
        newTracks.add(track);
        track.setMargin(0);
        return track;
    }

}
