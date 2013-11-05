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

package org.broad.igv.track;


import java.awt.*;
import java.util.Collection;

/**
 * Track to serve as a container for several tracks rendered on top of each other
 *
 * @author jacob
 * @date 2013-Nov-05
 */
public class MergedTracks extends AbstractTrack{

    private Collection<Track> trackList;

    public MergedTracks(String id, String name, Collection<Track> trackList){
        super(id, name);
        this.trackList = trackList;
    }

    @Override
    public void render(RenderContext context, Rectangle rect) {
        resetLastY();
        for(Track track: trackList){
            if(isRepeatY(rect)){
                track.overlay(context, rect);
            }else{
                track.render(context, rect);
                lastRenderY = rect.y;
            }

        }
    }

    @Override
    public int getHeight() {
        int height = super.getHeight();
        for(Track track: trackList){
            height = Math.max(height, track.getHeight());
        }
        return height;
    }


}
