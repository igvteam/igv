package org.broad.igv.feature.bedpe;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.*;
import java.util.List;


/**
 * Created by jrobinso on 6/29/18.
 */
public class InteractionTrack extends AbstractTrack {

    private Map<String, List<BedPEFeature>> featureMap;

    private PEArcRenderer renderer;

    public InteractionTrack(ResourceLocator locator, List<BedPEFeature> featureList, Genome genome) {
        super(locator);
        init(featureList, genome);
        renderer = new PEArcRenderer();
        setHeight(250, true);
        setColor(new Color(180, 25, 137));
    }

    private void init(List<BedPEFeature> featureList, Genome genome) {

        this.featureMap = new HashMap<>();

        // Sort feature lists by "start" (minimum of start1, start2)
       Collections.sort(featureList, (o1, o2) -> o1.getStart() - o2.getStart());

        for(BedPEFeature f : featureList) {

            String key;
            if(f.chr1.equals(f.chr2)) {
                key = genome == null ? f.chr1 : genome.getCanonicalChrName(f.chr1);
            }
            else {
                key = "OTHER";
            }

            List<BedPEFeature> features = featureMap.get(key);
            if(features == null) {
                features = new ArrayList<>();
                featureMap.put(key, features);
            }

            features.add(f);
        }



    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return true;
    }

    @Override
    public void load(ReferenceFrame frame) {
        // Nothing to do, this track is pre-loaded
    }

    @Override
    public void render(RenderContext context, Rectangle rect) {

        String chr = context.getReferenceFrame().getChrName();
        List<BedPEFeature> features = featureMap.get(chr);

        if(features != null) {
            renderer.render(features, context, rect, this);
        }


    }
}
