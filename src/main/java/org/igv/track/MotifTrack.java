package org.igv.track;

import org.igv.feature.CachingFeatureSource;
import org.igv.feature.Strand;
import org.igv.feature.genome.GenomeManager;
import org.igv.tools.motiffinder.MotifFinderSource;
import org.w3c.dom.Element;

public class MotifTrack extends FeatureTrack {

    private String pattern;
    private Strand strand;

    /**
     * Empty constructor for unmarshalling session
     */
    public MotifTrack() {
    }

    public MotifTrack(String name, String pattern, Strand strand) {
        super(null, name, name);
        this.pattern = pattern;
        this.strand = strand;
        init();
    }

    private void init() {

        MotifFinderSource src = new MotifFinderSource(pattern, strand, GenomeManager.getInstance().getCurrentGenome());
        CachingFeatureSource source = new CachingFeatureSource(src);
        super.init(null, source);
        setVisibilityWindow(1000000);  // Must be done after super.init()
        setSortable(false);
    }

    @Override
    public void setVisibilityWindow(int windowSize) {
        super.setVisibilityWindow(windowSize);
    }

    @Override
    public int getVisibilityWindow() {
        return super.getVisibilityWindow();
    }

    @Override
    public TrackType getType() {
        return TrackType.motif;
    }

    @Override
    public boolean supportsWholeGenome() {
        return false;
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {
        super.unmarshalXML(element, version);
        this.pattern = element.getAttribute("pattern");
        this.strand = Strand.fromString(element.getAttribute("strand"));
        init();
    }

    @Override
    public void marshalJSON(org.json.JSONObject jsonObject) {
        super.marshalJSON(jsonObject);
        jsonObject.put("pattern", pattern);
        jsonObject.put("strand", strand.toString());
    }

    @Override
    public void unmarshalJSON(org.json.JSONObject jsonObject) {
        super.unmarshalJSON(jsonObject);
        this.pattern = jsonObject.getString("pattern");
        this.strand = Strand.fromString(jsonObject.getString("strand"));
        init();
    }
}
