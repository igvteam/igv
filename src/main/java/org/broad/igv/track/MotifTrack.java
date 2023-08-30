package org.broad.igv.track;

import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.tools.motiffinder.MotifFinderSource;
import org.w3c.dom.Document;
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
        super.init(null, new MotifFinderSource(pattern, strand, GenomeManager.getInstance().getCurrentGenome()));
        setSortable(false);
    }

    @Override
    public void marshalXML(Document document, Element element) {
        element.setAttribute("pattern", pattern);
        element.setAttribute("strand", String.valueOf(strand));
        super.marshalXML(document, element);
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {
        super.unmarshalXML(element, version);
        this.pattern = element.getAttribute("pattern");
        this.strand = Strand.valueOf(element.getAttribute("strand"));
        init();
    }
}
