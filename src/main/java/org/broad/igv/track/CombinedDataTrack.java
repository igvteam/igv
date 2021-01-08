package org.broad.igv.track;

import org.broad.igv.data.CombinedDataSource;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

public class CombinedDataTrack extends DataSourceTrack {

    public CombinedDataTrack(String id, String name) {
        super(null, id, name, null);
    }

    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        element.setAttribute("track1", ((CombinedDataSource) dataSource).getTrackl().id);
        element.setAttribute("track2", ((CombinedDataSource) dataSource).getTrack2().id);
        element.setAttribute("op", ((CombinedDataSource) dataSource).getOperation().toString());


//        for (DataTrack track : memberTracks) {
//            Element trackElement = document.createElement(SessionElement.TRACK);
//            track.marshalXML(document, trackElement);
//            element.appendChild(trackElement);
//        }
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        // Un-marshalling handled in IGVSessionReader
    }

}

