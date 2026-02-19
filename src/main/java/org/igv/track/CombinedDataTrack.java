package org.igv.track;

import org.igv.data.CombinedDataSource;
import org.igv.data.DataSource;
import org.json.JSONObject;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

public class CombinedDataTrack extends DataSourceTrack {

    /**
     * Special constructor provided for session unmarshalling
     * @param id
     * @param name
     */
    public CombinedDataTrack(String id, String name) {
        super(null, id, name, null);
    }

    public CombinedDataTrack(DataSource dataSource, String id, String name) {
        super(null, id, name, dataSource);
    }


    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);
        // Un-marshalling handled in IGVSessionReader
    }

    @Override
    public void marshalJSON(JSONObject jsonObject) {
        super.marshalJSON(jsonObject);
        jsonObject.put("track1", ((CombinedDataSource) dataSource).getTrackl().getId());
        jsonObject.put("track2", ((CombinedDataSource) dataSource).getTrack2().getId());
        jsonObject.put("op", ((CombinedDataSource) dataSource).getOperation().toString());
    }

    @Override
    public void unmarshalJSON(JSONObject jsonObject) {
        super.unmarshalJSON(jsonObject);
        // Un-marshalling handled in IGVSessionReader
    }
}

