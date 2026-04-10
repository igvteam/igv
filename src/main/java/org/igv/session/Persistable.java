package org.igv.session;

import org.json.JSONObject;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Interface for a session persistable object.
 *
 * @author User: jrobinso
 * Date: Feb 23, 2010
 */
public interface Persistable {

    /**
     * Marshal object state in XML element
     * @return
     */

    default void marshalXML(Document document, Element element) {
    }

    /**
     * Restore object state from an XML element
     */

    default void unmarshalXML(Element element, Integer version) {
    }

    default void unmarshalJSON(JSONObject jsonObject) {
    }

}
