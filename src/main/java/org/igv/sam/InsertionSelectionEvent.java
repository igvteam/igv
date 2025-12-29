package org.igv.sam;

import org.igv.event.IGVEvent;

/**
 * Created by jrobinso on 1/12/17.
 */
public record InsertionSelectionEvent(InsertionMarker insertionMarker) implements IGVEvent {

}
