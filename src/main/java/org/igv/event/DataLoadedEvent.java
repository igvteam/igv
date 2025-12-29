package org.igv.event;

import org.igv.ui.panel.ReferenceFrame;

/**
 * User: jacob
 * Date: 2013-Feb-06
 */
public record DataLoadedEvent(ReferenceFrame referenceFrame) implements IGVEvent {}
