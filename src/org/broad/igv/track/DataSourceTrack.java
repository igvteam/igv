/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


package org.broad.igv.track;


import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.data.CombinedDataSource;
import org.broad.igv.data.CoverageDataSource;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlSeeAlso;
import javax.xml.bind.annotation.XmlType;
import javax.xml.namespace.QName;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;


/**
 * @author jrobinso
 */
@XmlType(factoryMethod = "getNextTrack")
@XmlSeeAlso({CombinedDataSource.class})
public class DataSourceTrack extends DataTrack {

    private static Logger log = Logger.getLogger(DataSourceTrack.class);

    //All tracks have label "Track", we need to specify the type sometimes
    //but still preserve backwards compatibility
    @XmlAttribute
    protected Class clazz = DataSourceTrack.class;

    private DataSource dataSource;

    /**
     * We often don't have data when the track is constructed,
     * need to rescale when data is first loaded
     */
    private boolean firstDataLoaded = false;

    private boolean rescaleOnFirst = false;

    private static final String COMBINED_DATA_SOURCE = "COMBINED_DATA_SOURCE";

    @SubtlyImportant
    private DataSourceTrack(){
        super(null, null, null);
    }

    public DataSourceTrack(ResourceLocator locator, String id, String name, DataSource dataSource) {
        super(locator, id, name);

        this.dataSource = dataSource;

        if(this.dataSource != null){
            setTrackType(dataSource.getTrackType());
            List<LocusScore> scores = this.dataSource.getSummaryScoresForRange(Globals.CHR_ALL, -1, -1, 0);

            if(scores.size() == 0) rescaleOnFirst = true;

            initScale(dataSource, scores);
        }
    }

    void initScale(DataSource dataSource, List<LocusScore> scores){

        float min = (float) dataSource.getDataMin();
        float max = (float) dataSource.getDataMax();
        float baseline = 0;

        // If the range is all + numbers set the min to zero
        if (min > 0) {
            min = 0;
        }
        for (LocusScore score : scores) {
            max = Math.max(max, score.getScore());
        }

        setDataRange(new DataRange(min, baseline, max));
    }

    public List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        List<LocusScore> tmp = dataSource.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
        if (tmp == null) tmp = Collections.EMPTY_LIST;
        if(!firstDataLoaded && rescaleOnFirst){
            initScale(dataSource, tmp);
            firstDataLoaded = true;
        }
        return tmp;
    }


    @Override
    public void setWindowFunction(WindowFunction statType) {
        clearCaches();
        if (dataSource != null) {
            dataSource.setWindowFunction(statType);
        }
    }

    public boolean isLogNormalized() {
        return dataSource.isLogNormalized();
    }


    public WindowFunction getWindowFunction() {
        //For JAXB session loading/unloading, dataSource might be null, and we need to guard against that
        return dataSource != null ? dataSource.getWindowFunction() : null;
    }


    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return dataSource != null ? dataSource.getAvailableWindowFunctions() : null;
    }

    @Override
    public void restorePersistentState(Node node, int version) throws JAXBException {
        super.restorePersistentState(node, version);
        if (node.hasChildNodes()) {
            NodeList childNodes = node.getChildNodes();
            for (int ii = 0; ii < childNodes.getLength(); ii++) {
                Node child = childNodes.item(ii);
                String nodeName = child.getNodeName();
                if (nodeName.contains("#text")) continue;

                if (nodeName.equalsIgnoreCase(COMBINED_DATA_SOURCE)) {
                    dataSource = IGVSessionReader.getJAXBContext().createUnmarshaller().unmarshal(child, CombinedDataSource.class).getValue();
                }
            }
        }
    }

    public void marshalSource(Marshaller m, Element trackElement) throws JAXBException {
        if (dataSource == null) return;
        DataSource rawSource = dataSource;

        if(rawSource instanceof CombinedDataSource){
            JAXBElement element = new JAXBElement<CombinedDataSource>(new QName("", COMBINED_DATA_SOURCE), CombinedDataSource.class,
                    (CombinedDataSource) rawSource);
            m.marshal(element, trackElement);
        }
    }

    public void updateTrackReferences(List<Track> allTracks) {
        if (dataSource instanceof CombinedDataSource) {
            ((CombinedDataSource) dataSource).updateTrackReferences(allTracks);
        }
    }

    @Override
    public void dispose() {
        super.dispose();
        if(dataSource != null) {
            dataSource.dispose();
        }
    }

    @SubtlyImportant
    @XmlAttribute
    private void setNormalize(boolean normalize){
        if (dataSource != null && dataSource instanceof CoverageDataSource) {
            ((CoverageDataSource) dataSource).setNormalize(normalize);
        }
    }

    @SubtlyImportant
    private boolean getNormalize(){
        if (dataSource != null && dataSource instanceof CoverageDataSource) {
            return ((CoverageDataSource) dataSource).getNormalize();
        }
        return false;
    }

    @SubtlyImportant
    private static DataSourceTrack getNextTrack(){
        DataSourceTrack out = (DataSourceTrack) IGVSessionReader.getNextTrack();
        if (out == null){
            out = new DataSourceTrack();
        }
        return out;
    }
}
