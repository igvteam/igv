/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.cli_plugin;

import org.apache.log4j.Logger;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.LineReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * User: jacob
 * Date: 2012-Sep-27
 */
public class AsciiDecoder<D extends Feature> implements FeatureDecoder<D> {

    private static Logger log = Logger.getLogger(AsciiDecoder.class);

    protected LineFeatureDecoder<D> lineFeatureDecoder;

    public AsciiDecoder() {
    }

    public AsciiDecoder(LineFeatureDecoder<D> lineFeatureDecoder) {
        this.lineFeatureDecoder = lineFeatureDecoder;
    }

    public Iterator<D> decodeAll(InputStream is, boolean strictParsing) throws IOException {

        List<D> featuresList = new ArrayList<D>();
        String line;
        D feat;


        BufferedReader bis = new BufferedReader(new InputStreamReader(is));
        QueuingLineReader lrw = new QueuingLineReader(bis);

        lrw.setQueueing(true);
        lineFeatureDecoder.readHeader(lrw);
        lrw.setQueueing(false);

        while ((line = lrw.readLine()) != null) {
            try {
                feat = decode(line);
                if (feat != null) {
                    featuresList.add(feat);
                }
            } catch (Exception e) {
                log.error(e.getMessage(), e);
                if (strictParsing) {
                    throw new RuntimeException(e);
                }
            }
        }

        is.close();
        return featuresList.iterator();
    }

    public D decode(String line) {
        return this.lineFeatureDecoder.decode(line);
    }

    @Override
    public void setAttributes(List<Map<String, Object>> attributes) {
    }

    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
    }

    /**
     * Wrap an AsciiFeatureCodec into implementing LineFeatureDecoder
     *
     * @param <T>
     */
    public static class DecoderWrapper<T extends Feature> extends AsciiDecoder<T> implements LineFeatureDecoder<T> {

        private AsciiFeatureCodec<T> wrappedCodec;

        public DecoderWrapper(AsciiFeatureCodec<T> wrappedCodec) {
            this.wrappedCodec = wrappedCodec;
            this.lineFeatureDecoder = this;
        }

        @Override
        public T decode(String line) {
            return wrappedCodec.decode(line);
        }

        @Override
        public Object readHeader(LineReader reader) throws IOException{
            return wrappedCodec.readHeader(reader);
        }

    }

    /**
     * The purpose of this class is 2-fold:
     * 0. Wrap a BufferedReader in a LineReader so the interface works
     * 1. Keep track of read lines so reading the header doesn't throw away potential data
     * When the QueuingLineReader is set to {@code queueing}, it queues each readLine.
     * Once {@code queueing} is turned off the queued lines are returned. The queue is
     * cleared when {@code queueing} is toggled on the theory that only blocks
     * will want to be read.
     * @author jacob
     * @since 3 Jul 2013
     */
    private static class QueuingLineReader implements org.broad.tribble.readers.LineReader{

        private BufferedReader wrappedReader;

        private Queue<String> lineBuffer = new ArrayDeque<String>();
        private boolean queueing = false;

        QueuingLineReader(BufferedReader lineReader){
            this.wrappedReader = lineReader;
        }

        void setQueueing(boolean queueing){
            this.queueing = queueing;
            if(this.queueing) lineBuffer.clear();
        }

        @Override
        public String readLine() throws IOException {
            String line;
            if(queueing){
                line = this.wrappedReader.readLine();
                if(line != null) lineBuffer.add(line);
            }else{
                line = lineBuffer.poll();
                if(line == null) line = this.wrappedReader.readLine();
            }
            return line;
        }

        @Override
        public void close() {
            this.lineBuffer = null;
            try {
                this.wrappedReader.close();
            } catch (IOException e) {
                log.error(e.getMessage(), e);
            }
        }
    }
}
