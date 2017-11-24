/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2017 Broad Institute
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

package org.broad.igv.ui.javafx;

import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.event.GenomeChangeEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.Main.IGVArgs;
import org.broad.igv.ui.panel.FrameManager;

import java.io.IOException;

import static org.broad.igv.prefs.Constants.DEFAULT_FONT_FAMILY;

public class IGVBackendPlaceholder {
    private static Logger log = Logger.getLogger(IGVBackendPlaceholder.class);

    public static void startupInit(IGVArgs igvArgs) {
        final IGVPreferences preferenceManager = PreferencesManager.getPreferences();
        Genome genome = null;
        if (igvArgs.getGenomeId() != null) {
            String genomeId = igvArgs.getGenomeId();
            try {
                genome = GenomeManager.getInstance().loadGenome(genomeId, null);
            } catch (IOException e) {
                log.error("Error loading genome: " + genomeId, e);
            }
        }

        // Ignoring Session for now; don't want to get into that yet.
        if (genome == null) { // && igvArgs.getSessionFile() == null) {
            String genomeId = preferenceManager.getDefaultGenome();
            try {
                genome = GenomeManager.getInstance().loadGenome(genomeId, null);
            } catch (Exception e) {
                log.error("Error loading genome: " + genomeId, e);
            }
        }

        if (genome != null) {
            IGVEventBus.getInstance().post(new GenomeChangeEvent(genome));
            FrameManager.getDefaultFrame().chromosomeNameProperty().set(Globals.CHR_ALL);
        }
    }

    public static final Font getFont(FontWeight fontWeight, double size) {
        final IGVPreferences prefManager = PreferencesManager.getPreferences();
        String fontFamily = prefManager.get(DEFAULT_FONT_FAMILY);

        return Font.font(fontFamily, fontWeight, size);
    }

    public static final Font getFont(double size) {
        final IGVPreferences prefManager = PreferencesManager.getPreferences();
        String fontFamily = prefManager.get(DEFAULT_FONT_FAMILY);

        return Font.font(fontFamily, size);
    }
}
