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
package org.igv.ui.toolbar;

import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import org.apache.log4j.Logger;
import org.broad.igv.ui.panel.FrameManager;

/**
 * @author jrobinso on 7/6/17.
 * @author eby JavaFX port
 */
// TODO: need to implement JavaFX version of SearchHints.
public class SearchTextField extends TextField {

    private static Logger log = Logger.getLogger(SearchTextField.class);

    public SearchTextField() {
        setTooltip(new Tooltip("Enter a gene or locus, e.f. EGFR,   chr1,   or chr1:100,000-200,000"));
        setOnAction(actionevent -> searchByLocus(getText()));
    }


    public void searchByLocus(final String searchText) {
        (new SearchCommand(FrameManager.getDefaultFrame(), searchText)).execute();
    }
}
