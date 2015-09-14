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

package org.broad.igv.cursor;

import javax.swing.*;
import java.awt.*;
import java.io.Serializable;

/**
 * @author jrobinso
 *         Date: 1/26/14
 *         Time: 5:11 PM
 */
public class CursorTrackPanelContainer extends JPanel implements Serializable {

    int height;
    int trackHeight = 60;
    int vmargin = 0;

    @Override
    public void doLayout() {
        height = 0;
        synchronized (this.getTreeLock()) {
            final int width = getWidth();
            for (Component c : this.getComponents()) {
                c.setBounds(0, height, width, trackHeight - 2*vmargin);
                height += trackHeight + vmargin;
            }
        }
    }

    @Override
    public int getHeight() {
        return height;
    }

    @Override
    public Component add(Component comp) {
        int sz = (this.getComponents().length + 1) * trackHeight;
        setPreferredSize(new Dimension(getWidth(), sz));
        return super.add(comp);
    }

    public void setVmargin(int vmargin) {
        this.vmargin = vmargin;
    }

}
