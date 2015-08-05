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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.util;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;

/**
 * @author eflakes
 */
abstract public class Filter {

    private String name;
    private LinkedHashSet<FilterElement> elements =
            new LinkedHashSet<FilterElement>();
    private boolean isEnabled = true;

    public Filter(String name, List<? extends FilterElement> elements) {
        this.name = name;

        if (elements != null) {
            this.elements.addAll(elements);
        }
    }

    public void removeAll() {
        elements.clear();
    }

    public void setEnabled(boolean value) {
        isEnabled = value;
    }

    public boolean isEnabled() {
        return isEnabled;
    }

    public boolean isEmpty() {
        return elements.isEmpty();
    }

    public Iterator getFilterElements() {
        return elements.iterator();
    }

    public void add(FilterElement element) {
        elements.add(element);
    }

    public void remove(FilterElement element) {
        elements.remove(element);
    }

    /**
     * Evaluate the FilterElement set.
     *
     * @return
     */
    abstract protected void evaluate();

    public String getName() {
        return name;
    }
}
