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

package org.broad.igv.dev.api.batch;

import java.util.List;

/**
 * Interface which must be implemented by
 * any class which is a custom batch script command
 * <p/>
 * @author jacob
 * @date 2012-Nov-30
 * @api
 */

public interface Command {

    /**
     * @param args The string arguments passed AFTER the command name. May be empty. <br/>
     *             Example:<br/>
     *             myCommand abba zabba
     *             <br/>
     *             would provide args = ['abba', 'zabba']
     *             <br/>
     *             myCommand
     *             <br/>
     *             would provide an empty list
     * @return A status string, which should be "OK" if everything went fine, and start with "ERROR" if something went wrong.
     *         Never null
     */
    String run(List<String> args);
}
