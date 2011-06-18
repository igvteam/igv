/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package edu.mit.broad.prodinfo.sequence;

import java.util.Iterator;

public class WindowSlider implements Iterator<SequenceRegion> {
	
	private int windowSize;
	private int overlap;
	private int atPosition = 0; 
	private int seqStart;
	private Sequence seq;

	public static WindowSlider getSlider(Sequence seq, int windowSize, int overlap) {
		WindowSlider ws = new WindowSlider(windowSize, overlap);
		ws.seq = seq;
		ws.seqStart = 0;
		return ws;
	}
	
	public static WindowSlider getSlider(SequenceRegion seq, int windowSize, int overlap) {
		//System.out.println("Overloaded method used");
		WindowSlider ws = new WindowSlider(windowSize, overlap);
		ws.seq = seq;
		ws.atPosition = seq.getStart();
		ws.seqStart = seq.getStart();
		return ws;
	}
	
	WindowSlider(int windowSize, int overlap) {
		this.windowSize = windowSize;
		this.overlap    = overlap;
	}
	
	public boolean hasNext() {
		//System.out.println("atPosition " + atPosition  + " seq_end " + seq.getEnd());
		return atPosition < seq.getEnd() - windowSize;
	}

	public SequenceRegion next() {
		SequenceRegion window = seq.getRegion(atPosition - seqStart, atPosition - seqStart + windowSize);
		window.setStart(atPosition);
		window.setEnd(atPosition + windowSize);
		atPosition = atPosition + windowSize - overlap;
		return window;
	}

	public void remove() {
		// DO NOTHING

	}

}
