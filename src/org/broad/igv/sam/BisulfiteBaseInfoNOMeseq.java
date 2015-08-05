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

/**
 * 
 */
package org.broad.igv.sam;

import org.broad.igv.sam.AlignmentTrack.BisulfiteContext;

import java.awt.*;

/**
 * @author benb
 *
 */
public class BisulfiteBaseInfoNOMeseq extends BisulfiteBaseInfo {

	
	public static Color GC_METH_COLOR = new Color(0,204,153);
	public static Color GC_UNMETH_COLOR = new Color(255,127,191);
	public static Color CG_METH_COLOR = Color.black;
    public static Color CG_UNMETH_COLOR = Color.white;
    public static final BisulfiteContext[] BISULFITE_CONTEXTS = new BisulfiteContext[]{BisulfiteContext.HCG, BisulfiteContext.GCH};


    /**
	 * @param inReference
     * @param baseAlignment
	 * @param block
	 * @param bisulfiteContext
	 */
	public BisulfiteBaseInfoNOMeseq(byte[] inReference, Alignment baseAlignment, AlignmentBlock block,
			BisulfiteContext bisulfiteContext) {
		super(inReference, baseAlignment, block, bisulfiteContext);
	}


	@Override
	protected Color getContextColor(byte readbase,
			BisulfiteContext bisulfiteContext) {
		
		Color out = null;
		if (bisulfiteContext.equals(BisulfiteContext.HCG))
		{
			out = (AlignmentUtils.compareBases((byte) 'T', readbase)) ? CG_UNMETH_COLOR : CG_METH_COLOR;
		}
		else
		{
			out = (AlignmentUtils.compareBases((byte) 'T', readbase)) ? GC_UNMETH_COLOR : GC_METH_COLOR;
		}
		return out;
	}


	@Override
	protected BisulfiteContext contextIsMatching(byte[] reference, byte[] read, int idx,
			BisulfiteContext bisulfiteContext) {

		for (BisulfiteContext context : BISULFITE_CONTEXTS)
		{
			if (super.contextIsMatching(reference, read, idx,context) != null) return context;
		}
		return null;
	}


	@Override
	protected double getBisulfiteSymmetricCytosineShift(BisulfiteContext item) {
		// TODO Auto-generated method stub
		return 0.0;
	}
	
	
	

}
