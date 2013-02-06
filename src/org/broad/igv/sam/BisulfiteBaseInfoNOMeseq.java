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
