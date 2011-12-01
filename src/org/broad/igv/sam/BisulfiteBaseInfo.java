package org.broad.igv.sam;

import java.awt.Color;

import org.broad.igv.sam.AlignmentTrack.BisulfiteContext;
import org.broad.igv.sam.AlignmentRenderer;
import org.broad.igv.sam.reader.GeraldParser;

/**
 * @author benb
 * Benjamin Berman, University of Southern California
 * 
 * Note that this is only currently supporting a single bisulfite protocol with the following assumptions:
 *  - The first end of paired end reads have C->T variants, while the second ends have G->A
 *  - Bisulfite only affects one strand.  So for the first end, you don't see G->A and for the
 *    second end you don't see C->T.  This allows us to detect true C->T genomic changes by examining
 *    the strand opposite the bisulfite event (we use the DEAMINATION_COLOR for these variants).
 * 
 * This is the Illumina protocol published by Joe Ecker lab (Lister et al. 2009, Lister et al. 2011)
 * and our lab (Berman et al. 2011)
 *
 */
public class BisulfiteBaseInfo {

	public enum DisplayStatus
	{
		NOTHING, COLOR, CHARACTER
	}

	// Constants
	public static Color CYTOSINE_MISMATCH_COLOR = Color.orange; //Color.black;
	public static Color NONCYTOSINE_MISMATCH_COLOR = Color.orange;
	public static Color DEAMINATION_COLOR = new Color(139,94,60);
	public static Color METHYLATED_COLOR = Color.red;
	public static Color UNMETHYLATED_COLOR = Color.blue;

	
	// Private vars
	private DisplayStatus[] displayStatus = null;
	private byte[] displayChars = null;
	private Color[] displayColors = null;
	private AlignmentBlock myBlock = null;
	private BisulfiteContext myContext = null;


	/**
	 * We require loading this once for the entire alignment so that we don't have to calculate reverse complements
	 * more than once.
	 * 
	 * @param inReference
	 * @param inRead
	 * @param alignmentLen
	 * @param block
	 * @param bisulfiteContext
	 */
	public BisulfiteBaseInfo(byte[] inReference, byte[] inRead, int alignmentLen, AlignmentBlock block, BisulfiteContext bisulfiteContext)
	{
		super();

		myBlock = block;
		myContext = bisulfiteContext;
//		System.err.printf("Block=%s, alignment=%s\n", block, block.getBaseAlignment());
		boolean isNegativeStrand = block.getBaseAlignment().isNegativeStrand();
		boolean readIsSecondOfPair = block.getBaseAlignment().isSecondOfPair();

		// We will only need reverse complement if the strand and paired end status don't match (2nd ends are G->A)
		boolean flipRead = isNegativeStrand ^ readIsSecondOfPair;
		byte[] read = (flipRead) ? GeraldParser.reverseComplementCopy(inRead) : inRead;
		byte[] reference = (flipRead) ? GeraldParser.reverseComplementCopy(inReference) : inReference;

		
		displayChars = new byte[alignmentLen];
		displayStatus = new DisplayStatus[alignmentLen];
		displayColors = new Color[alignmentLen];
		
		for (int idxFw = 0; idxFw < alignmentLen; idxFw++)
		{

			//// Everything comes in relative to the positive strand.
			// This is extremely inefficient to do it in here because we have to do it for every position.
			int idx = (flipRead) ? ((alignmentLen - 1) - idxFw) : idxFw;

			// The read base can be an equals sign, so change thet to the actual ref base
			byte refbase = reference[idx];
			byte readbase = read[idx];
			if (readbase == '=') readbase = refbase;

			Color out = null;
			boolean matchesContext = false;

			
			// The logic is organized according to the reference base.  If it's an A/T, it could only be a novel
			// cytosine.  I'm not sure if it's safe to do a switch on refbase here. Would be faster, but I'm
			// being safe (BPB).
			if (AlignmentRenderer.compareBases((byte)'T',refbase) || AlignmentRenderer.compareBases((byte)'A',refbase))
			{
				// Eventually, we could check if it's a novel cytosine.  For now, we show mismatches
				if (!AlignmentRenderer.compareBases(readbase,refbase)) out = NONCYTOSINE_MISMATCH_COLOR;
			}
			else if (AlignmentRenderer.compareBases((byte)'C',refbase))
			{
				// Check if the read is actually a cytosine. If the ref is one and we're not,
				// make a mismatch color
				if (!AlignmentRenderer.compareBases((byte)'C',readbase)  && !AlignmentRenderer.compareBases((byte)'T',readbase)  )
				{
					out = CYTOSINE_MISMATCH_COLOR;
				}
				else	
				{
					BisulfiteContext matchingContext = contextIsMatching(reference, read, idx, bisulfiteContext);
					matchesContext = (matchingContext != null);
					if (matchesContext)
					{
						out = getContextColor(readbase, matchingContext);
					}
				}


			}
			else if (AlignmentRenderer.compareBases((byte)'G',refbase))
			{
				// If it's a guanine in the reference, this could be a special case
				// of deamination of the opposite strand.  I guess we might want to
				// check the context of the opposite strand too.  But for now i'm
				// not.

				if (AlignmentRenderer.compareBases((byte)'A',readbase)) out = DEAMINATION_COLOR;
			}

			
			// Remember, the output should be relative to the FW strand (use idxFw)
			this.displayColors[idxFw] = out;
			if (out == null)
			{
				this.displayStatus[idxFw] = DisplayStatus.NOTHING;
			}
			else
			{
				if (matchesContext)
				{
					// Display the color
					this.displayStatus[idxFw] = DisplayStatus.COLOR;
				}
				else
				{
					// Display the character
					this.displayStatus[idxFw] = DisplayStatus.CHARACTER;
					this.displayChars[idxFw] = 'X';
				}
			}
//			System.err.printf("\tSeting displayStatus[%d] = %s\n", idx, displayStatus[idx]);
		}
		
		this.numDisplayStatus();
	}

	
	protected Color getContextColor(byte readbase,
			BisulfiteContext bisulfiteContext) 
	{
		Color out = null;
		if (AlignmentRenderer.compareBases((byte)'T',readbase))
		{
			out = UNMETHYLATED_COLOR;
		}
		else if (AlignmentRenderer.compareBases((byte)'C',readbase))
		{
			out = METHYLATED_COLOR;
		}
		// C and T should be the only options at this point
		//						else
		//						{
		//							out = Color.yellow;
		//						}
		
		return out;
	}


	/**
	 * @param reference
	 * @param read
	 * @param idx
	 * @param bisulfiteContext
	 * @return Returns the context that matched (in the case of the base class, this is always the same
	 *  as the context passed in.)  If we don't match, return null.
	 */
	protected BisulfiteContext contextIsMatching(byte[] reference, byte[] read, int idx,
			BisulfiteContext bisulfiteContext) {
		
		// Get the context and see if it matches our desired context.
		byte[] preContext = AlignmentTrack.getBisulfiteContextPreContext(bisulfiteContext);
		byte[] postContext = AlignmentTrack.getBisulfiteContextPostContext(bisulfiteContext);

		boolean matchesContext = true;

		// First do the "post" context
		int minLen = Math.min(reference.length, read.length);
		if ((idx+postContext.length)>=minLen)
		{
			matchesContext = false;
		}
		else
		{
			// Cut short whenever we don't match
			for (int posti = 0; matchesContext && (posti < postContext.length); posti++)
			{
				byte contextb = postContext[posti];
				int offsetidx = idx + 1 + posti;
				matchesContext &= positionMatchesContext(contextb, reference, read, offsetidx);

				//				System.err.printf("POST posMatchesContext(posti=%d, contextb=%c, refb=%c, readb=%c, offsetidx=%d) = %s\n",
				//						posti, contextb, reference[offsetidx], read[offsetidx], offsetidx, matchesContext);

			}
		}

		// Now do the pre context
		if ((idx-preContext.length)<0)
		{
			matchesContext = false;
		}
		else
		{
			// Cut short whenever we don't match
			for (int prei = 0; matchesContext && (prei < preContext.length); prei++)
			{
				byte contextb = preContext[prei];
				int offsetidx = idx - (preContext.length - prei);
				matchesContext &= positionMatchesContext(contextb, reference, read, offsetidx);
				//				System.err.printf("PRE posMatchesContext(prei=%d, contextb=%c, refb=%c, readb=%c, offsetidx=%d) = %s\n",
				//						prei, contextb, reference[offsetidx], read[offsetidx], offsetidx, matchesContext);
			}
		}
		
		return (matchesContext) ? bisulfiteContext : null;
	}

	/**
	 * @param contextb The residue in the context string (IUPAC)
	 * @param reference The reference sequence (already checked that offsetidx is within bounds)
	 * @param read The read sequence (already checked that offsetidx is within bounds)
	 * @param offsetidx The index of the position in both reference and read
	 * @return
	 */
	protected boolean positionMatchesContext(byte contextb, byte[] reference, byte[] read, int offsetidx)
	{
		boolean matchesContext = true;
		matchesContext &= AlignmentRenderer.compareBases(contextb, reference[offsetidx]);

		// For the read, we have to handle C separately
		boolean matchesReadContext = false;
		matchesReadContext |= AlignmentRenderer.compareBases(contextb, read[offsetidx]);
		if (AlignmentRenderer.compareBases((byte)'T', read[offsetidx]))
		{
			matchesReadContext |= AlignmentRenderer.compareBases(contextb, (byte)'C');
		}	
		matchesContext &= matchesReadContext;
		
		return matchesContext;
	}		
		
	public Color getDisplayColor(int idx)
	{
		return displayColors[idx];
	}
	
	public DisplayStatus getDisplayStatus(int idx)
	{
		return displayStatus[idx];
	}

	
	public int numDisplayStatus()
	{
		int len = this.displayStatus.length;
		boolean done = false;
		int i =0;
		for (i = 0; (!done && (i<len)); i++)
		{
			done = (this.displayStatus[i] == null);
		}
//		System.err.printf("numDisplayStatus, len=%d\ti=%d\n", len, i);
		return i;
	}

	public double getXaxisShift(int idx) {
		
		double offset = 0.0;
		// This gets CpGs on opposite strands to line up. Only do it for meth values though, not snps

		if (getDisplayStatus(idx).equals(DisplayStatus.COLOR))
		{
			double baseOffset = getBisulfiteSymmetricCytosineShift(myContext);
			offset = offset + ( ( (myBlock.getBaseAlignment().isNegativeStrand() ^ myBlock.getBaseAlignment().isSecondOfPair() ) ? -1 : 1) * baseOffset);
		}
		
		return offset;
	}
	

	/**
     * @param item
     * @return 0 if the cytosine context is not symmetric, else the number of
     * base pairs to shift the position by (caller must be careful to shift
     * relative to the strand the cytosine is on).
     */
    protected double getBisulfiteSymmetricCytosineShift(BisulfiteContext item)
    {
    	double out = 0.0;

    	// The following may be too non-intuitive? BPB
    	switch (item)
    	{
    		case CG:
    		case HCG:
    		case WCG:
    			out = 0.5;
    			break;
    		case CHG:
    			out = 1.0;
    			break;
    		default:
    			out = 0.0;
    			break;
    	}
    	
    	return out;
    }


}
