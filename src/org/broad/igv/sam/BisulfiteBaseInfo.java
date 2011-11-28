package org.broad.igv.sam;

import java.awt.Color;

import org.broad.igv.sam.AlignmentTrack.BisulfiteContext;
import org.broad.igv.sam.AlignmentRenderer;
import org.broad.igv.sam.reader.GeraldParser;

public class BisulfiteBaseInfo {

	public enum DisplayStatus
	{
		NOTHING, COLOR, CHARACTER
	}

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

			if (AlignmentRenderer.compareBases((byte)'T',refbase) || AlignmentRenderer.compareBases((byte)'A',refbase))
			{
				// Eventually, we could check if it's a novel cytosine
			}
			else if (AlignmentRenderer.compareBases((byte)'C',refbase))
			{
				// Check if we're actually a cytosine. If the ref is one and we're not,
				// make a mismatch color
				if ((readbase != 'C')  && (readbase != 'T'))
				{
					//out = Color.yellow;
				}
				else	
				{
					// Get the context and see if it matches our desired context.
					boolean matchesContext = true;

					// Change this to do it in an automatic way.  Requires the option constants to be simple
					// strings, where you unambiguously know where the central "C" is.
					// +1 base
					if ( ((idx+1) >= reference.length) && ((idx+1) >= read.length ))
					{
						matchesContext = false;
					}
					else
					{
						switch (bisulfiteContext)
						{
						case CG:
						case HCG:
						case WCG:
							if (!AlignmentRenderer.compareBases((byte)'G',reference[idx+1])) matchesContext = false;
							if (!AlignmentRenderer.compareBases((byte)'G',read[idx+1])) matchesContext = false;  // Be careful, must consider bisulfite C->Ts here if applicable
							break;
						case GCH:
							if (!AlignmentRenderer.compareBases((byte)'H',reference[idx+1])) matchesContext = false;
							if (!AlignmentRenderer.compareBases((byte)'H',read[idx+1])) matchesContext = false;  // Be careful, must consider bisulfite C->Ts here if applicable
							break;
						}
					}
					// -1 base
					if ((idx-1) < 0)
					{
						matchesContext = false;
					}
					else
					{
						switch (bisulfiteContext)
						{
						case WCG:
							if (!AlignmentRenderer.compareBases((byte)'W',reference[idx-1])) matchesContext = false;
							if (!AlignmentRenderer.compareBases((byte)'W',read[idx-1])) matchesContext = false;
							break;
						case HCG:
							if (!AlignmentRenderer.compareBases((byte)'H',reference[idx-1])) matchesContext = false;
							if (!AlignmentRenderer.compareBases((byte)'H',read[idx-1])) matchesContext = false;
							break;
						case GCH:
							if (!AlignmentRenderer.compareBases((byte)'G',reference[idx-1])) matchesContext = false;
							if (!AlignmentRenderer.compareBases((byte)'G',read[idx-1])) matchesContext = false;
							break;
						}
					}



					if (matchesContext)
					{
						if (AlignmentRenderer.compareBases((byte)'T',readbase))
						{
							out = Color.BLUE;
						}
						else if (AlignmentRenderer.compareBases((byte)'C',readbase))
						{
							out = Color.red;
						}
						else
						{
							out = Color.yellow;
						}
					}
				}


			}
			else if (AlignmentRenderer.compareBases((byte)'G',refbase))
			{
				// If it's a guanine in the reference, check if it indicates a variant.  Only
				// print a color if it indicates a variant

				if (AlignmentRenderer.compareBases((byte)'A',readbase)) out = new Color(139,94,60);
			}

			
			// Remember, the output should be relative to the FW strand
			this.displayColors[idxFw] = out;
			if (out == null)
			{
				this.displayStatus[idxFw] = DisplayStatus.NOTHING;
			}
			else
			{
				boolean isMethColor = ((out.equals(Color.blue)) || (out.equals(Color.red)));
				this.displayStatus[idxFw] = (isMethColor) ? DisplayStatus.COLOR : DisplayStatus.CHARACTER;
			}
//			System.err.printf("\tSeting displayStatus[%d] = %s\n", idx, displayStatus[idx]);
		}
		
		this.numDisplayStatus();
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

		if (myContext.equals(AlignmentTrack.BisulfiteContext.CG) && (getDisplayStatus(idx).equals(DisplayStatus.COLOR)) )
		{
			// Would be better if they could meet in the middle, more accurate.  But can't do it with ints.
			offset = offset + ( ( (myBlock.getBaseAlignment().isNegativeStrand() ^ myBlock.getBaseAlignment().isSecondOfPair() ) ? -1 : 1) * 0.5);
			//if (isNegativeStrand) locToUse--;
		}

		return offset;
	}
	
	
}
