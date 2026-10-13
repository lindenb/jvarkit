/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.illumina;

import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

/** constants for illumina */
public class IlluminaConstants {
	private static final Logger LOG = Logger.of(IlluminaConstants.class);
	/* https://help.dragen.illumina.com/dragen-v4.5-trupath/product-guides/dragen-v4.5/dragen-trupath-pipeline */
	public static enum TruePathPhase { Phase_1, Phase_2, Phase_undefined };
	
	/** Phasing probability in Phred-scale  */
	public static final String TRUEPATH_PP = "pp";
	/** Haplotype tag for all reads where abs(pp)>=10 often 1 or 2: HP:i:1 HP:i:2 */
	public static final String TRUEPATH_HP = "HP";
	/** Phase block tag */
	public static final String TRUEPATH_PS = "PS";
	
	/** return true if any program record in the header looks like a dragen command line with 'proximity ensabled' */
	public static boolean isDragenProximityEnabled(final SAMFileHeader header) {
		return header.getProgramRecords()
				.stream()
				.map(PG->PG.getCommandLine())
				.anyMatch(CL->CL!=null && CL.contains("dragen") && CL.contains("--enable-proximity true"));
		
		}
	
	/** extract truepath Phase from SamRecord */
	public static TruePathPhase getTruePathPhase(final SAMRecord rec) {
		 final Object hp_tag = rec.getAttribute(IlluminaConstants.TRUEPATH_HP);
		if(hp_tag==null ) {
			return TruePathPhase.Phase_undefined;
			}
		if( !(hp_tag instanceof Number)) {
			LOG.warn("unexpected phase (NaN) in "+rec);
			return TruePathPhase.Phase_undefined;
			}
		
		final int hp_phase= Number.class.cast(hp_tag).intValue();
		switch(hp_phase) {
			case 1: return TruePathPhase.Phase_1;
			case 2: return TruePathPhase.Phase_2;
			default: LOG.warn("unexpected phase in "+rec); return TruePathPhase.Phase_undefined;
			}
		}
}
