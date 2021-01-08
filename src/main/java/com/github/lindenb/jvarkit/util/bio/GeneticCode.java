/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.bio;

import com.github.lindenb.jvarkit.util.bio.AminoAcids.AminoAcid;

public interface GeneticCode
	{
	static class BasicGeneticCode implements GeneticCode {
		final String ncbiTable;
		BasicGeneticCode(final String ncbiTable) {
			this.ncbiTable = ncbiTable;
			}
		/** convert a base to index */
		private static int base2index(char c)
			{
			switch(c)
				{
				case 'T': case 't': return 0;
				case 'C': case 'c': return 1;
				case 'A': case 'a': return 2;
				case 'G': case 'g': return 3;
				default: return -1;
				}
			}

		private String getNCBITable() {
			return this.ncbiTable;
			}
		/** translate cDNA to aminoacid */
		public char translate(char b1,char b2,char b3)
			{
			final int base1= base2index(b1);
			final int base2= base2index(b2);
			final int base3= base2index(b3);
			if(base1==-1 || base2==-1 || base3==-1)
				{
				return '?';
				}
			else
				{
				return getNCBITable().charAt(base1*16+base2*4+base3);
				}
			}

		}
	
	/** the standard genetic code */
	static final GeneticCode STANDARD=new BasicGeneticCode("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
	
		
	/** mitochondrial genetic code */
	static final GeneticCode MITOCHONDRIAL=new BasicGeneticCode("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG");
		
		
	/** test for stop-codon amino acid*/
	public default boolean isStop(char aminoacid) {
		return aminoacid=='*';
	}
		

	/** test if translation is stop */
	public default boolean isStopCodon(char b1,char b2,char b3)
		{
		return isStop(translate(b1,b2,b3));
		}
	
	/** translate cDNA to aminoacid */
	public char translate(char b1,char b2,char b3);
	
	/** get the standard genetic code */
	public static GeneticCode getStandard()
		{
		return STANDARD;
		}
	
	/** get the mitochondrial genetic code */
	public static GeneticCode getMitochondrial()
		{
		return MITOCHONDRIAL;
		}
	
	/** returns aminoacid to the 3 letter code. Returns *** for stop, return null if to correspondance */
	@Deprecated
	public static String aminoAcidTo3Letters(final char c)
		{
		if(c=='*') return "***";
		AminoAcid aa= AminoAcids.getAminoAcidFromOneLetterCode(c);
		return aa==null?null:aa.getThreeLettersCode();
		}
	}
