/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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


public abstract class GeneticCode
	{
	/** the standard genetic code */
	private static final GeneticCode STANDARD=new GeneticCode()
		{
		@Override
		protected String getNCBITable()
			{
			return "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
			}
		};
	/** mitochondrial genetic code */
	private static final GeneticCode MITOCHONDRIAL=new GeneticCode()
		{
		@Override
		protected String getNCBITable()
			{
			return "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
			}
		};
		
	/** test for stop-codon amino acid*/
	public boolean isStop(char aminoacid) {
		return aminoacid=='*';
	}
		
	/** get the genetic-code table (NCBI data) */ 
	protected abstract String getNCBITable();
	
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
	
	/** test if translation is stop */
	public boolean isStopCodon(char b1,char b2,char b3)
		{
		return isStop(translate(b1,b2,b3));
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
	public static String aminoAcidTo3Letters(final char c)
		{
		switch(Character.toUpperCase(c))
			{
			case 'A': return "Ala";
			case 'R': return "Arg";
			case 'N': return "Asn";
			case 'D': return "Asp";
			case 'C': return "Cys";
			case 'E': return "Glu";
			case 'Q': return "Gln";
			case 'G': return "Gly";
			case 'H': return "His";
			case 'I': return "Ile";
			case 'L': return "Leu";
			case 'K': return "Lys";
			case 'M': return "Met";
			case 'F': return "Phe";
			case 'P': return "Pro";
			case 'S': return "Ser";
			case 'T': return "Thr";
			case 'W': return "Trp";
			case 'Y': return "Tyr";
			case 'V': return "Val";
			case '*': return "***";
			default: return null;
			}
		}
	}
