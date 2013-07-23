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
	/** translate cDNA to aminoacid */
	public char translate(char b1,char b2,char b3)
		{
		int base1= base2index(b1);
		int base2= base2index(b2);
		int base3= base2index(b3);
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
	
	}
