package com.github.lindenb.jvarkit.util.vcf.bdb;


import htsjdk.variant.variantcontext.Allele;

import com.github.lindenb.jvarkit.util.bdb.XTupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;

public class AlleleBinding extends XTupleBinding<Allele>
	{
    private final static byte REF_A = (byte)1;
    private final static byte ALT_A = (byte)2;
    
    private final static byte REF_C = (byte)3;
    private final static byte ALT_C = (byte)4;
    
    private final static byte REF_G = (byte)5;
    private final static byte ALT_G = (byte)6;
    
    private final static byte REF_T = (byte)7;
    private final static byte ALT_T = (byte)8;
    
    private final static byte REF_N = (byte)9;
    private final static byte ALT_N = (byte)10;
    
    private final static byte NO_CALL = (byte)11;

    private final static byte REF_NOT_SIMPLE = (byte)21;
    private final static byte ALT_NOT_SIMPLE = (byte)22;

	
    public final Allele  entryToAllele(TupleInput in)
    	{
    	return this.entryToObject(in);
    	}
    
	@Override
	public Allele entryToObject(TupleInput in)
		{
		byte b=in.readByte();
		switch(b)
			{
			case ALT_A: 
			case REF_A: return Allele.create("A", (b%2==1));
			case ALT_T: 
			case REF_T: return Allele.create("T", (b%2==1));
			case ALT_G: 
			case REF_G: return Allele.create("G", (b%2==1));
			case ALT_C: 
			case REF_C: return Allele.create("C", (b%2==1));
			case ALT_N: 
			case REF_N: return Allele.create("N", (b%2==1));
			case NO_CALL: return Allele.NO_CALL;
			case ALT_NOT_SIMPLE:
			case REF_NOT_SIMPLE:
				{
				String base=in.readString();
				return Allele.create(base,(b%2==1));
				}
			default: throw new IllegalArgumentException("opcode="+(int)b);
			}		
		}

	@Override
	public void objectToEntry(Allele a, TupleOutput out)
		{
		if(a.isNoCall())
			{
			out.writeByte(NO_CALL);
			return;
			}
		byte seq[]=a.getDisplayBases();
		if(seq.length==1)
			{
			switch(seq[0])
				{
				case 'a':case 'A': out.writeByte(a.isReference()?REF_A:ALT_A); return;
				case 't':case 'T': out.writeByte(a.isReference()?REF_T:ALT_T); return;
				case 'g':case 'G': out.writeByte(a.isReference()?REF_G:ALT_G); return;
				case 'c':case 'C': out.writeByte(a.isReference()?REF_C:ALT_C); return;
				case 'n':case 'N': out.writeByte(a.isReference()?REF_N:ALT_N); return;
				}
			}
		out.writeByte(a.isReference()?REF_NOT_SIMPLE:ALT_NOT_SIMPLE); 
		out.writeString(a.getBaseString());
		}
	}
