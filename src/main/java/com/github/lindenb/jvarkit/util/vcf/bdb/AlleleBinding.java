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
package com.github.lindenb.jvarkit.util.vcf.bdb;


import htsjdk.variant.variantcontext.Allele;

import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;

/** Allele Binding for BerkeleyDB */
public class AlleleBinding extends TupleBinding<Allele>
	{
    private final static byte REF_A = (byte)1;
    private final static byte ALT_A = (byte)2;
    private final static byte bytes_A[]=new byte[]{'A'};
    
    private final static byte REF_C = (byte)3;
    private final static byte ALT_C = (byte)4;
    private final static byte bytes_C[]=new byte[]{'C'};
    
    private final static byte REF_G = (byte)5;
    private final static byte ALT_G = (byte)6;
    private final static byte bytes_G[]=new byte[]{'G'};
    
    private final static byte REF_T = (byte)7;
    private final static byte ALT_T = (byte)8;
    private final static byte bytes_T[]=new byte[]{'T'};
    
    private final static byte REF_N = (byte)9;
    private final static byte ALT_N = (byte)10;
    private final static byte bytes_N[]=new byte[]{'N'};
    
    private final static byte NO_CALL = (byte)11;

    private final static byte REF_NOT_SIMPLE = (byte)21;
    private final static byte ALT_NOT_SIMPLE = (byte)22;

	
    public AlleleBinding()
    	{
    	}
    
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
			case REF_A: return Allele.create(bytes_A, (b%2==1));
			case ALT_T: 
			case REF_T: return Allele.create(bytes_T, (b%2==1));
			case ALT_G: 
			case REF_G: return Allele.create(bytes_G, (b%2==1));
			case ALT_C: 
			case REF_C: return Allele.create(bytes_C, (b%2==1));
			case ALT_N: 
			case REF_N: return Allele.create(bytes_N, (b%2==1));
			case NO_CALL: return Allele.NO_CALL;
			case ALT_NOT_SIMPLE:
			case REF_NOT_SIMPLE:
				{
				int n= in.readInt();
				byte bases[]=new byte[n];
				in.read(bases);
				return Allele.create(bases,  (b%2==1));
				}
			default: throw new IllegalArgumentException("opcode="+(int)b);
			}		
		}

	@Override
	public void objectToEntry(Allele a, TupleOutput out)
		{
		if(a==null || a.isNoCall())
			{
			out.writeByte(NO_CALL);
			return;
			}
		byte seq[]=a.getDisplayBases();
		if(seq.length==0) throw new RuntimeException("n=0 for "+a);
		
		if(seq.length==1)
			{
			switch(seq[0])
				{
				case 'a':case 'A': out.writeByte(a.isReference()?REF_A:ALT_A); return;
				case 't':case 'T': out.writeByte(a.isReference()?REF_T:ALT_T); return;
				case 'g':case 'G': out.writeByte(a.isReference()?REF_G:ALT_G); return;
				case 'c':case 'C': out.writeByte(a.isReference()?REF_C:ALT_C); return;
				case 'n':case 'N': out.writeByte(a.isReference()?REF_N:ALT_N); return;
				default:break;
				}
			}
		
		out.writeByte(a.isReference()?REF_NOT_SIMPLE:ALT_NOT_SIMPLE);
		out.writeInt(seq.length);
		out.write(seq);
		}
	}
