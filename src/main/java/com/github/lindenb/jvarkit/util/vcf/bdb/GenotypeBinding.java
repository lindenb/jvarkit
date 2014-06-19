package com.github.lindenb.jvarkit.util.vcf.bdb;


import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;

import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;

public class GenotypeBinding extends AbstractVCFBinding<Genotype>
	{	
	private final AlleleBinding alleleBinding=new AlleleBinding();
	
	@Override
	public Genotype entryToObject(TupleInput in)
		{
		GenotypeBuilder gb=new GenotypeBuilder(in.readString());
		if(in.readBoolean()) gb.DP(in.readInt());
		if(in.readBoolean()) gb.AD(arrayOfIntToEntry(in));
		if(in.readBoolean()) gb.GQ(in.readInt());
		if(in.readBoolean()) gb.PL(arrayOfIntToEntry(in));
		
		/* ALLELES ======================================== */
		int n=in.readInt();
		List<Allele> alleles=new ArrayList<Allele>(n);
		for(int i=0;i< n;++i)
			{
			alleles.add(this.alleleBinding.entryToObject(in));
			}
		gb.alleles(alleles);
		/* ATTRIBUTES ===================================== */
		n=in.readInt();
		for(int i=0;i< n;++i)
			{
			String key=in.readString();
			gb.attribute(key, super.readAttribute(in));
			}
		 
		return gb.make();
		}

	@Override
	public void objectToEntry(Genotype g, TupleOutput out)
		{
		out.writeString(g.getSampleName());
		/* ============================================= */
		if(g.hasDP())
			{	
			out.writeBoolean(true);
			out.writeInt(g.getDP());
			}
		else
			{
			out.writeBoolean(false);
			}
		/* ============================================= */
		if(g.hasAD())
			{	
			out.writeBoolean(true);
			arrayOfIntToEntry(g.getAD(),out);
			}
		else
			{
			out.writeBoolean(false);
			}
		/* ============================================= */
		if(g.hasGQ())
			{	
			out.writeBoolean(true);
			out.writeInt(g.getGQ());
			}
		else
			{
			out.writeBoolean(false);
			}
		/* ============================================= */
		if(g.hasPL())
			{	
			out.writeBoolean(true);
			arrayOfIntToEntry(g.getPL(),out);
			}
		else
			{
			out.writeBoolean(false);
			}
		/* ALLELES ======================================== */
		List<Allele> alleles=g.getAlleles();
		out.writeInt(alleles.size());
		for(Allele a:alleles)
			{
			this.alleleBinding.objectToEntry(a, out);
			}
		
		
		/* ATTRIBUTES ===================================== */
		Map<String,Object> xAtts=g.getExtendedAttributes();
		out.writeInt(xAtts.size());
		for(String attKey:xAtts.keySet())
			{
			out.writeString(attKey);
			super.writeAttribute(out, xAtts.get(attKey));
			}
		}
	
	private void arrayOfIntToEntry(int array[],TupleOutput out)
		{
		out.writeInt(array.length);
		for(int v:array) out.writeInt(v);
		}
	
	private int[] arrayOfIntToEntry(TupleInput in)
		{
		int array[]=new int[in.readInt()];
		for(int i=0;i< array.length;++i)
			{
			array[i]=in.readInt();
			}
		return array;
		}

	}
