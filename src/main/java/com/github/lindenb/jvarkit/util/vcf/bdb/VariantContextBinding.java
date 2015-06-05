package com.github.lindenb.jvarkit.util.vcf.bdb;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;

public class VariantContextBinding extends AbstractVCFBinding<VariantContext>
	{
	private AlleleBinding alleleBinding=new AlleleBinding();
	private GenotypeBinding genotypeBinding=new GenotypeBinding();
	
	@Override
	public VariantContext entryToObject(TupleInput in)
		{
		VariantContextBuilder vcb=new VariantContextBuilder();
		vcb.chr(in.readString());
		vcb.start(in.readInt());
		vcb.stop(in.readInt());
		if(in.readBoolean())
			{
			vcb.id(in.readString());
			}
		/* ALLELES ======================== */
		int n=in.readInt();
		List<Allele> alleles=new ArrayList<Allele>(n);
		for(int i=0;i< n;++i)
			{
			alleles.add(this.alleleBinding.entryToObject(in));
			}
		vcb.alleles(alleles);
		/* QUAL ======================== */
		if(in.readBoolean())
			{
			vcb.log10PError(in.readDouble());
			}
		/* FILTERS ======================== */
		int n_filters=in.readInt();
		Set<String> filters=new HashSet<String>(n_filters);
		for(int i=0;i< n_filters;++i)
			{
			filters.add(in.readString());
			}
		vcb.filters(filters);
		/* INFO ======================== */
		int n_infokeys=in.readInt();
		Map<String,Object> hash=new HashMap<String,Object>(n_infokeys);
		for(int i=0;i< n_infokeys;++i)
			{
			String key=in.readString();
			Object value=super.readAttribute(in);
			hash.put(key, value);
			}
		vcb.attributes(hash);
		
		/* GENOTYPES ======================== */
		
		int n_genotypes=in.readInt();
		List<Genotype> genotypes=new ArrayList<Genotype>(n_genotypes);
		for(int i=0;i< n_genotypes;++i)
			{
			genotypes.add(this.genotypeBinding.entryToObject(in));
			}
		vcb.genotypes(genotypes);
		
		return vcb.make();
		}

	@Override
	public void objectToEntry(VariantContext ctx, TupleOutput out)
		{
		out.writeString(ctx.getContig());
		out.writeInt(ctx.getStart());
		out.writeInt(ctx.getEnd());
		if(ctx.hasID())
			{
			out.writeBoolean(true);
			out.writeString(ctx.getID());
			}
		else
			{
			out.writeBoolean(false);
			}
		
		List<Allele> alleles=ctx.getAlleles();
		out.writeInt(alleles.size());
		for(Allele allele:alleles)
			{
			this.alleleBinding.objectToEntry(allele, out);
			}
		/* QUAL ======================== */
		if(ctx.hasLog10PError())
			{
			out.writeBoolean(true);
			out.writeDouble(ctx.getLog10PError());
			}
		else
			{
			out.writeBoolean(false);
			}
		/* FILTERS ======================== */
		Set<String> filters=ctx.getFilters();
		out.writeInt(filters.size());
		for(String filter:filters)
			{
			out.writeString(filter);
			}
		/* INFO ======================== */
		Map<String,Object> infoMap=ctx.getAttributes();
		out.writeInt(infoMap.size());
		for(String key:infoMap.keySet())
			{
			out.writeString(key);
			super.writeAttribute(out,infoMap.get(key));
			}
		/* Genotypes ======================== */
		List<Genotype> genotypes=ctx.getGenotypes();
		out.writeInt(genotypes.size());
		for(Genotype g:genotypes)
			{
			this.genotypeBinding.objectToEntry(g, out);
			}
		}
	
	}
