/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfMultiToOneAllele
	extends AbstractVcfMultiToOneAllele
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfMultiToOneAllele.class);

	 public VcfMultiToOneAllele()
		{
		}
	 
	@Override
	/* public for knime */
	public Collection<Throwable> doVcfToVcf(String inputName,
			VcfIterator in, VariantContextWriter out) throws IOException {
		final String TAG="VCF_MULTIALLELIC_SRC";
		final List<String> noSamples=Collections.emptyList();
	
		final VCFHeader header=in.getHeader();
		final List<String> sample_names=header.getSampleNamesInOrder();
		final Set<VCFHeaderLine> metaData=new HashSet<>(header.getMetaDataInInputOrder());
		addMetaData(metaData);		
		metaData.add(new VCFInfoHeaderLine(TAG, 1, VCFHeaderLineType.String,
				"The variant was processed with VcfMultiAlleleToOneAllele and contained the following alleles."));
		VCFHeader h2;
		
		if(!super.print_samples)
			{
			h2 = new VCFHeader(
					metaData,
					noSamples
					);
			}
		else
			{
			h2 = new VCFHeader(
					metaData,
					sample_names
					);
			}
		SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header);
		out.writeHeader(h2);
		while(in.hasNext())
			{
			final VariantContext ctx=progess.watch(in.next());
			final List<Allele> alleles = new ArrayList<>(ctx.getAlternateAlleles());
			if(alleles.isEmpty())
				{
				LOG.warn("Remove no ALT variant:"+ctx);
				continue;
				}
			else if(alleles.size()==1)
				{
				if(!print_samples)
					{
					final VariantContextBuilder vcb = super.getVariantContextBuilderFactory().newVariantContextBuilder(ctx);
					vcb.noGenotypes();
					out.add(vcb.make());
					}
				else
					{
					out.add(ctx);
					}
				}
			else
				{
				//Collections.sort(aioulleles); don't sort , for VCFHeaderLineCount.A
				final Map<String,Object> attributes = ctx.getAttributes();
				final StringBuilder sb=new StringBuilder();
				for(int i=0;i< alleles.size();++i)
					{
					if(sb.length()>0) sb.append("|");
					sb.append(alleles.get(i).getDisplayString());
					}
				final String altAsString= sb.toString();
				for(int i=0;i< alleles.size();++i)
					{
					final Allele the_allele = alleles.get(i);

					final VariantContextBuilder vcb = super.getVariantContextBuilderFactory().newVariantContextBuilder(ctx);
					vcb.alleles(Arrays.asList(ctx.getReference(),the_allele));
					
					for(final String attid:attributes.keySet())
						{
						final VCFInfoHeaderLine info = header.getInfoHeaderLine(attid);
						if(info==null) throw new IOException("Cannot get header INFO tag="+attid);
						if(info.getCountType()!=VCFHeaderLineCount.A) continue;
						final Object o = 	attributes.get(attid);
						if(!(o instanceof List)) {
							final String msg="For INFO tag="+attid+" got "+o.getClass()+" instead of List in "+ctx;
							if(super.rmErrorAttributes)
								{
								LOG.warn("remove this attribute : "+msg);
								vcb.rmAttribute(attid);
								continue;
								}
							else
								{
								throw new IOException(msg);
								}				
							}
						@SuppressWarnings("rawtypes")
						final List list = (List)o;
						if(alleles.size()!=list.size()) {
							final String msg= ctx.getContig()+":"+ctx.getStart()+" : For INFO tag="+attid+" got "+alleles.size()+" ALT, incompatible with "+list.toString();
							if(super.rmErrorAttributes)
								{
								LOG.warn("remove this attribute : "+msg);
								vcb.rmAttribute(attid);
								continue;
								}
							else
								{
								throw new IOException(msg);
								}
							}
						else
							{	
							vcb.attribute(attid, list.get(i));	
							}
						}
					
					vcb.attribute(TAG,altAsString);
					
					if(!print_samples)
						{
						vcb.noGenotypes();
						}
					else
						{
						final List<Genotype> genotypes=new ArrayList<>(sample_names.size());
						
						for(final String sampleName: sample_names)
							{							
							final Genotype g= ctx.getGenotype(sampleName);
							if(!g.isCalled() || g.isNoCall() )
								{
								genotypes.add(g);
								continue;
								}
							
							
							final GenotypeBuilder gb =new GenotypeBuilder(g);
							final List<Allele> galist = new ArrayList<>(g.getAlleles());
							
							if(galist.size()>0)
								{
								boolean replace=false;
								for(int y=0;y< galist.size();++y)
									{
									final Allele ga = galist.get(y);
									if(ga.isSymbolic()) throw new RuntimeException("How should I handle "+ga);
									if(!(ga.isNoCall() || 
										 ga.equals(ctx.getReference()) ||
										 ga.equals(the_allele)))
										{
										replace=true;
										galist.set(y, ctx.getReference());
										}
									}
								if(replace)
									{
									gb.reset(true);/* keep sample name */
									gb.alleles(galist);
									}
								}
							genotypes.add(gb.make());
							}
						
						vcb.genotypes(genotypes);
						}
					
					out.add(vcb.make());
					}
				}
			}
		progess.finish();
		return RETURN_OK;
		}

	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	
	public static void main(String[] args)
		{
		new VcfMultiToOneAllele().instanceMainWithExit(args);
		}
	
	}
