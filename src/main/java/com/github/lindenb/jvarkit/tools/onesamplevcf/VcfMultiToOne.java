/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.onesamplevcf;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/*
 * VcfMultiToOne
 */
@Program(name="vcfmulti2one",description=">Convert VCF with multiple samples to a VCF with one SAMPLE, duplicating variant and adding the sample name in the INFO column",keywords={"vcf","sample"})
public class VcfMultiToOne extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfMultiToOne.class).make();

	@Parameter(names={"-c","--discard_no_call"},description="discard if variant is no-call")
	private boolean discard_no_call = false;
	@Parameter(names={"-r","--discard_hom_ref"},description="discard if variant is hom-ref")
	private boolean discard_hom_ref = false;
	@Parameter(names={"-a","--discard_non_available"},description="discard if variant is not available")
	private boolean discard_non_available = false;
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	

	public static final String DEFAULT_VCF_SAMPLE_NAME="SAMPLE";
	public static final String DEFAULT_SAMPLE_TAGID="SAMPLENAME";
	public static final String DEFAULT_SAMPLE_FILETAGID="SAMPLESOURCE";
	public static final String SAMPLE_HEADER_DECLARATION="VcfMultiToOne.Sample";
	
	public VcfMultiToOne()
		{
		}
	
	/** general utility for program using VCFMulti2One:
	 *  Extract SampleNames
	 */
	static Set<String> extractSampleNames(final VCFHeader header)
		{
		List<String> sample_list =header.getSampleNamesInOrder();
		if(sample_list.size()!=1 || !sample_list.get(0).equals(DEFAULT_VCF_SAMPLE_NAME))
			{
			throw new IllegalArgumentException("Not a VCF produced by VcfMultiToOne");
			}
		Set<String> samples = new TreeSet<String>();
		for(VCFHeaderLine h:header.getMetaDataInInputOrder())
			{
			if(h.getKey().equals(SAMPLE_HEADER_DECLARATION))
				{
				sample_list.add(h.getValue());
				}
			}
		return samples;
		}
	@Override
	public int doWork(final List<String> arguments) {
		VariantContextWriter  out=null;
		Set<String> args= IOUtils.unrollFiles(arguments);
		List<VcfIterator> inputs=new ArrayList<>(args.size()+1);
		List<String> inputFiles=new ArrayList<>(args.size()+1);
		
		
		try
			{
			if(args.isEmpty() && arguments.isEmpty())
				{
				inputs.add(VCFUtils.createVcfIteratorStdin());
				inputFiles.add("stdin");
				}
			else if(args.isEmpty())
				{
				LOG.error("No vcf provided");
				return -1;
				}
			else
				{
				for(String vcfFile: args)
					{
					inputs.add(VCFUtils.createVcfIterator(vcfFile));
					inputFiles.add(VCFUtils.escapeInfoField(vcfFile));
					}
				}
			SAMSequenceDictionary dict=null;
			Set<String> sampleNames=new HashSet<String>();

			Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
			for(VcfIterator in:inputs)
				{
				VCFHeader header = in.getHeader();
				if(dict==null)
					{
					dict = header.getSequenceDictionary();
					}
				else if(header.getSequenceDictionary()==null)
					{
					LOG.error("No Dictionary in vcf");
					return -1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, header.getSequenceDictionary()))
					{
					LOG.error("Not the same dictionary between vcfs");
					return -1;
					}
				metaData.addAll(in.getHeader().getMetaDataInInputOrder());
				sampleNames.addAll(in.getHeader().getSampleNamesInOrder());
				}
			
			
			final Comparator<VariantContext> comparator = (dict==null?
					VCFUtils.createChromPosRefComparator():
					VCFUtils.createTidPosRefComparator(dict)
					);
			//addMetaData(metaData);
			metaData.add(new VCFInfoHeaderLine(
					DEFAULT_SAMPLE_TAGID,1,VCFHeaderLineType.String,
					"Sample Name from multi-sample vcf"
					));
			metaData.add(new VCFInfoHeaderLine(
					DEFAULT_SAMPLE_FILETAGID,1,VCFHeaderLineType.String,
					"Origin of sample"
					));
			
			for(String sample:sampleNames)
				{
				metaData.add(
					new VCFHeaderLine(
					SAMPLE_HEADER_DECLARATION,
					sample));
				}
			
			VCFHeader h2 = new VCFHeader(
					metaData,
					Collections.singleton(DEFAULT_VCF_SAMPLE_NAME)
					);
			
			out= super.openVariantContextWriter(this.outputFile);
			out.writeHeader(h2);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			for(;;)
				{
				if(out.checkError()) break;
				/* get 'smallest' variant */
				VariantContext smallest = null;
				int idx=0;
				int best_idx = -1;
				
				while(idx < inputs.size())
					{
					VcfIterator in= inputs.get(idx);
					if(!in.hasNext())
						{
						CloserUtil.close(in);
						inputs.remove(idx);
						inputFiles.remove(idx);
						}
					else
						{
						VariantContext ctx = in.peek();
						if( smallest==null ||
							comparator.compare(smallest,ctx)>0)
							{
							smallest = ctx;
							best_idx = idx;
							}
						++idx;
						}
					}
				
				if(smallest==null) break;
				
				VariantContext ctx = progress.watch(inputs.get(best_idx).next());
				
				
				if(ctx.getNSamples()==0)
					{
					if(!this.discard_no_call)
						{
						VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						vcb.attribute(DEFAULT_SAMPLE_FILETAGID,inputFiles.get(best_idx));
						vcb.genotypes(GenotypeBuilder.createMissing(DEFAULT_VCF_SAMPLE_NAME,2));
						out.add(vcb.make());
						}
					continue;
					}
				
				for(int i=0;i< ctx.getNSamples();++i)
					{
					Genotype g= ctx.getGenotype(i);
					String sample = g.getSampleName();
					
					if(!g.isCalled() && this.discard_no_call) continue;
					if(!g.isAvailable() && this.discard_non_available) continue;
					if(g.isHomRef() && this.discard_hom_ref) continue;
					
					
					GenotypeBuilder gb=new GenotypeBuilder(g);
					gb.name(DEFAULT_VCF_SAMPLE_NAME);
					
					
					VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.attribute(DEFAULT_SAMPLE_TAGID, sample);
					vcb.attribute(DEFAULT_SAMPLE_FILETAGID,inputFiles.get(best_idx));
					
					
					vcb.genotypes(gb.make());
					out.add(vcb.make());
					}
				}
			progress.finish();
			LOG.debug("done");
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(inputs);
			CloserUtil.close(out);
			}
		}
	
	
	
	public static void main(String[] args)
		{
		new VcfMultiToOne().instanceMainWithExit(args);
		}
	}
