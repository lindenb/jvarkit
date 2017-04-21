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
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;


/**
 * @author lindenb
 *
 */
@Program(name="vcfcomm",description="Equivalent of linux comm for VCF")
public class VCFComm extends AbstractVCFCompareBase {
	private final Logger LOG=Logger.build(VCFComm.class).make();

	public VCFComm() 
		{
		}
	@Parameter(names="-a",description="ignore variations present in ALL files")
	private boolean ignore_everywhere=false;
	@Parameter(names="-A",description="only print variations present in ALL files")
	private boolean only_everywhere=false;
	
	@Override
	public int doWork(List<String> args) {
	
		SortingCollection<LineAndFile> variants=null;
		try
			{
			if(args.isEmpty())
				{
				LOG.error("Illegal number of arguments");
				return -1;
				}
			
			Set<String> filenames=new HashSet<String>();
			for(final String filename:args)
				{
				
				if(!filename.endsWith(".list"))
					{
					filenames.add(filename);
					}
				else
					{
					LOG.info("Reading filenames from "+filename);
					BufferedReader in = IOUtils.openURIForBufferedReading(filename);
					String line;
					while((line=in.readLine())!=null)
						{
						if(line.trim().isEmpty() || line.startsWith("#")) continue;
						filenames.add(line);
						}
					in.close();
					}
				}

			
			Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
			
			final LineAndFileComparator posCompare=new LineAndFileComparator();
			
			
			variants=SortingCollection.newInstance(
					LineAndFile.class, new LineAndFileCodec(),
					posCompare,
					super.sortingCollectionArgs.getMaxRecordsInRam(),
					super.sortingCollectionArgs.getTmpDirectories()
					);
			variants.setDestructiveIteration(true);
			
			List<String> newSampleNames=new ArrayList<String>();
			Set<String> sampleSet=new HashSet<String>();
			
			int vcfindex=0;
			for(final String vcffilename:filenames)
				{
				++vcfindex;
				LOG.info("Reading from "+vcffilename);
				Input input=super.put(variants, vcffilename);
				final String sampleName="f"+(vcfindex);
				newSampleNames.add(sampleName);
				metaData.add(new VCFHeaderLine(sampleName,vcffilename));
				
				sampleSet.addAll(input.codecAndHeader.header.getSampleNamesInOrder());
				}
			
			variants.doneAdding();
			
			/** unique sample name, if any */
			String theSampleName=null;
			if(sampleSet.size()==1)
				{
				theSampleName=sampleSet.iterator().next();
				LOG.info("Unique sample name is "+theSampleName);
				}
			
			metaData.add(new VCFHeaderLine(getClass().getSimpleName(),"version:"+getVersion()+" command:"+getProgramCommandLine()));
			metaData.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String,
					"Genotype"+(theSampleName==null?"":" for Sample "+theSampleName)));
			metaData.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Depth"));
			metaData.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Integer, "Qual"));
			metaData.add(new VCFFormatHeaderLine("QUAL", 1, VCFHeaderLineType.Float, "VCF QUAL column"));
			metaData.add(new VCFFilterHeaderLine("diffCall", "Variant NOT called in all columns"));

			
			if(theSampleName!=null)
				{
				metaData.add(new VCFHeaderLine("UniqSample",theSampleName));
				metaData.add(new VCFFilterHeaderLine("diffGT", "Genotype difference for sample "+theSampleName));
				}
			VCFHeader header=new VCFHeader(
					metaData,
					newSampleNames
					);

			
			VariantContextWriter w= super.openVariantContextWriter(super.outputFile);
			w.writeHeader(header);
			List<LineAndFile> row=new ArrayList<LineAndFile>(super.inputs.size());
			
			
			CloseableIterator<LineAndFile> iter=variants.iterator();
			
			for(;;)
				{
				LineAndFile rec=null;
				if(iter.hasNext())
					{
					rec=iter.next();
					}
			
				if(rec==null || (!row.isEmpty() && posCompare.compare(row.get(0),rec)!=0))
					{
					if(!row.isEmpty())
						{
						Set<GenotypeType> typeGenotypes=new HashSet<GenotypeType>();
						VariantContext first=row.get(0).getContext();
						
						Set<Allele> alleles=new HashSet<Allele>();
						alleles.add(first.getReference());
						for(LineAndFile laf:row)
							{							
							alleles.addAll(laf.getContext().getAlleles());
							}
						
						
						VariantContextBuilder b=new VariantContextBuilder(
								getClass().getName(),
								first.getContig(),
								first.getStart(),
								first.getEnd(),
								alleles
								);
						Set<String> ids=new TreeSet<String>();
						//build genotypes
						List<Genotype> genotypes=new ArrayList<Genotype>();
						for(LineAndFile laf:row)
							{
							//alleles for this genotype
							List<Allele> galleles=new ArrayList<Allele>();
							if(theSampleName==null)
								{
								galleles.add(first.getReference());
								galleles.add(first.getReference());
								}
							else
								{
								Genotype g0=laf.getContext().getGenotype(theSampleName);
								if(g0==null) throw new IllegalStateException();
								for(Allele a:g0.getAlleles())
									{
									galleles.add(a);
									}
								}
							GenotypeBuilder gb=new GenotypeBuilder();
							//gb.DP(1);
							gb.alleles(galleles);
							gb.name(newSampleNames.get(laf.fileIdx));
							//gb.GQ(1);
							if(laf.getContext().hasLog10PError())
								{
								gb.attribute("QUAL", laf.getContext().getPhredScaledQual());
								}
							
							Genotype genotype=gb.make();
							typeGenotypes.add(genotype.getType());
							genotypes.add(genotype);
							
							if(laf.getContext().hasID())
								{
								ids.add(laf.getContext().getID());
								}
							}
						b.genotypes(genotypes);
						if(!ids.isEmpty())
							{
							StringBuilder sw=new StringBuilder();
							for(String s:ids)
								{
								if(sw.length()!=0) sw.append(",");
								sw.append(s);
								}
							b.id(sw.toString());
							}
						Set<String> filters=new HashSet<String>();
						
						
						if(theSampleName!=null && typeGenotypes.size()!=1)
							{
							filters.add("diffGT");
							}
						
						
						boolean print=true;
						if(row.size()==super.inputs.size() && ignore_everywhere)
							{
							print=false;
							}
						if(row.size()!=super.inputs.size())
							{
							filters.add("diffCall");
							if(only_everywhere)
								{
								print=false;
								}
							}
						
						b.filters(filters);
						
						if(print)
							{
							w.add(b.make());
							}
						row.clear();
						}
					if(rec==null) break;
					}
				
				row.add(rec);
				}
			iter.close();
			
			w.close();

			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			try
				{
				if(variants!=null) variants.cleanup();
				}
			catch(Exception err)
				{
				}
			}
		}
	public static void main(String[] args)
		{
		new VCFComm().instanceMainWithExit(args);
		}
	}
