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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class FixVcfMissingGenotypes extends AbstractCommandLineProgram
	{
	private Map<String,Set<File>> sample2bam=new HashMap<>();
	private int minDepth=10;
	
	private FixVcfMissingGenotypes()
		{
		
		}

	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/FixVcfMissingGenotypes";	
		}

	@Override
	public String getProgramDescription() {
		return "after a VCF-merge, read a VCF, look back at some BAMS to tells if the missing genotypes were homozygotes-ref or not-called. If the number of reads is greater than min.depth, then the missing genotypes is said hom-ref.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -f (file) list of indexed BAM path with read Groups");
		out.println(" -d (depth) min depth default:"+this.minDepth);
		super.printOptions(out);
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		Set<File> bamFiles=new HashSet<>();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"f:d:"))!=-1)
			{
			switch(c)
				{
				case 'd': this.minDepth=Integer.parseInt(opt.getOptArg());break;
				case 'f':
					{

					try {
						BufferedReader in=IOUtils.openURIForBufferedReading(opt.getOptArg());
						String line;
						while((line=in.readLine())!=null)
							{
							if(line.startsWith("#") || line.trim().isEmpty()) continue;
							bamFiles.add(new File(line));
							}
						in.close();
						} 
					catch (Exception e)
						{
						error(e);
						return -1;
						}
						
					break;
					}
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		VcfIterator in=null;
		try
			{
			SamReaderFactory srf=SamReaderFactory.makeDefault();
			srf.validationStringency(ValidationStringency.SILENT);
			
			for(File bam: bamFiles)
				{
				info("Reading header for "+bam);
				SamReader reader=srf.open(bam);
				SAMFileHeader header=reader.getFileHeader();
				for(SAMReadGroupRecord g:header.getReadGroups())
					{
					if(g.getSample()==null) continue;
					String sample=g.getSample();
					if(sample.isEmpty()) continue;
					Set<File> set=this.sample2bam.get(sample);
					if(set==null)
						{
						set=new HashSet<>();
						this.sample2bam.put(sample,set);
						}
					set.add(bam);
					}
				reader.close();
				}
			
			if(opt.getOptInd()==args.length)
				{
				in=VCFUtils.createVcfIteratorStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				in=VCFUtils.createVcfIterator(args[opt.getOptInd()]);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			File tmpFile1 = File.createTempFile("fixvcf", ".vcf",getTmpDirectories().get(0));
			File tmpFile2 = File.createTempFile("fixvcf", ".vcf",getTmpDirectories().get(0));
			VCFHeader header=in.getHeader();
			VCFHeader h2=new VCFHeader(header);
			final String FIXED_TAG="FXG";
			h2.addMetaDataLine(new VCFFormatHeaderLine(FIXED_TAG,1,VCFHeaderLineType.Integer,"Genotype was set as homozygous (min depth ="+this.minDepth+")"));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

			for(int i=0;i< header.getNGenotypeSamples();++i)
				{
				int countFixed=0;
				int countNonFixed=0;
				int countTotal=0;
				
				String sample= header.getSampleNamesInOrder().get(i);
				info("Sample: "+sample);
				Set<File> bams = this.sample2bam.get(sample);
				if(bams==null) bams = new HashSet<File>();
				if(bams.isEmpty())
					{
					warning("No bam to fix sample "+sample);
					//don't 'continue' for simplicity
					}
				List<SamReader> samReaders= new ArrayList<>(bams.size());
				for(File bam:bams)
					{
					info("Opening "+bam);
					samReaders.add(srf.open(bam));
					}
				VariantContextWriter w = VCFUtils.createVariantContextWriter(i%2==0?tmpFile1:tmpFile2);	
				w.writeHeader(h2);
				while(in.hasNext())
					{
					VariantContext ctx= in.next();
					countTotal++;
					if(samReaders.isEmpty()) 
						{
						w.add(ctx);
						continue;
						}
					Genotype genotype = ctx.getGenotype(sample);
					if(genotype!=null && genotype.isCalled())
						{	
						w.add(ctx);
						continue;
						}
					int depth=0;
					//get depth for this position
					for(SamReader sr: samReaders)
						{
						SAMRecordIterator iter=sr.query(ctx.getChr(), ctx.getStart(), ctx.getEnd(), false);
						while(iter.hasNext())
							{
							SAMRecord rec=iter.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(rec.getDuplicateReadFlag()) continue;
							if(rec.isSecondaryOrSupplementary()) continue;
							SAMReadGroupRecord rg=rec.getReadGroup();
							if(!sample.equals(rg.getSample())) continue;
							Cigar cigar=rec.getCigar();
							if(cigar==null) continue;
							int refPos=rec.getAlignmentStart();
							for(CigarElement ce:cigar.getCigarElements())
								{
								if(!ce.getOperator().consumesReferenceBases()) continue;
								if( ce.getOperator().consumesReadBases() &&
									refPos<= ctx.getStart() &&
									ctx.getStart() <= refPos+ce.getLength()
									)
									{
									depth++;
									break;
									}
								refPos= ce.getLength();
								}
							}
						iter.close();
						}
					if(depth< this.minDepth)
						{	
						countNonFixed++;
						w.add(ctx);
						continue;
						}
					List<Allele> homozygous=new ArrayList<>(2);
					homozygous.add(ctx.getReference());
					homozygous.add(ctx.getReference());
					GenotypeBuilder gb=new GenotypeBuilder(genotype);
					gb.alleles(homozygous);
					gb.attribute(FIXED_TAG, 1);
					if(header.getFormatHeaderLine("DP")!=null)
						{
						gb.DP(depth);
						}
					
					GenotypesContext gtx=GenotypesContext.copy(ctx.getGenotypes());
					gtx.replace(gb.make());
					
					VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.genotypes(gtx);
					w.add(vcb.make());
					countFixed++;
					}
				w.close();
				in.close();
				//closing BAMS
				for(SamReader r:samReaders) CloserUtil.close(r);
				
				info("done sample "+sample+
						" fixed="+countFixed+
						" not-fixed="+countNonFixed+
						" total="+countTotal+
						" genotypes");
				
				//reopen in
				in = VCFUtils.createVcfIteratorFromFile(i%2==0?tmpFile1:tmpFile2);
				h2= in.getHeader();
				}
			
			VariantContextWriter w = VCFUtils.createVariantContextWriterToStdout();
			w.writeHeader(h2);
			while(in.hasNext())
				{
				w.add(in.next());
				}
			in.close();
			w.close();
			tmpFile1.delete();
			tmpFile2.delete();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	public static void main(String[] args) {
		new FixVcfMissingGenotypes().instanceMainWithExit(args);

	}

}
