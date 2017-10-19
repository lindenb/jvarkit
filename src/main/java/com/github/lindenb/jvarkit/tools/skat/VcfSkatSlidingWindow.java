/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.skat;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

END_DOC

 */
@Program(
		name="vcfskatslidingwindow",
		description="SkatFactory Over genome using a sliding window.",
		keywords={"vcf","pedigree","skat","burden"}
		)
public class VcfSkatSlidingWindow extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfSkatSlidingWindow.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION+" If not defined, I will try to extract the pedigree from  the VCFheader.")
	private File pedigreeFile=null;
	@Parameter(names={"-C","--contig"},description="limit to this contig(s)")
	private Set<String> limit_contigs = new HashSet<>();
	@Parameter(names={"-j","--jobs"},description="When -exec is specified, use <n> jobs. A value lower than 1 means use all procs available. ")
	private int nJobs = 1;
	@Parameter(names={"--contigWinLength"},description="window size when splitting per contig")
	private int contigWinLength=1000;
	@Parameter(names={"--contigWinShift"},description="window shift when splitting per contig")
	private int contigWinShift=500;
	@ParametersDelegate
	private SkatFactory skat = new SkatFactory();
	
	
	
	private static class SkatCallerResult
		{
		Interval interval;
		String error_msg= null;
		double pvalue=1;
		int nVariants=0;
		
		@Override
		public String toString() {
			return 
					(error_msg==null?"":"#")+
					interval.getContig()+"\t"+(interval.getStart()-1)+"\t"+interval.getEnd()+
					"\t"+nVariants+"\t"+pvalue+
					(this.error_msg==null?"":"\t"+this.error_msg)
					;
			}
		}
	
	private class SkatCaller implements Callable<List<SkatCallerResult>>
		{
		private final File vcfFile;
		private final Interval fromTo;
		private final Set<Pedigree.Person> samples;
		private final SkatFactory.SkatExecutor skatExec;
		SkatCaller( final File vcfFile,
					final Interval fromTo,
					final Set<Pedigree.Person> samples,
					final SkatFactory.SkatExecutor skatExec
					) {
			this.vcfFile = vcfFile;
			this.fromTo = fromTo;
			this.samples= samples;
			this.skatExec = skatExec;
			}
		
		@Override
		public List<SkatCallerResult> call() throws Exception {
			VCFFileReader vcfFileReader = null;
			CloseableIterator<VariantContext> iter=null;
			final List<SkatCallerResult> results = new ArrayList<>();
			int x=1;
			try
				{
				final List<VariantContext> variants = new ArrayList<>();
				vcfFileReader=new VCFFileReader(this.vcfFile,true);
				
				while(x<this.fromTo.getEnd())
					{
					if(x < this.fromTo.getStart()) {
						x += VcfSkatSlidingWindow.this.contigWinShift;
						continue;
						}
					final SkatCallerResult result = new SkatCallerResult();
					
					
					result.interval = new Interval(
							this.fromTo.getContig(),
							x,
							x+VcfSkatSlidingWindow.this.contigWinLength
							);

					iter = vcfFileReader.query(result.interval.getName(), result.interval.getStart(), result.interval.getEnd());
					while(iter.hasNext())
						{
						final VariantContext ctx = iter.next();
						if(ctx.isFiltered()) continue;
						if(ctx.getNAlleles()!=2) continue;
						variants.add(ctx);
						}
					iter.close();
					iter=null;
					result.nVariants = variants.size();
					if(variants.isEmpty())
						{
						result.pvalue=1;
						}
					else
						{
						final SkatFactory.SkatResult skatResult = this.skatExec.execute(variants, this.samples);
						if(skatResult.isError())
							{
							result.error_msg=skatResult.getMessage();
							}
						else
							{
							result.pvalue = skatResult.getPValue();
							}
						}
					results.add(result);
					}
				vcfFileReader.close();
				vcfFileReader=null;
				return results;
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				throw new RuntimeException(err);
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(vcfFileReader);
				}
			}
		}
		
	@Override
	public int doWork(final List<String> args)
		{
		if(this.nJobs<1)
			{
			this.nJobs = Math.max(1, Runtime.getRuntime().availableProcessors());
			LOG.info("setting njobs to "+this.nJobs);
			}
		PrintWriter writer=null;
		VcfIterator r=null;
		try {	
			final VCFHeader header;
			final SAMSequenceDictionary dict;
			final File vcfFile= new File(oneAndOnlyOneFile(args));
			try (final VCFFileReader vr=new VCFFileReader(vcfFile, true))
				{
				header=vr.getFileHeader();
				dict=header.getSequenceDictionary();
				}
			if(dict==null || dict.isEmpty())
				{
				throw new JvarkitException.VcfDictionaryMissing(vcfFile);
				}
			if(!this.limit_contigs.isEmpty())
				{
				if(this.limit_contigs.stream().anyMatch(C->dict.getSequence(C)==null))
					{
					LOG.error("user contig missing in vcf dictionary.");
					return -1;
					}
				}
			final Pedigree pedigree;
			if(this.pedigreeFile!=null)
				{
				pedigree = new Pedigree.Parser().parse(this.pedigreeFile);
				}
			else
				{
				pedigree = new Pedigree.Parser().parse(header);
				}
			final Set<Pedigree.Person> samples= new HashSet<>(
					pedigree.getPersons()
					);
			samples.removeIf(I->!(I.isAffected() || I.isUnaffected()) || !header.getSampleNamesInOrder().contains(I.getId()));
			
			   
			writer = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			for(final SAMSequenceRecord ssr:dict.getSequences())
				{
				if(!this.limit_contigs.isEmpty() && !this.limit_contigs.contains(ssr.getSequenceName())) {
					LOG.warning("skipping " + ssr.getSequenceName());
					continue;
					}
				
				final ExecutorService executorService =  new ThreadPoolExecutor(
						   this.nJobs, this.nJobs,
				              0L, TimeUnit.MILLISECONDS,
				              new LinkedBlockingDeque<>(this.nJobs)
				              );				
				
				final List<Future<List<SkatCallerResult>>> results = new ArrayList<>(this.nJobs);
				for(int i=0;i< this.nJobs;i++)
					{
					final int winLen = Math.max(1,ssr.getSequenceLength()/this.nJobs);
					final SkatCaller caller = new SkatCaller(
							vcfFile,
							new Interval(ssr.getSequenceName(),i*winLen,(i+1)*winLen),
							samples,
							this.skat.build()
							)
							;
					results.add(executorService.submit(caller));
					}
				executorService.awaitTermination(365, TimeUnit.DAYS);
				executorService.shutdown();
				for(final Future<List<SkatCallerResult>> fl:results)
					{
					try {
						for(final SkatCallerResult scr: fl.get())
							{
							writer.println(scr.toString());
							}
						}
					catch(Exception err)
						{
						writer.close();
						throw new RuntimeException(err);
						}
					}				
				}
				
			writer.flush();
			writer.close();
			writer=null;
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(writer);
			}
		}
	public static void main(String[] args)
		{
		new VcfSkatSlidingWindow().instanceMainWithExit(args);
		}
	}
