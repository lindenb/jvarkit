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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
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
	
	
	private PrintWriter writer=null;
	
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
	
	private static class SkatCaller implements Callable<SkatCallerResult>
		{
		private final File vcfFile;
		private final Interval interval;
		private final Set<Pedigree.Person> samples;
		private final SkatFactory.SkatExecutor skatExec;
		SkatCaller( final File vcfFile,
					final Interval interval,
					final Set<Pedigree.Person> samples,
					final SkatFactory.SkatExecutor skatExec
					) {
			this.vcfFile = vcfFile;
			this.interval = interval;
			this.samples= samples;
			this.skatExec = skatExec;
			}
		
		@Override
		public SkatCallerResult call() throws Exception {
			final SkatCallerResult result = new SkatCallerResult();
			result.interval = this.interval;
			VCFFileReader vcfFileReader = null;
			CloseableIterator<VariantContext> iter=null;
			try
				{
				final List<VariantContext> variants = new ArrayList<>();
				vcfFileReader=new VCFFileReader(this.vcfFile,true);
				iter = vcfFileReader.query(this.interval.getName(), this.interval.getStart(), this.interval.getEnd());
				while(iter.hasNext())
					{
					final VariantContext ctx = iter.next();
					if(ctx.isFiltered()) continue;
					if(ctx.getNAlleles()!=2) continue;
					variants.add(ctx);
					}
				vcfFileReader.close();
				vcfFileReader=null;
				iter.close();
				iter=null;
				result.nVariants = variants.size();
				if(variants.isEmpty())
					{
					result.pvalue=1;
					return result;
					}
				
				final SkatFactory.SkatResult skatResult = this.skatExec.execute(variants, this.samples);
				if(skatResult.isError())
					{
					result.error_msg=skatResult.getMessage();
					return result;
					}
				else
					{
					result.pvalue = skatResult.getPValue();
					}
				return result;
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				result.error_msg = "Exception: "+ err.getMessage() ;
				return result;
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
			final LinkedBlockingQueue<Runnable> queue = new LinkedBlockingQueue<Runnable>(this.nJobs);
			final ExecutorService executorService =  new ThreadPoolExecutor(
					   this.nJobs, this.nJobs,
			              0L, TimeUnit.MILLISECONDS,
			              queue
			              );
			final List<Future<SkatCallerResult>> futureResults= new ArrayList<>();
			
			final Runnable dumpResults = () -> {
				try {
					if(!futureResults.isEmpty())
						{
						int i=0;
						while(i<futureResults.size())
							{
					 		final Future<SkatCallerResult> ft=futureResults.get(i);
							if(ft.isCancelled())
								{
								throw new RuntimeException("Task was canceled.");
								}
							else if(ft.isDone())
								{
								final SkatCallerResult rez= futureResults.remove(i).get();
								synchronized (this.writer) {
									this.writer.println(rez.toString());
									}
								}
							else
								{
								i++;
								}
							}
						}
					}
				catch(final ExecutionException|InterruptedException err)
					{
					throw new RuntimeException(err);
					}
				};

			final Set<Pedigree.Person> samples= new HashSet<>(
					pedigree.getPersons()
					);
			samples.removeIf(I->!(I.isAffected() || I.isUnaffected()) || !header.getSampleNamesInOrder().contains(I.getId()));
			
			   
			this.writer = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			for(final SAMSequenceRecord ssr:dict.getSequences())
				{
				if(!this.limit_contigs.isEmpty() && !this.limit_contigs.contains(ssr.getSequenceName())) {
					LOG.warning("skipping " + ssr.getSequenceName());
					continue;
					}
				int x=1;
				while(x< ssr.getSequenceLength())
					{
					while(queue.size()>= this.nJobs)
						{
   						dumpResults.run();
						}
					SkatCaller callable = new SkatCaller(
							vcfFile,
							new Interval(ssr.getSequenceName(),x,Math.min(ssr.getSequenceLength(), x+this.contigWinLength)),
							samples,
							this.skat.build()
							);
					final Future<SkatCallerResult> rez = executorService.submit(callable);
					futureResults.add(rez);

					
					x+=this.contigWinShift;
					}
				}
			executorService.shutdown();
			executorService.awaitTermination(365, TimeUnit.DAYS);
			LOG.warning("last dump");
			dumpResults.run();
			
			this.writer.flush();
			this.writer.close();
			this.writer=null;
			
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
			CloserUtil.close(this.writer);
			}
		}
	public static void main(String[] args)
		{
		new VcfSkatSlidingWindow().instanceMainWithExit(args);
		}
	}
