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


*/
package com.github.lindenb.jvarkit.tools.epistasis;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.tools.vcflist.VcfList;
import com.github.lindenb.jvarkit.tools.vcflist.VcfOffsetsIndexFactory;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
/**
BEGIN_DOC



END_DOC
 */
public class VcfEpistatis01 extends Launcher {
	private static final Logger LOG = Logger.build(VcfEpistatis01.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-p","--pedigree"},description=Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile = null;
	@Parameter(names={"--memory"},description="Load all variants in memory")
	private  boolean load_variants_in_memory=false;
	@Parameter(names={"-j","--jobs"},description="Number of parallel jobs.")
	private  int number_of_jobs =1;
	@Parameter(names={"-start","--start"},description="Specify start index in variant list. (for parallelisation)")
	private  int start_index_at=0;
	@Parameter(names={"-jexl","--jexl"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = (CTX)->true;
	@Parameter(names={"-score","--score"},description="[20171220] Output score. Default is printing the VCF.")
	private boolean output_score = false;
	
	
	private static final Function<Long,Integer> CTRLS_nAlt2score=(N)->{switch(N.intValue()){
		case 0: return 0;
		case 1 : return -10;
		case 2: //cont
		default: return -30;
		}};
	private static final Function<Long,Integer> CASES_nAlt2score=(N)->{switch(N.intValue()){
		case 0: return 0;
		case 1 : return 10;
		case 2: //cont
		default: return 30;
		}};
	
		
		
	private static class Result
		{
		final VariantContext ctx1;
		final Allele a1;
		final int idx1;
		final VariantContext ctx2;
		final Allele a2;
		final int idx2;
		final double score;
		Result(final VariantContext ctx1,final Allele a1,int idx1,
			   final VariantContext ctx2, final Allele a2,int idx2,
			   double score
				)
			{
			this.ctx1 = ctx1;
			this.a1 = a1;
			this.idx1 = idx1;
			this.ctx2 = ctx2;
			this.a2 = a2;
			this.idx2 = idx2;
			this.score=score;
			}
		@Override
		public String toString() {
			return 
					ctx1.getContig()+":"+ctx1.getStart()+":"+ctx1.getReference()+"/"+a1+"["+idx1+"] | "+
					ctx2.getContig()+":"+ctx2.getStart()+":"+ctx2.getReference()+"/"+a2+"["+idx2+"] | "+
					score;
			}
		}
	
	
	private static  class Runner implements Callable<Result>
		{
		private final List<VariantContext> variants;
		private final int caseIndexes[];
		private final int ctrlIndexes[];
		private final int startIndex;
		private Result result = null;
		private long duration=0L;
		Runner(
				final List<VariantContext> variants,
				final int startIndex,
				final int[] caseIndexes,
				final int[] ctrlIndexes
				)
			{
			this.variants = variants;
			this.startIndex = startIndex;
			this.caseIndexes = caseIndexes;
			this.ctrlIndexes = ctrlIndexes;
			}
		
		private double score(
				final VariantContext ctx,
				final Allele alt,
				final int samples_indexes[],
				final Function<Long,Integer> nAlt2score
				) 
			{
			double score_a1 = 0;
			
			for(final int  idx : samples_indexes) {
				final Genotype g = ctx.getGenotype(idx);
				if(g==null || g.isFiltered()) continue;
				score_a1 +=nAlt2score.apply( g.getAlleles().stream().filter(A->A.equals(alt)).count());
				}
			return score_a1;
			}
		@Override
		public Result call() throws Exception {
			final VariantContext ctx1 = this.variants.get(this.startIndex);
			final long startup = System.currentTimeMillis();

			int i = this.startIndex + 1;
			
			while(i< this.variants.size())
				{
				final VariantContext ctx2 = this.variants.get(i);
				for(final Allele a1 : ctx1.getAlleles())//faster then getAlternateAlleles
					{
					if(a1.isReference()) continue;
					double score_a1 = 0;
					score_a1 += score(ctx1,a1,this.ctrlIndexes,CTRLS_nAlt2score);
					score_a1 += score(ctx1,a1,this.caseIndexes,CASES_nAlt2score);
					
					for(final Allele a2 : ctx2.getAlleles())//faster then getAlternateAlleles
						{
						if(a2.isReference()) continue;
						double score_a2 = score_a1;
						score_a2 += score(ctx2,a2,this.ctrlIndexes,CTRLS_nAlt2score);
						score_a2 += score(ctx2,a2,this.caseIndexes,CASES_nAlt2score);
						if(this.result == null || this.result.score< score_a2)
							{
							final Result nr = new Result(ctx1,a1,this.startIndex,ctx2,a2,i,score_a2);
							
							
							this.result = nr;
							}
						
						}
					}
				i++;
				}
			if(this.variants instanceof VcfList)
				{
				CloserUtil.close(VcfList.class.cast(this.variants));
				}
			this.duration = System.currentTimeMillis() - startup;
			LOG.info("index ["+startIndex+"] That took "+(duration/1000f)+" seconds.");
			return this.result;
			}
		}
	
	public VcfEpistatis01()
		{
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.number_of_jobs<1) {
			LOG.error("bad number of jobs");
			return -1;
			}
		try
			{
			final int variantsCount;
			final List<VariantContext> inMemoryVariants;
			final File vcfFile = new File(oneAndOnlyOneFile(args));
			final File tmpIndexFile;
			
			if(vcfFile.equals(this.outputFile))
				{
				LOG.error("input == output");
				return -1;
				}
			
			VCFFileReader vcfFileReader = new VCFFileReader(vcfFile,false);
			final VCFHeader header =  vcfFileReader.getFileHeader();

			
			final Pedigree pedigree;
			if(this.pedigreeFile!=null)
				{
				pedigree = new Pedigree.Parser().parse(this.pedigreeFile);
				}
			else
				{
				pedigree = new Pedigree.Parser().parse(header);
				}
			
			
			pedigree.verifyPersonsHaveUniqueNames();
			final Map<String,Integer> sample2index = header.getSampleNameToOffset();
			
			final  int caseIndexes[] = pedigree.getAffected().stream().
					filter(P->sample2index.containsKey(P.getId())).
					mapToInt(P->sample2index.get(P.getId())).
					sorted().
					toArray();

			
			final int ctrlIndexes[] = new ArrayList<>(pedigree.getUnaffected()).stream().
					filter(P->sample2index.containsKey(P.getId())).
					mapToInt(P->sample2index.get(P.getId())).
					sorted().
					toArray();
			
				
			if( caseIndexes.length==0 || ctrlIndexes.length==0 )
					{
					LOG.error("empty ped or no case/ctrl");
					vcfFileReader.close();
					return -1;
					}

			
			if(this.load_variants_in_memory) {
				LOG.info("loading variants in memory");
				tmpIndexFile = null;
				final CloseableIterator<VariantContext> iter2=vcfFileReader.iterator();
				inMemoryVariants =  Collections.unmodifiableList(iter2.stream().
						filter(this.variantFilter).
						filter(V->V.getGenotypes().stream().filter(G->G.isCalled()).count()>0).//should fix https://github.com/samtools/htsjdk/issues/1026 ?
						collect(Collectors.toList())
						);
				variantsCount = inMemoryVariants.size();
				iter2.close();
				}
			else
				{
				tmpIndexFile = File.createTempFile("epistatsis",VcfOffsetsIndexFactory.INDEX_EXTENSION);
				tmpIndexFile.deleteOnExit();
				new VcfOffsetsIndexFactory().
					setLogger(LOG).
					setPredicate(variantFilter).
					indexVcfFile(vcfFile,tmpIndexFile);
				final VcfList tmpList = VcfList.fromFile(vcfFile, tmpIndexFile);
				variantsCount = tmpList.size();
				tmpList.close();
				inMemoryVariants = null;
				}
			

			
			vcfFileReader.close();
			LOG.info("Number of variants: "+variantsCount);
			
			
			Result bestResult =null;
			int x= this.start_index_at;
			while(x+1 < variantsCount)
				{
				final List<Runner> runners = new Vector<>(this.number_of_jobs);
				while(x+1 < variantsCount && runners.size() < this.number_of_jobs)
					{
					LOG.info("starting "+x+"/"+variantsCount);
					runners.add(new Runner(
							inMemoryVariants == null? 
									VcfList.fromFile(vcfFile,tmpIndexFile):
									new Vector<>(inMemoryVariants)
							,x,
							caseIndexes,
							ctrlIndexes
							)
							);
					++x;
					}
				final ExecutorService execSvc;
				if(this.number_of_jobs==1)
					{
					execSvc = null;
					}
				else
					{
					execSvc = Executors.newFixedThreadPool(this.number_of_jobs);;
					}

				if(this.number_of_jobs==1)
					{
					runners.get(0).call();
					}
				else
					{
					execSvc.invokeAll(runners);
					}
					
				if(execSvc!=null) {
					execSvc.shutdown();
					execSvc.awaitTermination(10000L, TimeUnit.DAYS);
					execSvc.shutdownNow();
				}
				runners.stream().mapToLong(R->R.duration).min().ifPresent(D->{
					LOG.info("That took "+ (D/1000f)+" seconds.");
					});
				
				
				for(final Runner r: runners)
					{
					final Result rez=r.result;
					if(rez==null) continue;
					if(bestResult==null || bestResult.score<rez.score)
						{
						bestResult =rez;
						
						if(this.output_score) {
							final PrintWriter pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
							pw.println(bestResult.score+ "\t"+bestResult.toString());
							pw.flush();
							pw.close();
							}
						else
							{
							final VariantContextWriter w = openVariantContextWriter(this.outputFile);
							final VCFHeader header2= new VCFHeader(header);
							header2.addMetaDataLine(new VCFHeaderLine(VcfEpistatis01.class.getName(),bestResult.toString()));
							w.writeHeader(header2);
							w.add(bestResult.ctx1);
							w.add(bestResult.ctx2);
							w.close();
							}
						}
					}
				LOG.info("best: "+bestResult);
				}
			if(tmpIndexFile!=null) tmpIndexFile.delete();
			
			return 0;
			}
		catch(final Exception err)
			{
			err.printStackTrace();
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	 	
	
	public static void main( String[] args)
		{
		//args=new String[] {"-jexl","vc.getStart()%105==0","-j","2","--memory","-o","/home/lindenb/jeter2.vcf","/home/lindenb/jeter.vcf"};
		new VcfEpistatis01().instanceMainWithExit(args);
		}

}
