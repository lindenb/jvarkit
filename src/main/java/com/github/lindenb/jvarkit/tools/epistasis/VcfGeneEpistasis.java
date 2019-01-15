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

*/
package com.github.lindenb.jvarkit.tools.epistasis;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;


import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.tools.skat.SkatFactory;
import com.github.lindenb.jvarkit.tools.skat.SkatFactory.SkatExecutor;
import com.github.lindenb.jvarkit.tools.skat.SkatFactory.SkatResult;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
/**

BEGIN_DOC

END_DOC
*/
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

@Program(name="vcfgeneepistasis",
		description="Burden: gene 1 vs gene 2",
		keywords={"vcf","burden"},
		generate_doc = false
		)
public class VcfGeneEpistasis
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfGeneEpistasis.class).make();

	@Parameter(names={"-o","--output"},description="output file. extension can be vcf or vcf.gz to produce VCF, else will be a plain text report",required=true)
	private File outputFile = null;
	@Parameter(names={"-B","-geneBed","--geneBed"},description="gene bed file. If doesn't exist, will be created.",required=true)
	private File geneBed = null;
	@Parameter(names={"-p","--pedigree"},description=Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile = null;
	@Parameter(names={"-bi","--begin-index"},description="Begin index inclusive. 0-based. For parallelisation.")
	private int user_begin_index=0;
	@Parameter(names={"-ei","--end-index"},description="End index exclusive. 0-based. For parallelisation. -1 : till the end ")
	private int user_end_index=-1;

	
	@ParametersDelegate
	private SkatFactory skatFactory = new SkatFactory();
	

	private VCFFileReader vcfFileReader = null;
	private VcfTools vcfTools = null;
	private final Map<String,Pedigree.Person> id2samples = new HashMap<>();
	
	private static class MergedList<T> extends AbstractList<T>
		{
		final List<T> L1;
		final List<T> L2;
		MergedList(final List<T> L1,final List<T> L2) {
			this.L1 = L1;
			this.L2 = L2;
			}
		public int size() { return L1.size()+L2.size();}
		@Override
		public T get(int index) {
			return index<L1.size()?
					L1.get(index):
					L2.get(index-L1.size())
					;
			}
		}
	
	
	public VcfGeneEpistasis()
		{
		}
	
	private Set<String> getGenes(final VariantContext ctx)
		{
		return this.vcfTools.getAnnPredictions(ctx).stream().
				map(P->P.getGeneId()).filter(P->!StringUtil.isBlank(P)).
				collect(Collectors.toSet());
		}
	
	
	
	private boolean accept(final VariantContext ctx) 
		{
		if(ctx.isFiltered()) return false;
		if(getGenes(ctx).isEmpty()) return false;
		return true;
		}
	
	private final Double eval(final SkatExecutor executor,final List<VariantContext> variants) {
		final SkatResult rez=executor.execute(variants, this.id2samples.values());
		if(rez.isError()) return null;
		return rez.getPValue();
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.geneBed==null) {
			LOG.error("gene file bed undefined");
			return -1;
			}
		if(this.outputFile==null) {
			LOG.error("output file undefined");
			return -1;
			}
		CloseableIterator<VariantContext> iter=null;
		try
			{
			final File vcfFile = new File(oneAndOnlyOneFile(args));
			this.vcfFileReader = new VCFFileReader(vcfFile,true);
			final VCFHeader header = this.vcfFileReader.getFileHeader();
			final Pedigree pedigree;
			if(this.pedigreeFile!=null)
				{
				pedigree = new Pedigree.Parser().parse(this.pedigreeFile);
				}
			else
				{
				pedigree = new Pedigree.Parser().parse(header);
				}
			
			if(pedigree==null || 
				pedigree.isEmpty() ||
				!pedigree.hasAffected() ||
				!pedigree.hasUnaffected()
				)
				{
				LOG.error("empty ped or no case/ctrl");
				return -1;
				}
			pedigree.verifyPersonsHaveUniqueNames();
			
			for(final Pedigree.Person p: pedigree.getPersons().stream().
					filter(P->P.isAffected() || P.isUnaffected()).
					filter(P->header.getSampleNamesInOrder().contains(P.getId())).
					collect(Collectors.toSet()))
					{
					this.id2samples.put(p.getId(), p);
					}
					
			
			this.vcfTools = new VcfTools(header);
			List<Interval> geneList;
			if(!this.geneBed.exists())
				{
				final Map<String, Interval> gene2interval = new HashMap<>(50000);
				LOG.info("building gene file"+this.geneBed);
				iter = this.vcfFileReader.iterator();
				//iter = this.vcfFileReader.query("chr3",1,300_000_000);
				final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
				while(iter.hasNext())
					{
					final VariantContext ctx= progress.watch(iter.next());
					if(!accept(ctx)) continue;
					for(final String geneName : getGenes(ctx))
						{
						final Interval old= gene2interval.get(geneName);
						if(old==null) {
							
							gene2interval.put(geneName, new Interval(
									ctx.getContig(),
									ctx.getStart(),
									ctx.getEnd(),
									false,
									geneName
									));
							LOG.info("adding "+geneName+". number of genes: "+gene2interval.size());
							}
						else if(!old.getContig().equals(ctx.getContig()))
							{
							LOG.error("boum :" +geneName +": on chrom "+ctx.getContig()+" vs "+old);
							return -1;
							}
						else
							{
							gene2interval.put(geneName, new Interval(ctx.getContig(),
									Math.min(ctx.getStart(),old.getStart()),
									Math.max(ctx.getEnd(),old.getEnd()),
									false,
									geneName
									));
							}
						}
					}
				iter.close();iter=null;
				progress.finish();
				
				geneList=new ArrayList<>(gene2interval.values());
				PrintWriter pw=new PrintWriter(this.geneBed);
				for(final Interval g:geneList)
					{
					pw.println(g.getContig()+"\t"+(g.getStart()-1)+"\t"+(g.getEnd())+"\t"+g.getName());
					}
				pw.flush();
				pw.close();
				pw=null;
				}
			else
				{
				BedLineCodec codec= new BedLineCodec();
				BufferedReader r= IOUtil.openFileForBufferedReading(geneBed);
				geneList = r.lines().
						map(L->codec.decode(L)).
						filter(B->B!=null).
						map(B-> new Interval(B.getContig(),B.getStart(),B.getEnd(),true,B.get(3))).
						collect(Collectors.toList());
				r.close();
				}
			if(geneList.isEmpty())
				{
				LOG.error("gene List is empty");
				return -1;
				}
			final Comparator<VariantContext> ctxSorter = VCFUtils.createTidPosRefComparator(header.getSequenceDictionary());
			
			final Function<Interval,List<VariantContext>> loadVariants = (R)->{
				List<VariantContext> L=new ArrayList<>();
				CloseableIterator<VariantContext> r = this.vcfFileReader.query(R.getContig(), R.getStart(), R.getEnd());
				while(r.hasNext()) 
					{
					final VariantContext ctx= r.next();
					if(!accept(ctx)) continue;
					if(!getGenes(ctx).contains(R.getName())) continue;
					L.add(ctx);
					}
				r.close();
				return L;
				};
			final SkatExecutor executor = this.skatFactory.build();
			Double bestSkat=null;
			LOG.info("number of genes : "+ geneList.size());
			final int list_end_index=(
					this.user_end_index<0?
					geneList.size():
					Math.min(geneList.size(),this.user_end_index)
					);
			for(int x=this.user_begin_index ; x <list_end_index;++x)
				{
				final Interval intervalx= geneList.get(x);
				final List<VariantContext> variantsx = loadVariants.apply(intervalx);
				if(variantsx.isEmpty()) continue;
				
				for(int y=x;y<geneList.size() /* pas list_end_index */;++y)
					{
					final Interval intervaly;
					final List<VariantContext> merge;
					if(y==x) { //we-re testing gene 1 only
						intervaly = intervalx;
						merge = variantsx;
						}
					else
						{
						intervaly = geneList.get(y);
						if(intervaly.intersects(intervalx)) continue;
						final List<VariantContext> variantsy = loadVariants.apply(intervaly);
						if(variantsy.isEmpty()) continue;
						merge = new MergedList<>(variantsx,variantsy);
						}
					LOG.info("testing : ["+x+"]"+intervalx+" ["+y+"]"+intervaly+" N:"+geneList.size()+" best: "+bestSkat);

					final Double skat = eval(executor,merge);
					if(skat==null) continue;
					if(bestSkat==null || skat.compareTo(bestSkat)<0) {
						bestSkat=skat;
						LOG.info("best "+bestSkat +" "+intervalx+" "+intervaly);
						
						if(
							this.outputFile.getName().endsWith(".vcf") ||
							this.outputFile.getName().endsWith(".vcf.gz")
							)
							{
							final VCFHeader header2= new VCFHeader(header);
							header2.addMetaDataLine(new VCFHeaderLine(VcfGeneEpistasis.class.getName(),intervalx.getName()+" "+intervaly.getName()+" "+bestSkat));
							final VariantContextWriter w = VCFUtils.createVariantContextWriter(outputFile);
							w.writeHeader(header2);
							merge.stream().sorted(ctxSorter).forEach(V->w.add(V));
							w.close();
							}
						else
							{
							final PrintWriter w = super.openFileOrStdoutAsPrintWriter(outputFile);
							w.println(String.valueOf(bestSkat)+"\t"+intervalx.getName()+"\t"+intervaly.getName());
							w.flush();
							w.close();
							}
						}
						
					}
				}
			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(this.vcfFileReader);
			}
		}
	 	
	
	public static void main(final String[] args)
		{
		new VcfGeneEpistasis().instanceMainWithExit(args);
		}
	}
