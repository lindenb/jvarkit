package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.regex.Pattern;

import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTree;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfCadd extends AbstractVCFFilter2
	{
	private TabixFileReader tabix=null;
	private int buffer_size=100000;
	private String ccaduri="http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz";
	
	private static class Record
		{
		int pos;
		Allele ref;
		Allele alt;
		float score;
		float phred;
		}
	
	private VcfCadd()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "Annotate VCF with  Combined Annotation Dependent Depletion (CADD) (Kircher & al. "+
				"A general framework for estimating the relative pathogenicity of human genetic variants. "+
				"Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892." +
				"PubMed PMID: 24487276.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfCadd";
		}
	
	

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header=in.getHeader();
		if(header.getSequenceDictionary()!=null)
			{
			SAMSequenceDictionary dict=header.getSequenceDictionary();
			Set<String> vcfchr=new HashSet<String>();
			for(SAMSequenceRecord ssr:dict.getSequences()) vcfchr.add(ssr.getSequenceName());
			if(!vcfchr.retainAll(this.tabix.getChromosomes()))//nothing changed
				{
				warning("#### !!!! NO common chromosomes between tabix and vcf file. Check chromosome 'chr' prefix ?");
				}
			}
		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		final String CCAD_CSCORE="CADD_CScore";
		final String CCAD_PHRED="CADD_PHRED";
		header.addMetaDataLine(new VCFInfoHeaderLine(CCAD_CSCORE,1,VCFHeaderLineType.Float,
				"CADD CScore. Suggests that that variant is likely to be  observed (negative values) vs simulated(positive values)."+
				"However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more likely to have deleterious effects"
						));
		header.addMetaDataLine(new VCFInfoHeaderLine(CCAD_PHRED,1,VCFHeaderLineType.Float,
				"CADD PHRED. Expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc"
				));
		
		Pattern tab=Pattern.compile("[\t]");
		out.writeHeader(header);
		String prevChromBuffer=null;
		long prevEndBuffer=-1;
		IntervalTree<Record> buffer=null;
		while(in.hasNext() )
			{	
			VariantContext ctx=in.next();
			
			progress.watch(ctx.getChr(), ctx.getStart());
			
			if(	buffer==null || 
				prevChromBuffer==null ||
				!ctx.getChr().equals(prevChromBuffer) ||
				prevEndBuffer<=ctx.getStart())
				{
				prevChromBuffer=ctx.getChr();
				prevEndBuffer=ctx.getEnd()+this.buffer_size;
				info("Fill buffer "+ctx.getChr()+":"+ctx.getStart()+"-"+prevEndBuffer);
				buffer=new IntervalTree<Record>();
				for(Iterator<String> iter=tabix.iterator(ctx.getChr(),
					(int)Math.max(1,ctx.getStart()-1),
					(int)Math.min(Integer.MAX_VALUE,prevEndBuffer)
					);
					iter.hasNext();
					)
					{
					String line=iter.next();
					String tokens[]=tab.split(line);
					if(tokens.length!=6) throw new IOException("Bad CADD line . Expected 6 fields:"+line);
					Record rec=new Record();
					rec.pos= Integer.parseInt(tokens[1]);
					rec.ref=Allele.create(tokens[2],true);
					rec.alt=Allele.create(tokens[3],false);
					rec.score=Float.parseFloat(tokens[4]);
					rec.phred=Float.parseFloat(tokens[5]);
					buffer.put(rec.pos,rec.pos, rec);
					}
				info("Fill: Done.");
				}
			
			boolean found=false;
			for( Iterator<IntervalTree.Node<Record>> reciter= buffer.iterator(ctx.getStart(), ctx.getEnd());
				!found && reciter.hasNext() ;
				)
				{
				Record rec=reciter.next().getValue();
				if(rec.pos!=ctx.getStart()) continue;
				if(!ctx.getReference().equals(rec.ref)) continue;
				for(Allele alt:ctx.getAlternateAlleles())
					{
					if(!alt.equals(rec.alt)) continue;				
					VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.attribute(CCAD_CSCORE,rec.score);
					vcb.attribute(CCAD_PHRED,rec.phred);
					out.add(vcb.make());
					found=true;
					break;
					}
				}
			
			if(!found)
				{
				out.add(ctx);
				}
			}
		progress.finish();
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -u (uri) Combined Annotation Dependent Depletion (CADD) Tabix file URI . Default:"+this.ccaduri);
		out.println(" -d (int) buffer size . Default:"+this.buffer_size);
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "u:b:"))!=-1)
			{
			switch(c)
				{
				case 'b': this.buffer_size=Math.max(1,Integer.parseInt(opt.getOptArg()));break;
				case 'u': this.ccaduri=opt.getOptArg(); break;
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(!this.ccaduri.endsWith(".gz"))
			{
			error("CCAD uri should end with gz. got "+this.ccaduri);
			return -1;
			}
		
		try
			{
			info("Loading index for "+this.ccaduri+". Please wait...");
			
		
			this.tabix=new TabixFileReader(this.ccaduri);
			info("End loading index");
			return doWork(opt.getOptInd(), args);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(tabix);
			}
		}

	public static void main(String[] args)
		{
		new VcfCadd().instanceMainWithExit(args);
		}
	}
