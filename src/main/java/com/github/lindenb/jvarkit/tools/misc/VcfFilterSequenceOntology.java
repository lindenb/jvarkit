package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;


public class VcfFilterSequenceOntology extends AbstractVCFFilter2
	{
	private Set<SequenceOntologyTree.Term> user_terms=new HashSet<SequenceOntologyTree.Term>();
	private boolean inverse_result=false;
	private final SequenceOntologyTree sequenceOntologyTree=SequenceOntologyTree.getInstance();
	private VcfFilterSequenceOntology()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "Filter a VCF file annotated with SNPEff or VEP with terms from Sequence-Ontology." +
				" Reasoning : Children of user's SO-terms will be also used.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfFilterSequenceOntology";
		}
	
	private boolean hasUserTem(Set<SequenceOntologyTree.Term> ctxTerms)
		{
		for(SequenceOntologyTree.Term ctxTerm:ctxTerms)
			{
			if(this.user_terms.contains(ctxTerm))
				{
				info("CONTAINS: "+ctxTerm.getAcn()+" "+ctxTerm.getLabel());
				return true;
				}
			}
		return false;
		}

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header=in.getHeader();
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		out.writeHeader(header);
		
		VepPredictionParser vepParser=new VepPredictionParser(header);
		SnpEffPredictionParser snpEffparser=new SnpEffPredictionParser(header);
		
		while(in.hasNext() )
			{	
			VariantContext ctx=in.next();
			boolean keep=false;
			
			for(SnpEffPredictionParser.SnpEffPrediction pred:snpEffparser.getPredictions(ctx))
				{
				if(hasUserTem(pred.getSOTerms())) { keep=true; break;}
				}
			if(!keep)
				{
				for(VepPredictionParser.VepPrediction pred:vepParser.getPredictions(ctx))
					{
					if(hasUserTem(pred.getSOTerms())) { keep=true; break;}
					}
				}
			
			if(inverse_result ) keep=!keep;
			if(keep) out.add(ctx);
			}
		
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -A (SO:ACN). add this SO:ACN");
		out.println(" -f (filename). Tab delimited file of SO accession numbers");
		out.println(" -S list the available SO accession and exit.");
		out.println(" -v invert selection.");
		out.println(" -d disable reasoning, don't use term's children.");
		super.printOptions(out);
		}
	
	private Set<SequenceOntologyTree.Term> parseAccessions(File f) throws IOException
		{
		Set<SequenceOntologyTree.Term>  set=new HashSet<SequenceOntologyTree.Term>();
		BufferedReader in=IOUtils.openFileForBufferedReading(f);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.startsWith("#")) continue;
			line=line.trim();
			if(line.trim().isEmpty()) continue;
			SequenceOntologyTree.Term t=sequenceOntologyTree.getTermByAcn(line);
			if(t==null)
				{
				throw new IllegalArgumentException("Unknown SO:Accession \""+line+"\"");
				}
			set.add(t);
			}
		in.close();
		return set;
		}
	
	@Override
	public int doWork(String[] args)
		{
		boolean reasoning=true;
		Set<SequenceOntologyTree.Term> setInit=new HashSet<SequenceOntologyTree.Term>();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "f:A:Svd"))!=-1)
			{
			switch(c)
				{
				case 'd': reasoning=false; break;
				case 'v': inverse_result=true; break;
				case 'S': 
					{
					for(SequenceOntologyTree.Term t:sequenceOntologyTree.getTerms())
						{
						System.out.println(t.getAcn()+"\t"+t.getLabel());
						}
					return 0;
					}
				case 'A':
					{
					String acn=opt.getOptArg().trim();
					SequenceOntologyTree.Term t=sequenceOntologyTree.getTermByAcn(acn);
					if(t==null)
						{
						System.err.println("Unknown SO:Accession \""+acn+"\"");
						return -1;
						}
					setInit.add(t);
					break;
					}
				case 'f':
					{
					try
						{
						setInit.addAll(parseAccessions(new File(opt.getOptArg())));
						}
					catch(IOException err)
						{
						error(err);
						return -1;
						}
					break;
					}
				default: 
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		for(SequenceOntologyTree.Term t:setInit)
			{
			this.user_terms.add(t);
			if(reasoning) this.user_terms.addAll(t.getAllDescendants());
			}
		info("Will be using :"+this.user_terms.toString());
		return doWork(opt.getOptInd(), args);
		}

	public static void main(String[] args)
		{
		new VcfFilterSequenceOntology().instanceMainWithExit(args);
		}
	}
