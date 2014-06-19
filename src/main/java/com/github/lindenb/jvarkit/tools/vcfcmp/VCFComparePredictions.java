/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser.SnpEffPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;


/**
 * @author lindenb
 *
 */
public class VCFComparePredictions extends AbstractVCFCompare {
	private VCFComparePredictions() {
	}
	static private class PredictionTuple
		{
		VepPredictionParser vepPredictionParser;
		SnpEffPredictionParser snpEffPredictionParser;
		}
	

	
	@Override
	public String getProgramDescription() {
		return "Compare predictions (SNPEff, VEP) for several VCFs";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VCFComparePredictions";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}
	
	private static void startLine(PrintWriter out,VariantContext ctx)
		{
		out.print(ctx.getChr());
		out.print("\t");
		out.print(ctx.getStart());
		out.print("\t");
		out.print(ctx.hasID()?ctx.getID():".");
		out.print("\t");
		out.print(ctx.getReference()==null?".":ctx.getReference().getDisplayString());
		}
	
	
	private static Set<SequenceOntologyTree.Term> unshared(
			Set<SequenceOntologyTree.Term> set1,
			Set<SequenceOntologyTree.Term> set2
			)
		{
		Set<SequenceOntologyTree.Term> so_terms=new HashSet<SequenceOntologyTree.Term>();
		for(SequenceOntologyTree.Term t:set1)
			{
			if(!set2.contains(t)) so_terms.add(t);
			}
		for(SequenceOntologyTree.Term t:set2)
			{
			if(!set1.contains(t)) so_terms.add(t);
			}
		return so_terms;
		}


	
	private static Set<SequenceOntologyTree.Term> getVepSoTerms(VepPredictionParser parser,VariantContext ctx)
		{
		Set<SequenceOntologyTree.Term> so_terms=new HashSet<SequenceOntologyTree.Term>();
		for(VepPrediction pred:parser.getPredictions(ctx))
			{
			so_terms.addAll(pred.getSOTerms());
			}
		return so_terms;
		}
	
	private static Set<SequenceOntologyTree.Term> getSnpEffSoTerms(SnpEffPredictionParser parser,VariantContext ctx)
		{
		Set<SequenceOntologyTree.Term> so_terms=new HashSet<SequenceOntologyTree.Term>();
		for(SnpEffPrediction pred:parser.getPredictions(ctx))
			{
			so_terms.addAll(pred.getSOTerms());
			}
		return so_terms;
		}

	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{	
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
		PrintWriter out=new PrintWriter(System.out);
		SortingCollection<LineAndFile> variants=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				error("Illegal number of arguments");
				return -1;
				}
			
			final AbstractVCFCompare.LineAndFileComparator posCompare=new AbstractVCFCompare.LineAndFileComparator();

			factory.setComponentType(AbstractVCFCompare.LineAndFile.class);
			factory.setComparator(posCompare);
			factory.setTmpDirs(this.getTmpDirectories());
			factory.setCodec(new AbstractVCFCompare.LineAndFileCodec());
			variants=this.factory.make();
			variants.setDestructiveIteration(true);
			
			
			for(int i=opt.getOptInd();i< args.length;++i)
				{
				String filename=args[i];
				info("Reading from "+filename);
				Input input=super.put(variants, filename);
				info("end reading "+input.filename);
				}
			List<PredictionTuple> predictionTuples=new ArrayList<PredictionTuple>(super.inputs.size());
			for(AbstractVCFCompare.Input input:this.inputs)
				{
				PredictionTuple predictionTuple=new PredictionTuple();
				predictionTuple.snpEffPredictionParser=new SnpEffPredictionParser(input.header);
				predictionTuple.vepPredictionParser=new VepPredictionParser(input.header);
				
				predictionTuples.add(predictionTuple);
				}


			
			List<AbstractVCFCompare.LineAndFile> row=new ArrayList<LineAndFile>(super.inputs.size());
			
			
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
						boolean printed=false;
						VariantContext ctx=row.get(0).getContext();
						if(row.size()!=this.inputs.size())
							{
							startLine(out,ctx);
							out.println("\tDiscordant number of variants");
							printed=true;
							}
						
						for(int i=0;i+1< row.size();++i)
							{
							VariantContext ctx1=row.get(i).getContext();
							PredictionTuple predtuple1=predictionTuples.get(row.get(i).fileIdx);
							Set<SequenceOntologyTree.Term> so_vep_1= getVepSoTerms(predtuple1.vepPredictionParser,ctx1);
							Set<SequenceOntologyTree.Term> so_snpeff_1= getSnpEffSoTerms(predtuple1.snpEffPredictionParser,ctx1);
							for(int j=i+1;j< row.size();++j)
								{
								VariantContext ctx2=row.get(j).getContext();
								PredictionTuple predtuple2=predictionTuples.get(row.get(j).fileIdx);
								Set<SequenceOntologyTree.Term> so_vep_2= getVepSoTerms(predtuple2.vepPredictionParser,ctx2);
								Set<SequenceOntologyTree.Term> so_snpeff_2= getSnpEffSoTerms(predtuple2.snpEffPredictionParser,ctx2);

								Set<SequenceOntologyTree.Term> discordant_so= unshared(so_vep_1, so_vep_2);
								if(!discordant_so.isEmpty())
									{
									startLine(out,ctx);
									out.print("\tVEP discordant SO:terms between ");
									out.print(this.inputs.get(row.get(i).fileIdx).filename);
									out.print(" and ");
									out.print(this.inputs.get(row.get(j).fileIdx).filename);
									out.print("\t");
									out.print(discordant_so);
									out.println();
									printed=true;
									}
								
								discordant_so= unshared(so_snpeff_1, so_snpeff_2);
								if(!discordant_so.isEmpty())
									{
									startLine(out,ctx);
									out.print("\tSNPEff discordant SO:terms between ");
									out.print(this.inputs.get(row.get(i).fileIdx).filename);
									out.print(" and ");
									out.print(this.inputs.get(row.get(j).fileIdx).filename);
									out.print("\t");
									out.print(discordant_so);
									out.println();
									printed=true;
									}
								
								}
														
							}
						if(!printed)
							{
							startLine(out,ctx);
							out.println("\tPASS");
							}
						
						row.clear();
						}
					if(rec==null) break;
					}
				
				row.add(rec);
				}
			iter.close();
			
			out.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
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
		new VCFComparePredictions().instanceMainWithExit(args);
		}
	}
