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
public class VCFComparePredictions extends AbstractVCFCompareBase {
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
		out.print(ctx.getContig());
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
	private static void printSOSet(PrintWriter out,
		Set<SequenceOntologyTree.Term> set)
		{
		for(SequenceOntologyTree.Term t:set)
			{
			out.print(" ");
			out.print(t.getAcn()+"["+t.getLabel()+"]");
			}
		}
	private static void printDiscordantSO(
		PrintWriter out,
		Input input1,
		Set<SequenceOntologyTree.Term> so1,
		Input input2,
		Set<SequenceOntologyTree.Term> so2
		)
		{
		out.print("\t");
		Set<SequenceOntologyTree.Term> set=new HashSet<>(so1);
		set.removeAll(so2);
		if(!set.isEmpty())
			{
			out.print("only in "+input1.filename+":");
			printSOSet(out,set);
			}
		out.print("\t");
		set=new HashSet<>(so2);
		set.removeAll(so1);
		if(!set.isEmpty())
			{
			out.print("only in "+input2.filename+":");
			printSOSet(out,set);
			}
		out.println();
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
			
			final AbstractVCFCompareBase.LineAndFileComparator posCompare=new AbstractVCFCompareBase.LineAndFileComparator();

			factory.setComponentType(AbstractVCFCompareBase.LineAndFile.class);
			factory.setComparator(posCompare);
			factory.setTmpDirs(this.getTmpDirectories());
			factory.setCodec(new AbstractVCFCompareBase.LineAndFileCodec());
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
			for(AbstractVCFCompareBase.Input input:this.inputs)
				{
				PredictionTuple predictionTuple=new PredictionTuple();
				predictionTuple.snpEffPredictionParser=new SnpEffPredictionParser(input.codecAndHeader.header);
				predictionTuple.vepPredictionParser=new VepPredictionParser(input.codecAndHeader.header);
				
				predictionTuples.add(predictionTuple);
				}


			
			List<AbstractVCFCompareBase.LineAndFile> row=new ArrayList<LineAndFile>(super.inputs.size());
			
			
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
							Input input1=this.inputs.get(row.get(i).fileIdx);
							VariantContext ctx1=row.get(i).getContext();
							PredictionTuple predtuple1=predictionTuples.get(row.get(i).fileIdx);
							List<VepPrediction> vepPredictions1= predtuple1.vepPredictionParser.getPredictions(ctx1);
							List<SnpEffPrediction> snpEffPredictions1= predtuple1.snpEffPredictionParser.getPredictions(ctx1);
							Set<SequenceOntologyTree.Term> so_vep_1= getVepSoTerms(predtuple1.vepPredictionParser,ctx1);
							Set<SequenceOntologyTree.Term> so_snpeff_1= getSnpEffSoTerms(predtuple1.snpEffPredictionParser,ctx1);
							for(int j=i+1;j< row.size();++j)
								{
								Input input2=this.inputs.get(row.get(j).fileIdx);
								VariantContext ctx2=row.get(j).getContext();
								PredictionTuple predtuple2=predictionTuples.get(row.get(j).fileIdx);
								List<VepPrediction> vepPredictions2= predtuple2.vepPredictionParser.getPredictions(ctx2);
								List<SnpEffPrediction> snpEffPredictions2= predtuple2.snpEffPredictionParser.getPredictions(ctx2);

								Set<SequenceOntologyTree.Term> so_vep_2= getVepSoTerms(predtuple2.vepPredictionParser,ctx2);
								Set<SequenceOntologyTree.Term> so_snpeff_2= getSnpEffSoTerms(predtuple2.snpEffPredictionParser,ctx2);
								
								if(vepPredictions1.size()!=vepPredictions2.size())
									{
									startLine(out,ctx);
									out.print("\tVEP discordant transcripts count");
									out.print("\t"+input1.filename+":"+vepPredictions1.size());
									out.print("\t"+input2.filename+":"+vepPredictions2.size());
									out.println();
									printed=true;
									}
								if(snpEffPredictions1.size()!=snpEffPredictions2.size())
									{
									startLine(out,ctx);
									out.print("\tSNPEFF discordant transcripts count");
									out.print("\t"+input1.filename+":"+snpEffPredictions1.size());
									out.print("\t"+input2.filename+":"+snpEffPredictions2.size());
									out.println();
									printed=true;
									}

								
								if(!unshared(so_vep_1, so_vep_2).isEmpty())
									{
									startLine(out,ctx);
									out.print("\tVEP discordant SO:terms");
									printDiscordantSO(out,
											input1,
											so_vep_1,
											input2,
											so_vep_2
											);
									printed=true;
									}
								
								if(!unshared(so_snpeff_1, so_snpeff_2).isEmpty())
									{
									startLine(out,ctx);
									out.print("\tSNPEFF discordant SO:terms");
									printDiscordantSO(out,
										input1,
										so_snpeff_1,
										input2,
										so_snpeff_2
										);
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
