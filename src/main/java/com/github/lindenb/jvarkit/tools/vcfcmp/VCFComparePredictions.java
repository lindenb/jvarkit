/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser.SnpEffPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;


/**
BEGIN_DOC

## Example

### VEP

```bash
$  java -jar dist/vcfcmppred.jar  f1.vcf f2.vcf 
(...)
7	8566286	rs2139	A	VEP discordant SO:terms between f1.vcf and f2.vcf	[SO:0001619, SO:0001632]
```


in f1.vcf (VEP 75) CSQ contains:

* intron_variant
* downstream_gene_variant
* nc_transcript_variant

in f2.vcf (VEP 71) CSQ contains:

* intron_variant

### SNPEFF

```bash
$  java -jar dist/vcfcmppred.jar  f1.vcf f2.vcf 
(...)
 8	1394127	.	G	SNPEff discordant SO:terms between between f1.vcf and f2.vcf	[SO:0001630]
```

in f1.vcf  (snpEff_3_6) EFF contains:

* downstream_gene_variant
* intron_variant
* splice_region_variant

in f2.vcf  (snpEff_3_4) EFF contains:

* downstream_gene_variant
* intron_variant




END_DOC
 */
@Program(name="vcfcmppred",description="Compare predictions (SNPEff, VEP) for several VCFs")
public class VCFComparePredictions extends AbstractVCFCompareBase {
	private VCFComparePredictions() {
	}
	private final Logger LOG=Logger.build(VCFComparePredictions.class).make();

	
	static private class PredictionTuple
		{
		VepPredictionParser vepPredictionParser;
		SnpEffPredictionParser snpEffPredictionParser;
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
	public int doWork(List<String> args) {
	
		PrintWriter out=null;
		SortingCollection<LineAndFile> variants=null;
		try
			{
			if(args.isEmpty())
				{
				LOG.error("Illegal number of arguments");
				return -1;
				}
			out= super.openFileOrStdoutAsPrintWriter(super.outputFile);

			variants=SortingCollection.newInstance(
					LineAndFile.class,
					new AbstractVCFCompareBase.LineAndFileCodec(),
					new AbstractVCFCompareBase.LineAndFileComparator(),
					super.sortingCollectionArgs.getMaxRecordsInRam(),
					super.sortingCollectionArgs.getTmpPaths()
					);
			variants.setDestructiveIteration(true);
			
			
			for(final String filename:args)
				{
				LOG.info("Reading from "+filename);
				Input input=super.put(variants, filename);
				LOG.info("end reading "+input.filename);
				}
			List<PredictionTuple> predictionTuples=new ArrayList<PredictionTuple>(super.inputs.size());
			for(AbstractVCFCompareBase.Input input:this.inputs)
				{
				PredictionTuple predictionTuple=new PredictionTuple();
				predictionTuple.snpEffPredictionParser=new SnpEffPredictionParserFactory(input.codecAndHeader.header).get();
				predictionTuple.vepPredictionParser=new VepPredictionParserFactory(input.codecAndHeader.header).get();
				
				predictionTuples.add(predictionTuple);
				}


			
			List<AbstractVCFCompareBase.LineAndFile> row=new ArrayList<LineAndFile>(super.inputs.size());
			
			
			CloseableIterator<LineAndFile> iter=variants.iterator();
			final Comparator<LineAndFile> posCompare = (A,B)->A.getContigPosRef().compareTo(B.getContigPosRef());

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
			out.close();out=null;
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
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
