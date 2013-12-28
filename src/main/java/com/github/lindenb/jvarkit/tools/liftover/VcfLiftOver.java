package com.github.lindenb.jvarkit.tools.liftover;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.picard.liftover.LiftOver;
import net.sf.picard.util.Interval;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;

import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfLiftOver extends AbstractVCFFilter2
	{
	private LiftOver liftOver=null;
	private File failedFile=null;
	private SAMSequenceDictionary newDict=null;

	@Override
	public String getProgramDescription() {
		return "Lift-over a VCF file.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfLiftOver";
		}
	
	protected Allele revcomp(Allele a)
		{
		if(a.isNoCall()) return a;
		if(a.isSymbolic()) return a;
		String seq=a.getBaseString();
		StringBuilder sb=new StringBuilder(seq.length());
		for(int i=seq.length()-1;i>=0;--i)
			{
			sb.append(AcidNucleics.complement(seq.charAt(i)));
			}
		return Allele.create(sb.toString(), a.isReference());
		}

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		final String TAG="LIFTOVER";
		VariantContextWriter failed=null;
		
		VCFHeader header=in.getHeader();
		

		if(this.failedFile!=null)
			{
			VCFHeader header2=new VCFHeader(header);
			header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			failed=VCFUtils.createVariantContextWriter(failedFile);
			failed.writeHeader(header2);
			}
		
		VCFHeader header3;
		
		if(newDict==null)
			{
			header3=new VCFHeader(header);
			warning("##contig files should be changed.");
			}
		else
			{
			Set<VCFHeaderLine> hh=new LinkedHashSet<VCFHeaderLine>();
			for(VCFHeaderLine h:header.getMetaDataInInputOrder())
				{
				if(!h.getKey().equals(VCFConstants.CONTIG_HEADER_KEY))
					{
					System.err.println(h);
					hh.add(h);
					}
				}
			for(SAMSequenceRecord ssr:newDict.getSequences())
				{
				Map<String,String> mapping=new HashMap<String,String>();
				mapping.put("ID", ssr.getSequenceName());
				mapping.put("length",String.valueOf(ssr.getSequenceLength()));
				VCFContigHeaderLine h=new VCFContigHeaderLine(mapping,ssr.getSequenceIndex());
				hh.add(h);
				}
			header3=new VCFHeader(
					hh,header.getGenotypeSamples()
					);
			}
		
		
		header3.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header3.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header3.addMetaDataLine(new VCFInfoHeaderLine(TAG,1,VCFHeaderLineType.String,"Chromosome|Position before liftOver."));
		out.writeHeader(header3);
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			progress.watch(ctx.getChr(),ctx.getStart());
			
			Interval lifted=liftOver.liftOver(
					new Interval(ctx.getChr(),ctx.getStart(),ctx.getEnd(),
					false,//negative strand
					""));
			if(lifted==null )
				{
				if(failed!=null) failed.add(ctx);
				}
			else
				{
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.chr(lifted.getSequence());
				vcb.start(lifted.getStart());
				vcb.stop(lifted.getEnd());
				vcb.attribute(TAG,ctx.getChr()+"|"+ctx.getStart() );
				
				if(lifted.isNegativeStrand())//strandess has changed
					{
					List<Genotype> genotypes=new ArrayList<Genotype>(header.getSampleNamesInOrder().size());
					Set<Allele> alleles=new HashSet<Allele>();
					alleles.add(revcomp(ctx.getReference()));
					for(String sample:header.getSampleNamesInOrder())
						{
						Genotype g=ctx.getGenotype(sample);
						GenotypeBuilder gb=new GenotypeBuilder(g);
						List<Allele> alleles2=new ArrayList<Allele>();
						for(Allele a0:g.getAlleles())
							{
							alleles2.add(revcomp(a0));
							}
						alleles.addAll(alleles2);
						gb.alleles(alleles2);
						genotypes.add(gb.make());
						}
					vcb.genotypes(genotypes);
					vcb.alleles(alleles);
					}
				out.add(vcb.make());
				}
			}
		CloserUtil.close(failed);	
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -f (chain-file) LiftOver file. Required.");
		out.println(" -m (double) lift over min-match. default:"+LiftOver.DEFAULT_LIFTOVER_MINMATCH);
		out.println(" -X (file.vcf) write variants failing the liftOver here. Optional.");
		out.println(" -D (reference) indexed reference file with the new sequence dictionary. Optional.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		double minMatch=LiftOver.DEFAULT_LIFTOVER_MINMATCH;
		File liftOverFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "f:m:X:D:"))!=-1)
			{
			switch(c)
				{
				case 'D': 
					{
					try
						{
						this.newDict=new SAMSequenceDictionaryFactory().load(new File(opt.getOptArg()));
						}
					catch(Exception err)
						{
						error(err);
						return -1;
						}
					break;
					}
				case 'X': this.failedFile=new File(opt.getOptArg()); break;
				case 'f': liftOverFile=new File(opt.getOptArg()); break;
				case 'm': minMatch=Double.parseDouble(opt.getOptArg()); break;
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
		if(liftOverFile==null)
			{
			error("LiftOver file is undefined.");
			return -1;
			}
		this.liftOver=new LiftOver(liftOverFile);
		this.liftOver.setLiftOverMinMatch(minMatch);
		return doWork(opt.getOptInd(), args);
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfLiftOver().instanceMainWithExit(args);
		}

	}
