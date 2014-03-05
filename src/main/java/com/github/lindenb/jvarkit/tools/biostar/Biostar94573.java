package com.github.lindenb.jvarkit.tools.biostar;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import net.sf.samtools.util.CloserUtil;

import org.broad.tribble.readers.LineIterator;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;
import org.broadinstitute.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public class Biostar94573 extends AbstractCommandLineProgram
	{
	private static final char CLIPPING=' ';
	private static final char DELETION='-';
	private static final char MISMATCH='*';
	private int align_length=0;
	
	private class Sequence
		{
		String name;
		StringBuilder seq=new StringBuilder();
		char at(int index)
			{
			return(index< 0 || index >=seq.length()?CLIPPING:Character.toUpperCase(seq.charAt(index)));
			}
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/Biostar94573";
		}
	@Override
	public String getProgramDescription() {
		return "Getting a VCF file from a CLUSTAW alignment. See also http://www.biostars.org/p/94573/";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-R (name) reference name. Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		String REF="chrUn";
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:"))!=-1)
			{
			switch(c)
				{
				case 'R': REF=opt.getOptArg();break;
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
		VariantContextWriter w=null;
		LineIterator r=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				r=IOUtils.openStdinForLineIterator();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				r=IOUtils.openURIForLineIterator(filename);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			Sequence consensus=new Sequence();
			consensus.name="";
		
			Map<String, Sequence> sample2seq=new HashMap<String, Biostar94573.Sequence>();
			int columnStart=-1;
			while(r.hasNext())
				{
				String line=r.next();
				
				if( line.trim().isEmpty() || line.startsWith("CLUSTAL W"))
					{
					columnStart=-1;
					continue;
					}
				Sequence curr;
				if(line.charAt(0)==' ')
					{
					if(columnStart==-1)
						{
						error("illegal consensus line for "+line);
						return -1;
						}	
					curr=consensus;
					}
				else
					{
					 if(columnStart==-1)
						 {
						columnStart=line.indexOf(' ');
						if(columnStart==-1)
							{
							error("no whithespace in "+line);
							return -1;
							}
						while(columnStart< line.length() && line.charAt(columnStart)==' ')
							{
							columnStart++;
							}
						}
					String seqname=line.substring(0, columnStart).trim();
					curr=sample2seq.get(seqname);
					if(curr==null)
						{
						curr=new Sequence();
						curr.name=seqname;
						sample2seq.put(curr.name, curr);
						}
					}
				curr.seq.append(line.substring(columnStart));
				this.align_length=Math.max(align_length, curr.seq.length());
				}
			CloserUtil.close(r);
			
			Set<VCFHeaderLine> vcfHeaderLines=new HashSet<VCFHeaderLine>();

			vcfHeaderLines.add(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth."));
			vcfHeaderLines.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
			vcfHeaderLines.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth"));
			vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			Map<String,String> mapping=new HashMap<String,String>();
			mapping.put("ID", REF);
			mapping.put("length",String.valueOf(this.align_length));
			vcfHeaderLines.add(new VCFContigHeaderLine(mapping,1));

			Set<String> samples=new TreeSet<String>(sample2seq.keySet());
			VCFHeader vcfHeader=new VCFHeader(vcfHeaderLines,samples);
			
			w=VCFUtils.createVariantContextWriterToStdout();
			w.writeHeader(vcfHeader);
			int pos1=0;
			while(pos1< align_length)
				{
				if(consensus.at(pos1)==MISMATCH) {++pos1;continue;}
				int pos2=pos1+1;
				while(pos2<align_length && consensus.at(pos2)!='*')
					{
					++pos2;
					}
				boolean is_subsitution=(pos1+1==pos2);
				if(is_subsitution && pos1!=0)//need pos1>0 because ALT contains prev base.
					{
					for(Sequence seq: sample2seq.values())
						{
						if(seq.at(pos1)=='-')
							{
							is_subsitution=false;
							break;
							}
 						}
					}
				
				Set<Allele> alleles=new HashSet<Allele>();
				
				VariantContextBuilder vcb=new VariantContextBuilder();
				List<Genotype> genotypes=new ArrayList<Genotype>(samples.size());
				
				String longest=null;
				Counter<String> countAlleles=new Counter<String>();
				Map<String,String> sample2genotype=new HashMap<String,String>(samples.size());
				for(String sample:samples)
					{
					Sequence seq=sample2seq.get(sample);
					String al=null;
					if(is_subsitution)
						{
						if(seq.at(pos1)==CLIPPING) continue;
						al=String.valueOf(seq.at(pos1));
						}
					else
						{
						StringBuilder sb=new StringBuilder(pos2-pos1);
						for(int i=pos1-1;//yes -1
								i<pos2;
								++i)
							{
							if(seq.at(i)==CLIPPING) continue;
							if(seq.at(i)==DELETION) continue;
							sb.append(seq.at(i));
							}
						if(sb.length()==0) continue;
						al=sb.toString();
						}
					if(longest==null || longest.length()< al.length())
						{
						longest=al;
						}
					countAlleles.incr(al);
					sample2genotype.put(sample,al);
					}
				
				if(countAlleles.isEmpty()) continue;
				String refAllStr=(is_subsitution?countAlleles.getMostFrequent():longest);
				Allele refAllele=Allele.create(refAllStr, true);
				alleles.add(refAllele);
				
				
				for(String sample:sample2genotype.keySet())
					{
					GenotypeBuilder gb=new GenotypeBuilder(sample);
					Allele al=null;
					if( sample2genotype.get(sample).equals(refAllStr))
						{
						al=refAllele;
						}
					else
						{
						al=Allele.create(sample2genotype.get(sample), false);
						alleles.add(al);
						}
					
					gb.alleles(Collections.singletonList(al));
					gb.DP(1);
					
					genotypes.add(gb.make());
					}
				int start=pos1+(is_subsitution?1:0);//got to 1-based ref if subst, for indel with use pos(base)-1
				vcb.start(start);
				vcb.stop(start+(refAllStr.length()-1));
				vcb.chr(REF);
				HashMap<String, Object> atts=new HashMap<String,Object>();
				atts.put(VCFConstants.DEPTH_KEY, genotypes.size());
				vcb.attributes(atts);
				vcb.alleles(alleles);
				vcb.genotypes(genotypes);
				w.add(vcb.make());
				pos1=pos2;
				}
			w.close();
			info("Done");
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(w);
			}
		}
	
	public static void main(String[] args) {
		new Biostar94573().instanceMainWithExit(args);
	}
}
