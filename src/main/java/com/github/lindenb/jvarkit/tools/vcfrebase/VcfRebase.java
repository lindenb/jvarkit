/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcfrebase;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.util.CloserUtil;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.Rebase;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfRebase extends AbstractVCFFilter2 {
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private Rebase rebase=Rebase.createDefaultRebase();
	private final float DEFAULT_WEIGHT=5f;

	private VcfRebase() {
		}
	

	@Override
	public String getProgramDescription() {
		return "Finds restriction sites overlaping variants";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfRebase";
		}
	
	
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		final String ATT="ENZ";
		GenomicSequence genomicSequence=null;
		VCFHeader header=in.getHeader();
		header.addMetaDataLine(new VCFInfoHeaderLine(ATT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Enzyme overlapping: Format: (Name,Site,Sequence,pos-1,strand)"));
		out.writeHeader(header);
		while(in.hasNext())
			{
			VariantContext var=in.next();
			if(genomicSequence==null || !genomicSequence.getChrom().equals(var.getChr()))
				{
				info("Loading sequence "+var.getChr());
				genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile,var.getChr());
				}
			
			Set<String> hits=new HashSet<String>();
			for(Rebase.Enzyme enz:this.rebase)
				{
				int start0=Math.max(0, var.getStart() - enz.size());
				for(int y=start0;y<=var.getStart();++y)
					{
					//run each strand
					for(int strand=0;strand<2;++strand)
						{
						int x=0;
						//loop over bases of the enzyme
						for(x=0;x< enz.size() && y+x < genomicSequence.length() ;++x )
							{
							char c=(strand==0?
									enz.at(x):
									AcidNucleics.complement(enz.at((enz.size()-1)-x))
									);
							if(!Rebase.compatible(genomicSequence.charAt(y+x),c)) break;
							}
						// match found
						if(x==enz.size())
							{
							StringBuilder b=new StringBuilder("(");
							b.append(enz.getName());
							b.append("|");
							b.append(enz.getDecl());
							b.append("|");
							for(x=0;x < enz.size();++x)
								{
								char c=genomicSequence.charAt(y+x);
								if(y+x>=var.getStart()-1 && y+x<=var.getEnd()-1)
									{
									c=Character.toLowerCase(c);
									}
								b.append(c);
								}
							b.append("|");
							b.append(y+1);
							b.append("|");
							b.append(strand==0?"+":"-");
							b.append(")");
							hits.add(b.toString());
							break;
							}
						if(enz.isPalindromic()) break;
						}
					}
				}
			if(hits.isEmpty())
				{
				out.add(var);
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(var);
			vcb.attribute(ATT, hits.toArray(new String[hits.size()]));
			out.add(vcb.make());
			}
		}
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of java.util.logging.Level . currently:"+getLogger().getLevel());
		out.println(" -E (name) restrict to that enzyme. Can be called multiple times. Optional.");
		out.println(" -R (fasta) path to reference sequence indexed with picard. Required.");
		out.println(" -w (float) min enzyme weight (e.g: A=1, N=0.0, S=0.5) default="+DEFAULT_WEIGHT+" . Optional.");
		}
	
	@Override
	public int doWork(String[] args)
		{
		float weight=DEFAULT_WEIGHT;
		File fasta=null;
		Set<String> onlyEnz=new HashSet<String>();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args, "hvL:E:R:w:"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(java.util.logging.Level.parse(opt.getOptArg()));break;
				case 'E': onlyEnz.add(opt.getOptArg()); break;
				case 'R': fasta=new File(opt.getOptArg()); break;
				case 'w': weight=Float.parseFloat(opt.getOptArg()); break;
				case ':': System.err.println("Missing argument for option -"+opt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+opt.getOptOpt());return -1;
				}
			}
		if(!onlyEnz.isEmpty())
			{
			Rebase rebase2=new Rebase();
			for(String e:onlyEnz)
				{
				Rebase.Enzyme enz=this.rebase.getEnzymeByName(e);
				if(enz==null)
					{
					System.err.println("Cannot find enzyme "+enz +" in RE list.");
					System.err.println("Current list is:");
					for(Rebase.Enzyme E: this.rebase)
						{
						System.err.println("\t"+E);
						}
					return -1;
					}
				rebase2.getEnzymes().add(enz);
				}
			this.rebase=rebase2;
			}
		
		int i=0;
		while(i< rebase.size())
			{
			if(rebase.get(i).getWeight()< weight)
				{
				rebase.getEnzymes().remove(i);
				}
			else
				{
				++i;
				}
			}
		
		if(rebase.size()==0)
			{
			warning("REBASE IS EMPTY");
			}
		
		if(fasta==null)
			{
			super.error("Undefined fasta sequence");
			return -1;
			}
		try
			{
			info("opening "+fasta);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(fasta);
			return doWork(opt.getOptInd(), args);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		}

	public static void main(String[] args)
		{
		new VcfRebase().instanceMainWithExit(args);
		}
	}
