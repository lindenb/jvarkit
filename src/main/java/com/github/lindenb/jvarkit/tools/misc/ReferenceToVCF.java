package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public class ReferenceToVCF extends AbstractCommandLineProgram
	{
	private IntervalTreeMap<Boolean> limitBed=null;
	
	@Override
	public String getProgramDescription() {
		return "Creates a VCF containing all the possible substitutions from a Reference Genome.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-L (file) limit to this BED file.");
		out.println("-d (int) generate deletions (size).");
		out.println("-i (int) generate insertions (size).");
		out.println("-g group ALT for same REF.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		boolean group_alts=false;
		int insertion_size = 0;
		int deletion_size = 0;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"L:i:d:g"))!=-1)
			{
			switch(c)
				{
				case 'g': group_alts=true;break;
				case 'i': insertion_size = Integer.parseInt(opt.getOptArg()); break;
				case 'd': deletion_size = Integer.parseInt(opt.getOptArg()); break;
				case 'L':
					{
					if(limitBed==null) limitBed=new IntervalTreeMap<Boolean>();
					try
						{
						info("reading "+opt.getOptArg());
						Pattern tab=Pattern.compile("[\t]");
						LineIterator r=IOUtils.openURIForLineIterator(opt.getOptArg());
						while(r.hasNext())
							{
							String line=r.next();
							if(line.startsWith("#") || line.isEmpty()) continue;
							String tokens[]=tab.split(line,4);
							limitBed.put(new Interval(
									tokens[0],
									1+Integer.parseInt(tokens[1]),
									1+Integer.parseInt(tokens[2])
									), true);
							}
						CloserUtil.close(r);
						}
					catch(Exception err)
						{
						error(err);
						return -1;
						}
				
					break;
					}
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
		Random random=new Random(0L);
		VariantContextWriter out=null;
		try
			{
			if(opt.getOptInd()+1!=args.length)
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			IndexedFastaSequenceFile fasta=new IndexedFastaSequenceFile(new File(args[opt.getOptInd()]));
			SAMSequenceDictionary dict=fasta.getSequenceDictionary();
			out=VCFUtils.createVariantContextWriterToOutputStream(System.out);
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			VCFHeader header=new VCFHeader();
			header.setSequenceDictionary(dict);
			out.writeHeader(header);
			
			final List<List<Allele>> combination=new ArrayList<List<Allele>>(4); 
			//always keep REF as first allele please
			combination.add(Arrays.asList(Allele.create("A", true), Allele.create("C", false),Allele.create("G", false),Allele.create("T", false)));
			combination.add(Arrays.asList(Allele.create("C", true ),Allele.create("A", false),Allele.create("G", false),Allele.create("T", false)));
			combination.add(Arrays.asList(Allele.create("G", true ),Allele.create("A", false),Allele.create("C", false),Allele.create("T", false)));
			combination.add(Arrays.asList(Allele.create("T", true ),Allele.create("A", false),Allele.create("C", false),Allele.create("G", false)));
				
			
			for(SAMSequenceRecord ssr: dict.getSequences())
				{
				GenomicSequence genome=new GenomicSequence(fasta, ssr.getSequenceName());
				for(int n=0;n< genome.length();++n)
					{
					progress.watch(ssr.getSequenceIndex(), n);
					List<Allele> alleles=null;
					byte ref=(byte)genome.charAt(n);
					switch(ref)
						{
						case 'a': case 'A':alleles = combination.get(0);break;
						case 'c': case 'C':alleles = combination.get(1);break;
						case 'g': case 'G':alleles = combination.get(2);break;
						case 't': case 'T':alleles = combination.get(3);break;
						default:break;
						}
					if(alleles==null) continue;
					
					if(this.limitBed!=null)
						{
						Interval interval=new Interval(ssr.getSequenceName(), n+1, n+1);
						if(!this.limitBed.containsOverlapping(interval)	) continue;
						}
					if(group_alts)
						{
						VariantContextBuilder vcb=new VariantContextBuilder();
						vcb.chr(ssr.getSequenceName());
						vcb.start(n+1);
						vcb.stop(n+1);
						vcb.alleles(alleles);
						vcb.log10PError(-5);
						out.add(vcb.make());
						}
					else
						{
						for(int a=1;a< 4;++a)//index 0 is always REF
							{
							VariantContextBuilder vcb=new VariantContextBuilder();
							vcb.chr(ssr.getSequenceName());
							vcb.start(n+1);
							vcb.stop(n+1);
							vcb.alleles(Arrays.asList(alleles.get(0),alleles.get(a)));//index 0 is always REF
							vcb.log10PError(-5);
							out.add(vcb.make());
							}
						}
					
					if(insertion_size>0 &&
						n+1 < 	 genome.length() )
						{
						alleles=new ArrayList<Allele>(2);
						//REFERENCE
						alleles.add(Allele.create(""+genome.charAt(n)+genome.charAt(n+1),true));
						
						StringBuilder sb=new StringBuilder(insertion_size+2);
						sb.append(genome.charAt(n));
						for(int n2=0;n2<insertion_size;++n2)
							{	
							switch(random.nextInt(4))
								{
								case 0:sb.append('A');break; 
								case 1:sb.append('C');break; 
								case 2:sb.append('G');break; 
								case 3:sb.append('T');break; 
								}
							}
						sb.append(genome.charAt(n+1));
						alleles.add(Allele.create(sb.toString(),false));
						
						VariantContextBuilder vcb=new VariantContextBuilder();
						vcb.chr(ssr.getSequenceName());
						vcb.start(n+1);
						vcb.alleles(alleles);
						vcb.computeEndFromAlleles(alleles, n+1);
						vcb.log10PError(-5);
						out.add(vcb.make());
						}
					
					if(deletion_size>0 &&
							n+deletion_size+1 < 	 genome.length() )
						{
						
						alleles=new ArrayList<Allele>(2);
						
						//REF
						StringBuilder sb=new StringBuilder(deletion_size+2);
						sb.append(genome.charAt(n));
						int lastpos=n+1;
						for(int n2=0;n2<deletion_size;++n2,lastpos++)
							{	
							sb.append(genome.charAt(lastpos));
							
							}
						sb.append(genome.charAt(lastpos));
						alleles.add(Allele.create(sb.toString(),true));
						
						
						alleles.add(Allele.create(""+genome.charAt(n)+genome.charAt(lastpos),false));
						
						VariantContextBuilder vcb=new VariantContextBuilder();
						vcb.chr(ssr.getSequenceName());
						vcb.start(n+1);
						vcb.alleles(alleles);
						vcb.computeEndFromAlleles(alleles, n+1);
						vcb.log10PError(-5);
						out.add(vcb.make());
						}

					
					
					if(System.out.checkError()) break;
					
					
					}
				if(System.out.checkError()) break;
				}
			progress.finish();
			
			
			
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
			}
		}
	public static void main(String[] args) {
		new ReferenceToVCF().instanceMainWithExit(args);
	}
	}
