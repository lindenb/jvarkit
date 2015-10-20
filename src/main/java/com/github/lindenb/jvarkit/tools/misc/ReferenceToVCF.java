/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class ReferenceToVCF extends AbstractReferenceToVCF
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(ReferenceToVCF.class);

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractReferenceToVCF.AbstractReferenceToVCFCommand
	 	{

		private IntervalTreeMap<Boolean> limitBed=null;

	@Override
	public Collection<Throwable> call() throws Exception
		{
		final List<String> args = getInputFiles();
		Random random=new Random(0L);
		VariantContextWriter out=null;
		IndexedFastaSequenceFile fasta= null;
		
		try
			{
			if(args.size()!=1)
				{
				return wrapException(getMessageBundle("illegal.number.of.arguments"));
				}
			if(super.bedFile!=null)
				{
				BufferedReader r = null;
				this.limitBed=new IntervalTreeMap<Boolean>();
				try
					{
					LOG.info("reading "+super.bedFile);
					final Pattern tab=Pattern.compile("[\t]");
					r = IOUtils.openFileForBufferedReading(super.bedFile);
					String line;
					while((line=r.readLine())!=null)
						{
						if(line.startsWith("#") || line.isEmpty()) continue;
						String tokens[]=tab.split(line,4);
						this.limitBed.put(new Interval(
								tokens[0],
								1+Integer.parseInt(tokens[1]),
								1+Integer.parseInt(tokens[2])
								), true);
						}
					}
				catch(Exception err)
					{
					return wrapException(err);
					}
				finally
					{
					CloserUtil.close(r);
					}
				}
			fasta=new IndexedFastaSequenceFile(new File(args.get(0)));
			SAMSequenceDictionary dict=fasta.getSequenceDictionary();
			out= super.openVariantContextWriter();
			
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
				
				if(this.limitBed!=null)
					{
					Interval interval=new Interval(ssr.getSequenceName(),1,genome.length());
					if(!this.limitBed.containsOverlapping(interval)	) continue;
					}
				
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

					
					
					if(out.checkError()) break;
					
					
					}
				if(out.checkError()) break;
				}
			progress.finish();
			
			
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(fasta);
			}
		}
	 	}
	public static void main(String[] args) {
		new ReferenceToVCF().instanceMainWithExit(args);
	}
	}
