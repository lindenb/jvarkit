/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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


*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
/**

## Example

```bash
$   java -jar dist/referencetovcf.jar -L regions.bed -i 1 -d 1  human_g1k_v37.fasta |\
grep -v "##" | cut -f 1-5

#CHROM	POS	ID	REF	ALT
1	10874	.	A	C
1	10874	.	A	G
1	10874	.	A	T
1	10874	.	AC	AGC
1	10874	.	ACA	AA
1	10875	.	C	A
1	10875	.	C	G
1	10875	.	C	T
1	10875	.	CA	CTA

```

 */
@Program(name="referencetovcf",
	description="Creates a VCF containing all the possible substitutions from a Reference Genome.",
	keywords={"vcf","reference","fasta"}
	)
public class ReferenceToVCF extends Launcher
	{
	private static final Logger LOG = Logger.build(ReferenceToVCF.class).make();
	

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-L","--bed"},description="limit to this BED")
	private File  bedFile = null;
	@Parameter(names={"-i","--insertions"},description="generate insertions")
	private int  insertion_size = 0;
	@Parameter(names={"-d","--deletions"},description="generate deletions")
	private int  deletion_size = 0;
	@Parameter(names={"-A","--disjoint"},description="disjoint ALT")
	private boolean  disjoint_alts=false;

	private IntervalTreeMap<Boolean> limitBed=null;
	
	
	@Override
	public int doWork(final List<String> args) {

		if(this.bedFile!=null)
			{
			if(limitBed==null) limitBed=new IntervalTreeMap<Boolean>();
			try
				{
				final BedLineCodec codec = new BedLineCodec();
				BufferedReader r=IOUtils.openFileForBufferedReading(this.bedFile);
				String line;
				while((line=r.readLine())!=null)
					{
					if(BedLine.isBedHeader(line)) continue;
					final BedLine record = codec.decode(line);
					limitBed.put(record.toInterval(), true);
					}
				CloserUtil.close(r);
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			}
		final Random random=new Random(0L);
		VariantContextWriter out=null;
		try
			{
			final IndexedFastaSequenceFile fasta=new IndexedFastaSequenceFile(new File(oneAndOnlyOneFile(args)));
			SAMSequenceDictionary dict=fasta.getSequenceDictionary();
			out= super.openVariantContextWriter(this.outputFile);
			
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
					if(!disjoint_alts)
						{
						VariantContextBuilder vcb=new VariantContextBuilder();
						vcb.chr(ssr.getSequenceName());
						vcb.start(n+1);
						vcb.stop(n+1);
						vcb.alleles(alleles);
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
						out.add(vcb.make());
						}
					
					if(deletion_size>0 &&
							n+deletion_size+1 < 	 genome.length() )
						{
						
						alleles=new ArrayList<Allele>(2);
						
						//REF
						final  StringBuilder sb=new StringBuilder(deletion_size+2);
						sb.append(genome.charAt(n));
						int lastpos=n+1;
						for(int n2=0;n2<deletion_size;++n2,lastpos++)
							{	
							sb.append(genome.charAt(lastpos));
							
							}
						sb.append(genome.charAt(lastpos));
						alleles.add(Allele.create(sb.toString(),true));
						
						
						alleles.add(Allele.create(""+genome.charAt(n)+genome.charAt(lastpos),false));
						
						final  VariantContextBuilder vcb=new VariantContextBuilder();
						vcb.chr(ssr.getSequenceName());
						vcb.start(n+1);
						vcb.alleles(alleles);
						vcb.computeEndFromAlleles(alleles, n+1);
						out.add(vcb.make());
						}

					if(out.checkError()) break;
					}
				if(out.checkError()) break;
				}
			progress.finish();
			
			
			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	
	public static void main(final String[] args) {
		new ReferenceToVCF().instanceMainWithExit(args);
	}
	}
