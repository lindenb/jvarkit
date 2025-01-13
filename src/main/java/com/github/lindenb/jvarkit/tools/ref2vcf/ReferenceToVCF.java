/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.ref2vcf;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
/**

## Example

```bash
$   java -jar dist/jvarkit.jar fasta2vcf -L regions.bed -i 1 -d 1  human_g1k_v37.fasta |\
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
@Program(name="fasta2vcf",
	description="Creates a VCF containing all the possible substitutions from a Reference Genome.",
	keywords={"vcf","reference","fasta"},
	modificationDate="20240711",
	creationDate="20140910",
	jvarkit_amalgamion = true
	)
public class ReferenceToVCF extends Launcher
	{
	private static final Logger LOG = Logger.build(ReferenceToVCF.class).make();
	

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-L","--bed"},description="limit to this BED")
	private Path  bedFile = null;
	@Parameter(names={"-i","--insertions"},description="generate insertions")
	private int  insertion_size = 0;
	@Parameter(names={"-d","--deletions"},description="generate deletions")
	private int  deletion_size = 0;
	@Parameter(names={"-A","--disjoint"},description="disjoint ALT. Write one variant per ALT for SNVs")
	private boolean  disjoint_alts=false;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();
	
	
	
	@Override
	public int doWork(final List<String> args) {
		IntervalTreeMap<Boolean> limitBed=null;

		if(this.bedFile!=null)
			{
			if(limitBed==null) limitBed=new IntervalTreeMap<Boolean>();
			try
				{
				try(final BedLineReader r = new BedLineReader(this.bedFile)) {
					while(r.hasNext())
						{
						final BedLine record = r.next();
						limitBed.put(record.toInterval(), true);
						}
					}
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			}
		final Random random=new Random(0L);
		try
			{
			final ReferenceSequenceFile fasta= ReferenceSequenceFileFactory.getReferenceSequenceFile(Paths.get(oneAndOnlyOneFile(args)));
			final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(fasta);
			try(VariantContextWriter out= this.writingVariants.dictionary(dict).open(this.outputFile)) {
				final VCFHeader header=new VCFHeader();
				JVarkitVersion.getInstance().addMetaData(this, header);
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
					final GenomicSequence genome=new GenomicSequence(fasta, ssr.getSequenceName());
					
					if(limitBed!=null)
						{
						if(!limitBed.containsOverlapping(ssr)	) continue;
						}
					
					for(int n=0;n< genome.length();++n)
						{
						
						if(limitBed!=null && !limitBed.containsOverlapping(new SimplePosition(ssr.getSequenceName(), n+1))) {
							continue;
							}
						
						List<Allele> alleles=null;
						final byte ref=(byte)genome.charAt(n);
						switch(ref)
							{
							case 'a': case 'A':alleles = combination.get(0);break;
							case 'c': case 'C':alleles = combination.get(1);break;
							case 'g': case 'G':alleles = combination.get(2);break;
							case 't': case 'T':alleles = combination.get(3);break;
							default:break;
							}
						if(alleles==null) continue;
						
						
						if(!disjoint_alts)
							{
							final VariantContextBuilder vcb=new VariantContextBuilder();
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
								final VariantContextBuilder vcb=new VariantContextBuilder();
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
							
							final StringBuilder sb=new StringBuilder(insertion_size+2);
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
							
							final VariantContextBuilder vcb=new VariantContextBuilder();
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
						}
					}
				}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
	
		}
	
	public static void main(final String[] args) {
		new ReferenceToVCF().instanceMainWithExit(args);
	}
	}
