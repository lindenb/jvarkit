/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter.OnNotFound;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFStandardHeaderLines;


/**
BEGIN_DOC

 ## Building dgv bed:

```
$ wget -O - -q "http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2016-05-15.txt" | 
	grep -v '^variantaccession' | 
	awk -F $'\t' '{printf("%s\t%s\t%s\t%s\n",$2,$3,$4,$1);}' | 
	sort -t $'\t' -k1,1V -k2,2n > dgv.bed
```

END_DOC

 */
@Program(name="validatecnv",
	description="experimental CNV detection. Look depths before/after putative known CNV.",
	keywords= {"cnv","bam","sam","vcf"},
	generate_doc=false
	)
public class ValidateCnv extends Launcher
	{
	private static final Logger LOG = Logger.build(ValidateCnv.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-B","--bed"},description="Bed file of CNV.",required=true)
	private File bedFile=null;	
	@Parameter(names={"-x","--extend"},description="Search the boundaries in a region that is 'x'*(CNV-length)")
	private double extendFactor=0.1;	
	@Parameter(names={"-md","--min-dp"},description="At least one of the bounds must have a median-depth greater than this value.")
	private int min_depth = 20;
	@Parameter(names={"-E","-del","--del","--deletion"},description="Deletion Treshold. Which fraction of the median depth is considered as aa deletion. Must be <1.0" )
	private double deletion_treshold = 0.5;
	@Parameter(names={"-U","-dup","--dup","--duplication"},description="Duplication Treshold. Which fraction of the median depth is considered as a duplication. Must be >1.0" )
	private double duplication_treshold = 1.5;

	
	private class Input implements Closeable
		{
		SamReader samReader;
		SAMFileHeader header;
		SAMSequenceDictionary dict;
		ContigNameConverter ctgNameConverter;
		String sampleName;
		Input(final String uri) throws IOException {
			final SamInputResource sri = SamInputResource.of(uri);
			this.samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(sri);
			this.header  = samReader.getFileHeader();
			this.dict = this.header.getSequenceDictionary();
			this.ctgNameConverter = ContigNameConverter.fromOneDictionary(this.dict);
			this.ctgNameConverter.setOnNotFound(OnNotFound.SKIP);
			this.sampleName = this.header.getReadGroups().stream().map(R->R.getSample()).filter(S->!StringUtil.isBlank(S)).findFirst().orElse(uri);
			}
		@Override
		public void close() throws IOException {
			this.samReader.close();
			}
		public OptionalDouble medianCov(final Locatable interval)
			{
			final String contig = this.ctgNameConverter.apply(interval.getContig());
			if(StringUtil.isBlank(contig)) return OptionalDouble.empty();
			final int array_size = 1+(interval.getEnd()-interval.getStart());
			if(array_size<=0) {
				LOG.warn("skipping negative interval "+interval);
				return OptionalDouble.empty();
			}
			CloseableIterator<SAMRecord> iter = this.samReader.queryOverlapping(contig,interval.getStart(),interval.getEnd());
			final double array[]=new double[array_size];
			Arrays.fill(array, 0);
			while(iter.hasNext()) {
				final SAMRecord rec= iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getReadFailsVendorQualityCheckFlag()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				int ref=rec.getStart();
				for(final CigarElement ce:cigar) {
					final CigarOperator op =ce.getOperator();
					if(op.consumesReferenceBases())
						{
						if(op.consumesReadBases()) {
							for(int x=0;x < ce.getLength();++x)
								{
								final int p = ref+x - interval.getStart();
								if(p<0 || p>=array.length) continue;
								array[p]++;
								}
							}
						ref+=ce.getLength();
						}
					}
				}
			iter.close();
			final Median median = new Median();
			return OptionalDouble.of(median.evaluate(array));
			}
		
		}
	

	
	@Override
	public int doWork(final List<String> args) {		
		if(this.extendFactor<=0)
			{
			LOG.error("bad extend factor "+this.extendFactor);
			return -1;
			}
		if(this.deletion_treshold>=1.0)
			{
			LOG.error("bad deletion treshold . Must be <1.0 but got "+this.deletion_treshold);
			return -1;
			}
		if(this.duplication_treshold<=1.0)
			{
			LOG.error("bad dup treshold . Must be >1.0 but got "+this.duplication_treshold);
			return -1;
			}
		final List<Input> inputs = new ArrayList<>();
		VariantContextWriter out = null;
		BufferedReader bedIn = null;
		try
			{
			final BedLineCodec bedCodec = new BedLineCodec();
			bedIn = IOUtils.openFileForBufferedReading(this.bedFile);
			
			final List<String> urls;
			if(args.size()==1 && args.get(0).endsWith(".list")) {
				urls = IOUtil.slurpLines(new File(args.get(0)));
				}
			else
				{
				urls = args;
				}
			
			for(final String url:urls)
				{
				if(StringUtil.isBlank(url)) continue;
				inputs.add(new Input(url));
				}
			final Set<VCFHeaderLine> metadata = new HashSet<>();
			
			final VCFFormatHeaderLine leftMedianDepth = 
					new VCFFormatHeaderLine("LDP",1,VCFHeaderLineType.Integer,"Left median depth or -1");
			metadata.add(leftMedianDepth);
			final VCFFormatHeaderLine rightMedianDepth = 
					new VCFFormatHeaderLine("RDP",1,VCFHeaderLineType.Integer,"Right median depth or -1");
			metadata.add(rightMedianDepth);
			
			
			VCFStandardHeaderLines.addStandardFormatLines(metadata, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.GENOTYPE_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(metadata, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.END_KEY
					);

			
			final VCFHeader header = new VCFHeader(metadata,
					inputs.stream().map(S->S.sampleName).collect(Collectors.toList())
					);
			
			
			out =  super.openVariantContextWriter(this.outputFile);
			out.writeHeader(header);
			final Allele refAllele = Allele.create("N", true);
			final Allele delAllele = Allele.create("<DEL>", false);
			final Allele dupAllele = Allele.create("<DUP>", false);
			String line;
			
			while((line=bedIn.readLine())!=null)
				{
				final BedLine bedLine = bedCodec.decode(line);
				if(bedLine==null) continue;
				final int length = bedLine.getEnd()-bedLine.getStart();
				final int extend = 1+(int)(length/this.extendFactor);
				
				
				
				
				
				final List<Genotype> genotypes = new ArrayList<>(inputs.size());
				for(final Input input:inputs)
					{
					OptionalDouble mid = input.medianCov(bedLine);
					final OptionalDouble left = input.medianCov(new Interval(
							bedLine.getContig(),
							Math.max(0, bedLine.getStart()-extend),
							Math.max(0, bedLine.getStart()-1)
							));
					final OptionalDouble right = input.medianCov(new Interval(
							bedLine.getContig(),
							bedLine.getEnd()+1,
							bedLine.getEnd()+extend
							));
					
					
					final Genotype gt;
					
					if(mid.isPresent())
						{
						final double midV = mid.getAsDouble();
						final GenotypeBuilder gb = new GenotypeBuilder(input.sampleName);
						gb.DP((int)midV);
						gb.attribute(leftMedianDepth.getID(),left.isPresent()?(int)left.getAsDouble():-1);
						gb.attribute(rightMedianDepth.getID(),left.isPresent()?(int)right.getAsDouble():-1);
						
						if(left.isPresent() &&  right.isPresent())
							{
							double leftV = left.getAsDouble();
							double rightV = right.getAsDouble();
							if(leftV < this.min_depth || rightV<this.min_depth)
								{
								gt = GenotypeBuilder.createMissing(input.sampleName, 2);
								}
							else if(midV <= leftV*this.deletion_treshold &&
									midV <= rightV*this.deletion_treshold )
								{
								System.err.println(leftV+" / "+midV+" / "+rightV);
								gb.alleles(Arrays.asList(refAllele,delAllele));
								gt = gb.make();
								}
							else if(midV >= leftV*this.duplication_treshold &&
									midV>= rightV*this.duplication_treshold &&
									midV >= this.min_depth
									)
								{
								gb.alleles(Arrays.asList(refAllele,dupAllele));
								gt = gb.make();
								}
							else
								{
								gb.alleles(Arrays.asList(refAllele,refAllele));
								gt = gb.make();
								}
							}
						else
							{
							gt = GenotypeBuilder.createMissing(input.sampleName, 2);
							}
						
						}
					else
						{
						gt = GenotypeBuilder.createMissing(input.sampleName, 2);
						}
					
					genotypes.add(gt);
					}
				
				final Set<Allele> alleles = new HashSet<>();
				alleles.add(refAllele);
				alleles.addAll(genotypes.stream().flatMap(G->G.getAlleles().stream()).filter(A->!A.isNoCall()).collect(Collectors.toSet()));
				
				final VariantContextBuilder vcb = new VariantContextBuilder(
						this.bedFile.getPath(),
						bedLine.getContig(),
						bedLine.getStart(),
						bedLine.getEnd(),
						alleles
						);
				if(!StringUtil.isBlank(bedLine.get(3)))
					{
					vcb.id(bedLine.get(3));
					}
				vcb.attribute(VCFConstants.END_KEY, bedLine.getEnd());
				vcb.genotypes(genotypes);
				
				out.add(vcb.make());
				}
			
			bedIn.close();
			out.close();
			CloserUtil.close(inputs);

			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(bedIn);
			CloserUtil.close(inputs);
			CloserUtil.close(out);
			}
		}
	
	

	
	public static void main(final String[] args) {
		new ValidateCnv().instanceMainWithExit(args);
		}
	}
