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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;



/**

BEGIN_DOC

### Example

```
$ find ~/src/gatk-ui/testdata/ -name "*.bam" > input.list

$ tail -2 input.vcf
rotavirus	1064	.	G	A	21.5606	.	DP=250;VDB=2.70971e-16;SGB=8.40135;RPB=0.935144;MQB=1;BQB=0.683886;MQ0F=0;AF1=0.25;G3=0.75,2.37734e-17,0.25;HWE=0.033921;AC1=2;DP4=0,219,0,31;MQ=60;FQ=22.8019;PV4=1,1.22605e-06,1,1	GT:PL	0/0:0,244,70	0/0:0,199,65	0/0:0,217,68	1/1:69,84,0
rotavirus	1064	.	G	A	21.5606	.	DP=250;VDB=2.70971e-16;SGB=8.40135;RPB=0.935144;MQB=1;BQB=0.683886;MQ0F=0;AF1=0.25;G3=0.75,2.37734e-17,0.25;HWE=0.033921;AC1=2;DP4=0,219,0,31;MQ=60;FQ=22.8019;PV4=1,1.22605e-06,1,1	GT:PL	./.	./.	./.	./.

$ java -jar dist/fixvcfmissinggenotypes.jar -d 50 --fixDP --filtered zz -B input.list input.vcf | tail -2
rotavirus	1064	.	G	A	21.56	.	AC1=2;AF1=0.25;BQB=0.683886;DP=188;DP4=0,219,0,31;FQ=22.8019;G3=0.75,2.37734e-17,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1.22605e-06,1,1;RPB=0.935144;SGB=8.40135;VDB=2.70971e-16	GT:DP:PL	0/0:48:0,244,70	0/0:63:0,199,65	0/0:53:0,217,68	1/1:24:69,84,0
rotavirus	1064	.	G	A	21.56	.	AC1=2;AF1=0.25;BQB=0.683886;DP=72;DP4=0,219,0,31;FQ=22.8019;G3=0.75,2.37734e-17,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1.22605e-06,1,1;RPB=0.935144;SGB=8.40135;VDB=2.70971e-16	GT:DP:FT:FXG	./.:48:PASS	0/0:63:zz:1	0/0:53:zz:1	./.:24:PASS
```

### Example


```
$ yourtool-mergingvcf 1.vcf 2.vcf 3.vcf > merged.vcf
$ find ./ -name "*.bam" > bams.list
$  java -jar dist/fixvcfmissinggenotypes.jar -f bams.list < merged.vcf > out.vcf
```

```
$ find DIR1 -name "PREFIX_*final.bam"  | grep -E '(S1|S2|S3|S4)' ) > bams.list

$ find DIR1 -name "PREFIX_*_variations.gatk.annotations.vcf.gz" |\
grep -E '(S1|S2|S3|S4)' |\
xargs perl  vcftools_0.1.12b/perl vcftools_0.1.12b/bin/vcf-merge |\
java -jar dist/fixvcfmissinggenotypes.jar -d 10 -f  bams.list |\
gzip --best > out.vcf.gz

```



### History

 * 2017-07-24 : rewrite whole program 
 * 2014: Creation


END_DOC
*/


@Program(name="fixvcfmissinggenotypes",
description="After a VCF-merge, read a VCF, look back at some BAMS to tells if the missing genotypes were homozygotes-ref or not-called. If the number of reads is greater than min.depth, then the missing genotypes is said hom-ref.",
biostars={119007,263309},
keywords={"sam","bam","vcf"}
)
public class FixVcfMissingGenotypes extends Launcher
	{
	private static final Logger LOG = Logger.build(FixVcfMissingGenotypes.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--depth"},description="minimal depth before setting a genotype to HOM_REF")
	private int minDepth = 10 ;
	@Parameter(names={"-B","--bams"},description="path of indexed BAM path with read Groups. You can put those paths in a text file having a *.list sufffix")
	private List<String> bamList=new ArrayList<>();
	@Parameter(names={"-filter","--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter filter  = SamFilterParser.buildDefault();
	@Parameter(names={"-T","--tag"},description="FORMAT 'Tag' for fixed genotype")
	private String fixedTag = "FXG";
	@Parameter(names={"--fixDP"},description="Update/create DP field even if genotype is called but there is no DP")
	private boolean fixDP=false;
	@Parameter(names={"--filtered"},description="Mark fixed genotypes as FILTERED with this FILTER")
	private String fixedGenotypesAreFiltered=null;

	
	/** return DP at given position */
	private int fetchDP(final VariantContext ctx,final String sample,List<SamReader> samReaders)
		{
		int depth=0;
		if(samReaders==null) samReaders=Collections.emptyList();
		for(final SamReader sr: samReaders)
			{
			final SAMRecordIterator iter=sr.query(ctx.getContig(), ctx.getStart(), ctx.getEnd(), false);
			while(iter.hasNext())
				{
				final SAMRecord rec=iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(filter.filterOut(rec)) continue;
				final SAMReadGroupRecord rg=rec.getReadGroup();
				if(!sample.equals(rg.getSample())) continue;
				final Cigar cigar=rec.getCigar();
				if(cigar==null) continue;
				int refPos=rec.getAlignmentStart();
				for(final CigarElement ce:cigar.getCigarElements())
					{
					if( refPos > ctx.getEnd() ) break;
					if(!ce.getOperator().consumesReferenceBases()) continue;
					if( ce.getOperator().consumesReadBases())
						{
						for(int n=0;n< ce.getLength();++n )
							{
							if( refPos+n < ctx.getStart() ) continue;
							if( refPos+n > ctx.getEnd()) break;
							depth++;
							}
						}
					refPos+= ce.getLength();
					}
				}
			iter.close();
			}
	depth /= ( 1 + ctx.getEnd() - ctx.getStart() );
	return depth;
	}
	
	@Override
	protected int doVcfToVcf(final String inputName,final VcfIterator in,final VariantContextWriter out) {
		final Set<String> bamFiles=  IOUtils.unrollFiles(bamList);
		final Map<String,List<SamReader>> sample2bam=new HashMap<>(bamFiles.size());

		try {
			final VCFHeader header=in.getHeader();
			for(final String bamFile: bamFiles)
				{
				LOG.info("Reading header for "+bamFile);
				final SamReader reader=super.openSamReader(bamFile);
				if(!reader.hasIndex())
					{
					LOG.error("No BAM index available for "+bamFile);
					return -1;
					}		
				final SAMFileHeader samHeader=reader.getFileHeader();
				for(final SAMReadGroupRecord g:samHeader.getReadGroups())
					{
					if(g.getSample()==null) continue;
					final String sample=g.getSample();
					if(StringUtil.isBlank(sample)) continue;
					List<SamReader> readers = sample2bam.get(sample);
					if(readers==null)
						{
						readers=new ArrayList<>();
						sample2bam.put(sample,readers);
						}
					readers.add(reader);
					}
				}
			
			
			final VCFHeader h2=new VCFHeader(header);
			if(h2.getFormatHeaderLine(VCFConstants.DEPTH_KEY)==null)
				{
				h2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
				}
			if(h2.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)==null)
				{
				h2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
				}
			if(h2.getFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY)==null)
				{
				h2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
				}
			
			h2.addMetaDataLine(new VCFFormatHeaderLine(this.fixedTag,1,VCFHeaderLineType.Integer,"Genotype was set as homozygous (min depth ="+this.minDepth+")"));
			
			
			SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
			
			out.writeHeader(h2);
			while(in.hasNext())
				{
				final VariantContext ctx = progress.watch(in.next());
				boolean somethingWasChanged=false;
				final List<Genotype> genotypes = new ArrayList<>(ctx.getNSamples());
				for(int i=0;i< ctx.getNSamples();++i)
					{
					Genotype genotype = ctx.getGenotype(i);
					final String sample = genotype.getSampleName();
					if(!genotype.isCalled() || (!genotype.hasDP() && this.fixDP))
						{
						List<SamReader> samReaders = sample2bam.get(sample);
						if(samReaders==null) samReaders=Collections.emptyList();
						
						final int depth=fetchDP(ctx,sample,samReaders);
						
						if(genotype.isCalled() && !genotype.hasDP())
							{
							genotype = new GenotypeBuilder(genotype).DP(depth).make();
							somethingWasChanged=true;
							}
						else // genotype was not called
							{
							if(depth>= this.minDepth)
								{	
								final List<Allele> homozygous=new ArrayList<>(2);
								homozygous.add(ctx.getReference());
								homozygous.add(ctx.getReference());
								final GenotypeBuilder gb=new GenotypeBuilder(genotype);
								gb.alleles(homozygous);
								gb.attribute(this.fixedTag, 1);
								gb.DP(depth);
								if(!StringUtil.isBlank(this.fixedGenotypesAreFiltered)) gb.filter(this.fixedGenotypesAreFiltered);
								somethingWasChanged=true;
								genotype = gb.make();
								}
							else if(!genotype.hasDP()) // cannot fix but we can update DP
								{
								genotype = new GenotypeBuilder(genotype).DP(depth).make();
								somethingWasChanged=true;
								}
							}
						}
					genotypes.add(genotype);
					}//end of for-each genotype
				
				if(somethingWasChanged)
					{
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.genotypes(genotypes);
					out.add(VCFUtils.recalculateAttributes(vcb.make()));
					}
				else
					{
					out.add(ctx);
					}
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
			sample2bam.values().stream().flatMap(L->L.stream()).forEach(R->CloserUtil.close(R));
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	public static void main(final String[] args) {
		new FixVcfMissingGenotypes().instanceMainWithExit(args);
	}

}
