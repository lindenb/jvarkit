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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.liftover;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.ucsc.LiftOverChain;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Warning

you'd better have a look at gatk/picard LiftOverVcf


## Example

Lift over a hg19 file to hg18

Before liftover:

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" | grep -v "#

chr1	145273345	.	T	C	289.85	.	.	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
chr1	156011444	.	T	C	2523.46	.	.	GT:AD:DP:GQ:PL	0/1:24,15:40:99:214,0,443	0/1:32,36:68:99:702,0,794	1/1:1,59:61:99:1656,132,0	0/0:34,1:35:69.10:0,69,717
chr5	64982321	.	T	C	61.12	.	.	GT:AD:DP:GQ:PL	1/1:0,2:2:6:58,6,0	1/1:0,1:1:3.01:37,3,0	./.	./.
chr10	1142208	.	T	C	3404.30	.	.	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
chr10	126678092	.	G	A	89.08	.	.	GT:AD:DP:GQ:PL	0/0:64,3:67:99:0,165,1505	0/0:11,1:12:7.31:0,7,240	0/0:52,10:62:54.97:0,55,1263	0/1:35,9:44:99:125,0,693
chr10	135210791	.	T	C	65.41	.	.	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
chr13	48873835	.	G	A	58.95	.	.	GT:AD:DP:GQ:PL	./.	./.	1/1:0,2:2:6.01:62,6,0	1/1:0,1:1:3.01:31,3,0
chr20	36779424	.	G	A	128.76	.	.	GT:AD:DP:GQ:PL	0/0:49,1:52:63.68:0,64,969	0/0:17,0:17:30.05:0,30,320	0/0:93,0:94:99:0,216,2384	0/1:24,9:33:99:165,0,505
chrX	17819377	.	T	C	7515.25	.	.	GT:AD:DP:GQ:PL	1/1:0,125:126:99:2343,237,0	1/1:0,26:26:78.14:837,78,0	1/1:0,90:92:99:2640,244,0	1/1:0,74:75:99:1695,171,0
```

Running the liftover:

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  java -jar dist/jvarkit.jar vcfliftover \
      -f hg19ToHg18 \
      -X failing.vcf |\
grep -v "#"

chr1	143984702	.	T	C	289.85	.	LIFTOVER=chr1|145273345	GT:AD:DP:GQ:PL0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
chr1	154278068	.	T	C	2523.46	.	LIFTOVER=chr1|156011444	GT:AD:DP:GQ:PL0/1:24,15:40:99:214,0,443	0/1:32,36:68:99:702,0,794	1/1:1,59:61:99:1656,132,0	0/0:34,1:35:69.10:0,69,717
chr5	65018077	.	T	C	61.12	.	LIFTOVER=chr5|64982321	GT:AD:DP:GQ:PL1/1:0,2:2:6:58,6,0	1/1:0,1:1:3.01:37,3,0	./.	./.
chr10	1132208	.	T	C	3404.30	.	LIFTOVER=chr10|1142208	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
chr10	126668082	.	G	A	89.08	.	LIFTOVER=chr10|126678092	GT:AD:DP:GQ:PL	0/0:64,3:67:99:0,165,1505	0/0:11,1:12:7.31:0,7,240	0/0:52,10:62:54.97:0,55,1263	0/1:35,9:44:99:125,0,693
chr10	135060781	.	T	C	65.41	.	LIFTOVER=chr10|135210791	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
chr13	47771836	.	G	A	58.95	.	LIFTOVER=chr13|48873835	GT:AD:DP:GQ:PL./.	./.	1/1:0,2:2:6.01:62,6,0	1/1:0,1:1:3.01:31,3,0
chr20	36212838	.	G	A	128.76	.	LIFTOVER=chr20|36779424	GT:AD:DP:GQ:PL0/0:49,1:52:63.68:0,64,969	0/0:17,0:17:30.05:0,30,320	0/0:93,0:94:99:0,216,2384	0/1:24,9:33:99:165,0,505
chrX	17729298	.	T	C	7515.25	.	LIFTOVER=chrX|17819377	GT:AD:DP:GQ:PL1/1:0,125:126:99:2343,237,0	1/1:0,26:26:78.14:837,78,0	1/1:0,90:92:99:2640,244,0	1/1:0,74:75:99:1695,171,0
```




END_DOC
*/
@Program(
		name="vcfliftover",
		description="Lift-over a VCF file",
		keywords={"vcf","liftover"},
		modificationDate="20210603",
		creationDate="20240114",
		jvarkit_amalgamion = true
		)
public class VcfLiftOver extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfLiftOver.class);
	@Parameter(names={"-f","--chain"},description=LiftOverChain.OPT_DESC,required=true)
	private String liftOverFile = null;
	@Parameter(names={"-x","--failed"},description="(file.vcf) write variants failing the liftOver here. Optional.")
	private Path failedFile = null;
	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;
	@Parameter(names={"--adaptivematch"},description="Use adapative liftover minmatch using the ratio between the min allele size and the longest allele size")
	private boolean adaptivematch = false ;
	@Parameter(names={"-D","-R","-r","--reference"},description="Destination Reference. "+INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"-T","--tag"},description="INFO tag")
	private String infoTag = "LIFTOVER";
	@Parameter(names={"-failtag","--failtag"},description="failed INFO tag")
	private String failedinfoTag = "LIFTOVER_FAILED";
	@Parameter(names={"-check","--check"},description="Check variant allele sequence is the same on REF")
	private boolean checkAlleleSequence = false;
	@Parameter(names={"--ignore-chain-validation"},description="Ignore LiftOver chain validation")
	private boolean ignoreLiftOverValidation=false;
	@Parameter(names={"--indel","--indels"},description="do not LiftOver indels")
	private boolean ignoreIndels=false;
	@Parameter(names={"--info"},description="remove attribute from INFO on the fly")
	private Set<String> removeInfo=new HashSet<>();
	

	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter out) {
		VariantContextWriter failed=null;
		GenomicSequence genomicSequence = null;
		try (ReferenceSequenceFile indexedFastaSequenceFile =  ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx)) {
			final VCFHeader inputHeader= in.getHeader();
			
			final LiftOver liftOver;
			try {
				liftOver = LiftOverChain.load(
						this.liftOverFile,
						new SequenceDictionaryExtractor().extractRequiredDictionary(inputHeader),
						new SequenceDictionaryExtractor().extractRequiredDictionary(indexedFastaSequenceFile)
						);
				liftOver.setLiftOverMinMatch(this.userMinMatch);
				if(!this.ignoreLiftOverValidation) {
					liftOver.validateToSequences(indexedFastaSequenceFile.getSequenceDictionary());
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				return -1;
				}
			
			
			final Set<VCFHeaderLine> headerLines=inputHeader.getMetaDataInInputOrder().
					stream().filter(V->{
				if(!(V instanceof VCFInfoHeaderLine)) return true;
				final VCFInfoHeaderLine vih = VCFInfoHeaderLine.class.cast(V);
				if(removeInfo.contains(vih.getID())) return false;
				return true;
				}).collect(Collectors.toSet());
			
			
			if(this.failedFile!=null)
				{
				final VCFHeader header2=new VCFHeader(headerLines,inputHeader.getSampleNamesInOrder());
				header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
				header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
				header2.addMetaDataLine(new VCFInfoHeaderLine(this.failedinfoTag,1,VCFHeaderLineType.String,"Why the liftOver failed."));

				failed= super.writingVariantsDelegate.
						dictionary(header2).
						open(failedFile);
				failed.writeHeader(header2);
				}
			
			final VCFHeader header3=new VCFHeader(headerLines,inputHeader.getSampleNamesInOrder());
			header3.setSequenceDictionary(indexedFastaSequenceFile.getSequenceDictionary());
			JVarkitVersion.getInstance().addMetaData(this, header3);
			header3.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY, true));
			header3.addMetaDataLine(new VCFInfoHeaderLine(this.infoTag,1,VCFHeaderLineType.String,"Chromosome|Position before liftOver."));
			out.writeHeader(header3);
			while(in.hasNext())
				{
				VariantContext ctx= in.next();
				if(!this.removeInfo.isEmpty())
					{
					VariantContextBuilder vcb= new VariantContextBuilder(ctx);
					for(final String tag:this.removeInfo) vcb.rmAttribute(tag);
					ctx = vcb.make();
					}
				
				if(ctx.isIndel() && this.ignoreIndels)
					{
					if(failed!=null) failed.add(new VariantContextBuilder(ctx).attribute(this.failedinfoTag, "Indel").make());
					continue;
					}
				
				if(adaptivematch)
					{
					double minAlleleLength = Math.min(0,ctx.getAlleles().stream().mapToInt(A->A.length()).min().orElse(0));
					double maxAlleleLength =Math.max(1,ctx.getAlleles().stream().mapToInt(A->A.length()).max().orElse(1));
					liftOver.setLiftOverMinMatch(minAlleleLength /maxAlleleLength);
					}
				else
					{
					liftOver.setLiftOverMinMatch(this.userMinMatch);
					}
				
				final Interval lifted=liftOver.liftOver(
						new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd(),
						false,//negative strand
						String.join("|",ctx.getContig(),String.valueOf(ctx.getStart()),ctx.getReference().toString()))
						);
				if(lifted==null )
					{
					if(failed!=null) failed.add(new VariantContextBuilder(ctx).attribute(this.failedinfoTag, "LiftOverFailed").make());
					}
				else if(indexedFastaSequenceFile.getSequenceDictionary().getSequence(lifted.getContig())==null)
					{
					if(failed!=null) failed.add(new VariantContextBuilder(ctx).attribute(this.failedinfoTag, "ContigMissingDictionary|"+lifted.getContig()).make());
					}
				else
					{
					boolean alleleAreValidatedVsRef=true;
					//part of the code was copied from picard/liftovervcf
					final Map<Allele, Allele> reverseComplementAlleleMap = new HashMap<>();
					final List<Allele> alleles = new ArrayList<Allele>();

	                for (final Allele oldAllele : ctx.getAlleles()) {
	                	final Allele fixedAllele;
	                	if( oldAllele.isSymbolic() || oldAllele.isNoCall() || oldAllele.equals(Allele.SPAN_DEL))
	                		{
	                		alleles.add(oldAllele);
	                		continue;
	                		}
	                	else if (lifted.isPositiveStrand()) {
	                		fixedAllele = oldAllele;
	                        alleles.add(oldAllele);
	                    	}
	                    else {
	                        fixedAllele = Allele.create(SequenceUtil.reverseComplement(oldAllele.getBaseString()), oldAllele.isReference());
	                        alleles.add(fixedAllele);
	                        reverseComplementAlleleMap.put(oldAllele, fixedAllele);
	                    	}
	                    
                        if(this.checkAlleleSequence) {
                        	if(genomicSequence==null || !genomicSequence.getChrom().equals(lifted.getContig())) {
                        		genomicSequence = new GenomicSequence(indexedFastaSequenceFile, lifted.getContig());
                        		}
                        	final String alleleStr = fixedAllele.getBaseString();
                        	int x=0;
                        	while(x<alleleStr.length() && lifted.getStart()-1+x < genomicSequence.length())
                        		{
                        		final char refChar= genomicSequence.charAt(lifted.getStart()-1+x);
                        		if(Character.toLowerCase(refChar)!=Character.toLowerCase(alleleStr.charAt(x)))
                        			{
                        			alleleAreValidatedVsRef=false;
                        			break;
                        			}
                        		++x;
                        		}
                        	if(x!=alleleStr.length())
                        		{
                        		alleleAreValidatedVsRef=false;
                        		break;
                        		}
                        	}
	                	}
	                
	                if(!alleleAreValidatedVsRef)
	                	{
	                	if(failed!=null) failed.add(new VariantContextBuilder(ctx).attribute(this.failedinfoTag, "AlleleMismatchRef").make());
	                	continue;
	                	}
	                
	                if( lifted.getEnd() - lifted.getStart() != ctx.getEnd() - ctx.getStart())
	                	{
	                	if(failed!=null) failed.add(new VariantContextBuilder(ctx).attribute(this.failedinfoTag, "AlleleBadLength|"+lifted.length()).make());
	                	continue;
	                	}
				
					final VariantContextBuilder vcb=new VariantContextBuilder(
							ctx.getSource(),
							lifted.getContig(),
							lifted.getStart(),
							lifted.getEnd(),
							alleles
							);
					vcb.id(ctx.getID());
					vcb.attributes(ctx.getAttributes());
					vcb.attribute(this.infoTag,ctx.getContig()+"|"+ctx.getStart()+"|"+ctx.getReference().getDisplayString());
					vcb.filters(ctx.getFilters());
					vcb.log10PError(ctx.getLog10PError());
					if(lifted.getStart()!=lifted.getEnd()) {
						vcb.attribute(VCFConstants.END_KEY,lifted.getEnd());
					}
					  
					final GenotypesContext genotypeContext = ctx.getGenotypes();
					final GenotypesContext fixedGenotypes = GenotypesContext.create(genotypeContext.size());
			        for ( final Genotype genotype : genotypeContext ) 
			        	{
			            final List<Allele> fixedAlleles = new ArrayList<Allele>();
			            for ( final Allele allele : genotype.getAlleles() ) {
			                final Allele fixedAllele = reverseComplementAlleleMap.containsKey(allele) ?
			                		reverseComplementAlleleMap.get(allele) : 
			                		allele;
			                fixedAlleles.add(fixedAllele);
			            	}
			            fixedGenotypes.add(new GenotypeBuilder(genotype).alleles(fixedAlleles).make());
			        	}
			        vcb.genotypes(fixedGenotypes);
				    out.add(vcb.make());
					}
				}
			if(failed!=null)
				{
				failed.close();
				failed=null;
				}
			LOG.warning("you'd better have a look at gatk/picard LiftOverVcf");
			return 0;
		}catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(failed!=null) {
				try {failed.close();}
				catch(Throwable err)  {}
				}
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeVcf() {
		if(StringUtils.isBlank(this.liftOverFile))
			{
			LOG.error("LiftOver file is undefined.");
			return -1;
			}
		return super.beforeVcf();
		}
	
	@Override
	protected void afterVcf() {
		super.afterVcf();
		}

	public static void main(String[] args)
		{
		new VcfLiftOver().instanceMainWithExit(args);
		}

	}
