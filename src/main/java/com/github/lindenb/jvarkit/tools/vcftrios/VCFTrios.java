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
package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## Using GATK VariantAnnotator ? 

GATK VariantAnnotator doesn't work if GQ is low or if there is no GQ.

## Example

a pedigree file:

```
$  cat pedigree.txt 

A	SAMPLE_P	0	0	0
A	SAMPLE_M	0	0	0
A	SAMPLE_E	SAMPLE_P	SAMPLE_M	0
```


find mendelian incompatibilities:

```
$  gunzip -c input.vcf.gz |\
   java -jar dist/vcftrio.jar -p pedigree.txt | grep -E '(#CHROM|MENDEL=SAMPLE_E)' |\
   verticalize 

(...)
>>> 23
$1	#CHROM	X
$2	POS	0573
$3	ID	rs358
$4	REF	G
$5	ALT	A
$6	QUAL	85.60
$7	FILTER	PASS
$8	INFO	MENDEL=SAMPLE_E
$9	FORMAT	GT:DP:DP4:GP:GQ:PL
$10	SAMPLE_E	0/1:11:6,0,5,0:97,0,122:97:96,0,118
$11	SAMPLE_M	1/1:5:0,0,5,0:134,19,0:19:120,15,0
$12	SAMPLE_P	1/1:6:0,0,6,0:136,22,0:22:121,18,0
<<< 23
(...)
>>> 59
$1	#CHROM	Y
$2	POS	19
$3	ID	rs5678
$4	REF	CA
$5	ALT	C,CAA
$6	QUAL	31.86
$7	FILTER	PASS
$8	INFO	MENDEL=SAMPLE_E
$9	FORMAT	GT:DP:DP4:GP:GQ
$10	SAMPLE_E	2/2:80:3,0,43,34:.,.,108,.,203,0:99
$11	SAMPLE_M	.
$12	SAMPLE_P	1/1:53:0,0,27,26:81,99,0,.,.,.:81
<<< 59

```
## History

  * [20180907] moved to the new DeNovoDetector
  * [20180704] changing the arguments that are not really clear.

 
END_DOC

 */
@Program(
		name="vcftrio",
		description="Find mendelian incompatibilitie / denovo variants in a VCF",
		keywords={"vcf","mendelian","pedigree","denovo"}
		)
public class VCFTrios
	extends Launcher
	{
	private static final  Logger LOG = Logger.build(VCFTrios.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-p","--ped","--pedigree"},description="Pedigree file. "+Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile = null;
	@Parameter(names={"-fo","--filter-out","--filter-no-denovo"},description="FILTER name if there is NO mendelian violation.")
	private String filterNoIncompat = null;
	@Parameter(names={"-fi","--filter-in"},description="FILTER name if there is ANY mendelian violation.")
	private String filterAnyIncompat = null;
	@Parameter(names={"-gtf","--gtfilter"},description="GENOTYPE FILTER name. Create a filter in the GENOTYPE column when there is NO mendelian violation")
	private String genotypeFilterNameNoIncompat = null;
	@Parameter(names={"-A","--attribute"},description="INFO Attribute name containing the name of the affected samples.")
	private String attributeName = "MENDEL";
	@Parameter(names={"-d","--dicard"},description="Discard the variant if there is NO mendelian violation.")
	private boolean discard_variants_without_mendelian_incompat=false;	
	@Parameter(names={"-hr","--hom-ref"},description="[20180705] treat NO_CALL genotypes as HOM_REF (when individual VCF/Sample have been merged).")
	private boolean nocall_to_homref = false;
		
	private static class Trio
		{
		int father_id = -1;
		int mother_id = -1;
		int child_id = -1;
		}
	
	@Override
	public int doVcfToVcf(final String inputName, VCFIterator r, final VariantContextWriter w) {
		long count_incompats=0L;
		final Set<String> sampleNotFoundInVcf = new HashSet<>();
		Pedigree pedigree=null;
		final List<Trio> trios = new ArrayList<>();
		
		try		{
				final DeNovoDetector detector = new DeNovoDetector();
				detector.setConvertingNoCallToHomRef(this.nocall_to_homref);
				
				final VCFHeader header = r.getHeader(); 
				if(this.pedigreeFile!=null) {
					LOG.info("reading pedigree "+this.pedigreeFile);
					pedigree = Pedigree.newParser().parse(this.pedigreeFile);
					}
				else
					{
					pedigree = Pedigree.newParser().parse(header);
					}
				final VCFHeader h2=new VCFHeader(header);
				final Set<VCFHeaderLine> meta = new HashSet<>();
				meta.add(new VCFInfoHeaderLine(
						this.attributeName,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Samples with mendelian incompatibilities." +
							(this.pedigreeFile==null?"":" Pedigree File was : " + this.pedigreeFile)
						));
				meta.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY, true));
				
				if(!StringUtil.isBlank(this.filterAnyIncompat)) {
					meta.add(new VCFFilterHeaderLine(this.filterAnyIncompat,
							"Variant contains at least one mendelian incompatibilities"));
					}
				if(!StringUtil.isBlank(this.filterNoIncompat)) {
					meta.add(new VCFFilterHeaderLine(this.filterNoIncompat,
							"Variant does not contain any mendelian incompatibilities"));
					}
				
				meta.stream().forEach(H->h2.addMetaDataLine(H));
				
				JVarkitVersion.getInstance().addMetaData(this, h2);
				
				
				for(final String sampleName:h2.getSampleNamesInOrder())
					{
					if(pedigree==null) continue;
					final Pedigree.Person p= pedigree.getUniqPersonById(sampleName);
					if(p==null)
						{
						sampleNotFoundInVcf.add(sampleName);
						LOG.info("Cannot find "+sampleName+" in "+pedigreeFile);
						}
					else
						{
						final Trio trio = new Trio();
						trio.child_id = header.getSampleNameToOffset().getOrDefault(sampleName,-1);
						if(trio.child_id<0) throw new IllegalStateException(sampleName);
						if(p.hasFather()) {
							final Pedigree.Person parent = p.getFather();
							trio.father_id = header.getSampleNameToOffset().getOrDefault(parent.getId(),-1);
							}
						if(p.hasMother()) {
							final Pedigree.Person parent = p.getMother();
							trio.mother_id = header.getSampleNameToOffset().getOrDefault(parent.getId(),-1);
							}
						if(!(trio.father_id==-1 && trio.mother_id==-1)) {
							trios.add(trio);
							}
						}
					}
			
				LOG.info("trios(s) in pedigree: "+trios.size());
				final ProgressFactory.Watcher<VariantContext> progress = 
						ProgressFactory.newInstance().
						dictionary(header).
						logger(LOG).
						build();
				w.writeHeader(h2);
				while(r.hasNext())
					{
					final VariantContext ctx = progress.apply(r.next());
								
					final Set<String> incompatibilities = new HashSet<String>();
					
					
					for(final Trio trio : trios)
						{
						final Genotype gChild = ctx.getGenotype(trio.child_id);
						if(gChild==null) throw new IllegalStateException();
						final Genotype gFather = trio.father_id<0?null:ctx.getGenotype(trio.father_id);
						final Genotype gMother = trio.mother_id<0?null:ctx.getGenotype(trio.mother_id);
						final DeNovoDetector.DeNovoMutation mut = detector.test(ctx, gFather, gMother, gChild);	
						
						if(mut!=null)
							{
							incompatibilities.add(gChild.getSampleName());
							}
						}	
					
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					if(!incompatibilities.isEmpty()) {
						//set filter for samples that are not a mendelian violation
						if(!StringUtil.isBlank(this.genotypeFilterNameNoIncompat))
							{
							vcb.genotypes( ctx.getGenotypes().stream().map(G->
									incompatibilities.contains(G.getSampleName())?
											G:new GenotypeBuilder(G).
											filters(this.genotypeFilterNameNoIncompat).
									make()
									).collect(Collectors.toList()));
							}
						
						
						++count_incompats;
						// set INFO attribute
						vcb.attribute(attributeName, incompatibilities.toArray());
						
						// set FILTER 
						if(!StringUtil.isBlank(this.filterAnyIncompat))
							{
							vcb.filter(this.filterAnyIncompat);
							}
						else if(!ctx.isFiltered())
							{
							vcb.passFilters();
							}
						}
					else//No denovo
						{
						// dicard variant
						if(this.discard_variants_without_mendelian_incompat) {
							continue;
							}
						// set filters
						if(!StringUtil.isBlank(this.filterNoIncompat))
							{
							vcb.filter(this.filterNoIncompat);
							}
						else if(!ctx.isFiltered())
							{
							vcb.passFilters();
							}
						}
					w.add(vcb.make());
					}
				progress.close();
				
				
				LOG.info("incompatibilitie(s) N="+count_incompats);
				if(! sampleNotFoundInVcf.isEmpty())
					{
					LOG.info("SAMPLE(S) not found: "+String.join(" / ",sampleNotFoundInVcf));
					}
				return 0;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			}
	
	
    public VCFTrios()
    	{
    	}

	@Override
	public int doWork(final List<String> args) {
		if(!StringUtil.isBlank(this.filterAnyIncompat) && !StringUtil.isBlank(this.filterNoIncompat))
			{
			LOG.error("Filters no/any incompatibilities both defined.");
			return -1;
			}
		return doVcfToVcf(args, this.outputFile);
		}
	
	
	public static void main(final String[] args)
		{
		new VCFTrios().instanceMainWithExit(args);
		}

	}
