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
package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlType;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
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

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

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

  * [20180704] changing the arguments that are not really clear.

 
END_DOC

 */
@Program(
		name="vcftrio",
		description="Find mendelian incompatibilitie in a VCF",
		keywords={"vcf","mendelian","pedigree"}
		)
public class VCFTrios
	extends Launcher
	{
	private static final  Logger LOG = Logger.build(VCFTrios.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	@XmlType(name="vcftrios")
	@XmlRootElement(name="vcftrios")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
	implements VariantContextWriterFactory
		{
		private class CtxWriter extends DelegateVariantContextWriter
			{
			private long count_incompats=0L;
			private final Map<String,Pedigree.Person> samplename2person = new HashMap<String,Pedigree.Person>();
			private final Set<String> sampleNotFound = new HashSet<>();
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				}
			 
			@Override
			public void writeHeader(final VCFHeader header) {
				
				final VCFHeader h2=new VCFHeader(header);
				final Set<VCFHeaderLine> meta = new HashSet<>();
				meta.add(new VCFInfoHeaderLine(
						CtxWriterFactory.this.attributeName,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Samples with mendelian incompatibilities. Pedigree was : "+CtxWriterFactory.this.pedigreeFile
						));
				meta.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY, true));
				
				if(!StringUtil.isBlank(CtxWriterFactory.this.filterAnyIncompat)) {
					meta.add(new VCFFilterHeaderLine(CtxWriterFactory.this.filterAnyIncompat,
							"Variant contains at least one mendelian incompatibilities"));
					}
				if(!StringUtil.isBlank(CtxWriterFactory.this.filterNoIncompat)) {
					meta.add(new VCFFilterHeaderLine(CtxWriterFactory.this.filterNoIncompat,
							"Variant does not contain any mendelian incompatibilities"));
					}
				
				meta.stream().forEach(H->h2.addMetaDataLine(H));
				
				
				
				for(final String sampleName:h2.getSampleNamesInOrder())
					{
					Pedigree.Person p=null;
					for(final Pedigree.Family f:CtxWriterFactory.this.pedigree.getFamilies())
						{
						for(final Pedigree.Person child:f.getIndividuals())
							{
							if(child.getId().equals(sampleName))
								{
								if(p!=null)
									{
									throw new IllegalArgumentException(sampleName+" found twice in pedigree !");
									}
								p=child;
								}
							}
						}
					if(p==null)
						{
						this.sampleNotFound.add(sampleName);
						LOG.info("Cannot find "+sampleName+" in "+pedigreeFile);
						}
					else
						{
						this.samplename2person.put(sampleName, p);
						}
					}
			
				LOG.info("person(s) in pedigree: "+samplename2person.size());
				
				super.writeHeader(h2);				
				}
			
			@Override
			public void add(final VariantContext ctx) {
				final Function<Genotype,Genotype> convertGT;
				
				if(nocall_to_homref) {
					final List<Allele> homref_alleles = Arrays.asList(ctx.getReference(),ctx.getReference());
					convertGT = G1 -> G1==null?null:(G1.isNoCall()?new GenotypeBuilder(G1.getSampleName(), homref_alleles).make():G1);
					}
				else
					{
					convertGT = G1 -> G1;
					}
				
				
				final VariantContextBuilder vcb= new VariantContextBuilder(ctx);
				final Map<String,Genotype> sample2genotype = new HashMap<>( 
						ctx.getGenotypes().stream().
						collect(Collectors.toMap(G->G.getSampleName(), G->G)));
				
							
				final Set<String> incompatibilities = new HashSet<String>();
				
				
				for(final Pedigree.Person child:this.samplename2person.values())
					{
					final Genotype gChildOriginal = sample2genotype.get(child.getId());
					final Genotype gChild = convertGT.apply(gChildOriginal);
					
					if(gChild==null)
						{
						if(this.sampleNotFound.add(child.getId()))
							{
							LOG.debug("cannot get genotype for child  "+child.getId());
							}
						continue;
						}
					
					
					if(gChild.getPloidy()!=2)
						{
						LOG.warn(getClass().getSimpleName()+" only handle two alleles child:"+ allelesToString(gChild));
						continue;
						}
					
					Pedigree.Person parent=child.getFather();
					Genotype gFather=convertGT.apply(parent==null?null:sample2genotype.get(parent.getId()));
					if(gFather==null && parent!=null)
						{
						if(this.sampleNotFound.add(parent.getId()))
							{
							LOG.debug("cannot get genotype for father  "+parent.getId());
							}
						}
					if(gFather!=null && gFather.isNoCall()) gFather=null;

					if(gFather!=null && gFather.getPloidy()!=2)
						{
						LOG.warn(getClass().getSimpleName()+" only handle two alleles father: "+ allelesToString(gFather));
						gFather=null;
						}
					parent=child.getMother();
					
					Genotype gMother = convertGT.apply(parent==null?null:sample2genotype.get(parent.getId()));
					
					if(gMother==null && parent!=null)
						{
						if(this.sampleNotFound.add(parent.getId()))
							{
							LOG.debug("cannot get genotype for mother  "+parent.getId());
							}
						}
					
					if(gMother!=null && gMother.isNoCall()) gMother=null;
					if(gMother!=null && gMother.getPloidy()!=2)
						{
						LOG.debug(getClass().getSimpleName()+" only handle two alleles mother:"+ allelesToString(gMother));
						gMother=null;
						}
					
					boolean is_ok=true;
					
					
					if(gFather!=null && gMother!=null)
						{
						// could be a deletion of the child
						if(gChild.isNoCall() && gFather.isCalled() && gMother.isCalled())
							{
							is_ok = false;
							}
						else
							{
							is_ok=trio(gChild,gFather,gMother);
							}
						}
					else if(gFather!=null)
						{
						is_ok=duo(gChild,gFather);
						}
					else if(gMother!=null)
						{
						is_ok=duo(gChild,gMother);
						}
					
					if(!is_ok)
						{
						incompatibilities.add(child.getId());
						}
					
					}
				vcb.genotypes(sample2genotype.values());
				
				
				if(!incompatibilities.isEmpty()) {
					//set filter for samples that are not a mendelian violation
					if(!StringUtil.isBlank(CtxWriterFactory.this.genotypeFilterNameNoIncompat))
						{
						vcb.genotypes( ctx.getGenotypes().stream().map(G->
								incompatibilities.contains(G.getSampleName())?
										G:new GenotypeBuilder(G).
										filters(CtxWriterFactory.this.genotypeFilterNameNoIncompat).
								make()
								).collect(Collectors.toList()));
						}
					
					
					++this.count_incompats;
					// set INFO attribute
					vcb.attribute(attributeName, incompatibilities.toArray());
					
					// set FILTER 
					if(!StringUtil.isBlank(CtxWriterFactory.this.filterAnyIncompat))
						{
						vcb.filter(CtxWriterFactory.this.filterAnyIncompat);
						}
					}
				else
					{
					if(CtxWriterFactory.this.discard_variants_without_mendelian_incompat) return;
					if(!StringUtil.isBlank(CtxWriterFactory.this.filterNoIncompat))
						{
						vcb.filter(CtxWriterFactory.this.filterNoIncompat);
						}
					}
				super.add(vcb.make());				
				}
			
			@Override
			public void close() {
				LOG.info("incompatibilitie(s) N="+this.count_incompats);
				if(!this.sampleNotFound.isEmpty())
					{
					LOG.info("SAMPLE(S) not found: "+String.join(" / ",this.sampleNotFound));
					}
				super.close();
				}
			
			
			private String allelesToString(final Genotype g)
		    	{
		    	if(!g.isCalled()) return g.getSampleName()+" not called";
		    	if(!g.isAvailable()) return g.getSampleName()+" not available";
		    	return g.getSampleName()+":"+g.getAlleles().
		    			stream().map(A->A.getDisplayString()).
		    			collect(Collectors.joining(" "));
		    	}
		    
			
			private boolean trio(
					final Genotype gChild,
					final List<Allele> fathers,
					final List<Allele> mothers
					)
				{		
				for(int f=0;f< fathers.size();++f)
					{
					for(int m=0;m< mothers.size();++m)
						{
						final Genotype gt=GenotypeBuilder.create(
								gChild.getSampleName(),
								Arrays.asList(fathers.get(f),mothers.get(m))
								);
						if(gt.sameGenotype(gChild,true)) return true;
						}
					}
				return false;
				}
			
			private boolean trio(final Genotype child,final Genotype father,final Genotype mother)
				{
				return	trio(
						child,
						father.getAlleles(),
						mother.getAlleles()
						);
				}
			
			private boolean duo(
					final Allele child1,final Allele child2,
					final Allele parent1,final Allele parent2
					)
				{
				return	 child1.equals(parent1) ||
						 child1.equals(parent2) ||
						 child2.equals(parent1) ||
						 child2.equals(parent2)
						 ;
				}
			
			private boolean duo(final Genotype child,final Genotype parent)
				{
				return	duo(
						child.getAllele(0),child.getAllele(1),
						parent.getAllele(0),parent.getAllele(1)
						);
				}
			}
		
		@XmlElement(name="pedigree")
		@Parameter(names={"-p","--ped","--pedigree"},description="Pedigree file. "+Pedigree.OPT_DESCRIPTION,required=true)
		private File pedigreeFile = null;
	
		@Parameter(names={"-fo","--filter-out"},description="FILTER name if there is NO mendelian violation.")
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
			
		@XmlTransient
		private Pedigree pedigree=null;
		
		@Override
		public VariantContextWriter open(final VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}
		@Override
		public int initialize() {
			if(this.pedigreeFile==null)
				{
				LOG.error("Pedigree undefined.");
				return -1;
				}
			if(!StringUtil.isBlank(this.filterAnyIncompat) && !StringUtil.isBlank(this.filterNoIncompat))
				{
				LOG.error("Filters no/any incompatibilities both defined.");
				return -1;
				}
			
			try {
				LOG.info("reading pedigree "+this.pedigreeFile);
				this.pedigree=Pedigree.newParser().parse(this.pedigreeFile);
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			return 0;
			}

	}
	

	
    public VCFTrios()
    	{
    	}
	

	
	
	@Override
	public int doVcfToVcf(final String inputName, VcfIterator r, final VariantContextWriter delegate)
		{
		final VariantContextWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(r.getHeader()).logger(LOG);
		out.writeHeader(r.getHeader());
		while(r.hasNext())
			{
			out.add(progress.watch(r.next()));
			}
		out.close();
		progress.finish();
		return 0;
		}
	
	

	@Override
	public int doWork(final List<String> args) {
		try {
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
	
	
	public static void main(String[] args)
		{
		new VCFTrios().instanceMainWithExit(args);
		}

	}
