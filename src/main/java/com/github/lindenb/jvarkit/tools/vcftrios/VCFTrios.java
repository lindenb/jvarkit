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
package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.BufferedReader;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

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

	@Parameter(names={"-p","--ped","--pedigree"},description="Pedigree file. "+Pedigree.OPT_DESCRIPTION,required=true)
	private File pedigreeFile = null;

	@Parameter(names={"-f","--filter"},description="filter name. create a filter in the FILTER column")
	private String filterName = null;

	@Parameter(names={"-if","--inversefilter"},description="inverse FILTER, flag variant having NO mendelian incompat.")
	private boolean inverseFilter = false;

	@Parameter(names={"-gf","--gfilter"},description="genotype filter name. create a filter in the GENOTYPE column")
	private String genotypeFilterName = null;
	
	@Parameter(names={"-A","--attribute"},description="INFO Attribute name")
	private String attributeName = "MENDEL";

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"--discard"},description="Discard variants without mendelian incompatibilities")
	private boolean discard_variants_without_mendelian_incompat=false;
	
	private Pedigree pedigree=null;
	
    public VCFTrios()
    	{
    	}
	
    private static String allelesToString(final Genotype g)
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
	
	
	@Override
	public int doVcfToVcf(final String inputName, VcfIterator r, VariantContextWriter w)
		{
		int count_incompats=0;
		final VCFHeader header=r.getHeader();
		final VCFHeader h2=new VCFHeader(header);
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				this.attributeName,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				"mendelian incompatibilities"
				));
		
		addMetaData(h2);

		if(!StringUtil.isBlank(this.filterName)) {
			h2.addMetaDataLine(new VCFFilterHeaderLine(filterName, "data filtered with VCFTrios"));
		}
		if(!StringUtil.isBlank(this.genotypeFilterName))
			{
			h2.addMetaDataLine(new VCFFormatHeaderLine(
					genotypeFilterName,1,
					VCFHeaderLineType.String,
					"Genotype with mendelian incompatibilities"
					));
			}
		w.writeHeader(h2);
		
		final Map<String,Pedigree.Person> samplename2person=new HashMap<String,Pedigree.Person>(h2.getSampleNamesInOrder().size());
		for(final String sampleName:h2.getSampleNamesInOrder())
			{
			Pedigree.Person p=null;
			for(final Pedigree.Family f:this.pedigree.getFamilies())
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
				LOG.info("Cannot find "+sampleName+" in "+pedigreeFile);
				}
			else
				{
				samplename2person.put(sampleName, p);
				}
			}
		
		LOG.info("persons in pedigree: "+samplename2person.size());
		
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header).logger(LOG);
		while(r.hasNext())
			{
			final VariantContext ctx= progress.watch(r.next());
			
			final VariantContextBuilder vcb= new VariantContextBuilder(ctx);
			final Map<String,Genotype> sample2genotype = new HashMap<>( 
					ctx.getGenotypes().stream().
					collect(Collectors.toMap(G->G.getSampleName(), G->G)));
			
						
			final Set<String> incompatibilities=new HashSet<String>();
			
			
			for(final Pedigree.Person child:samplename2person.values())
				{
				final Genotype gChild=sample2genotype.get(child.getId());
				if(gChild==null)
					{
					LOG.debug("cannot get genotype for child  "+child.getId());
					continue;
					}
				if(gChild.isNoCall())
					{
					continue;
					}
				
				if(gChild.getPloidy()!=2)
					{
					LOG.warn(getClass().getSimpleName()+" only handle two alleles child:"+ allelesToString(gChild));
					continue;
					}
				
				Pedigree.Person parent=child.getFather();
				Genotype gFather=(parent==null?null:sample2genotype.get(parent.getId()));
				if(gFather==null && parent!=null)
					{
					LOG.warn("cannot get genotype for father  "+parent.getId());
					}
				if(gFather!=null && gFather.isNoCall()) gFather=null;

				if(gFather!=null && gFather.getPloidy()!=2)
					{
					LOG.warn(getClass().getSimpleName()+" only handle two alleles father: "+ allelesToString(gFather));
					gFather=null;
					}
				parent=child.getMother();
				
				Genotype gMother=(parent==null?null:sample2genotype.get(parent.getId()));
				
				if(gMother==null && parent!=null)
					{
					LOG.debug("cannot get genotype for mother  "+parent.getId());
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
					is_ok=trio(gChild,gFather,gMother);
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
					if(genotypeFilterName!=null)
						{
						sample2genotype.put(child.getId(),
							new GenotypeBuilder(gChild).filters(genotypeFilterName).make()
							);
						}
					}
				}
			vcb.genotypes(sample2genotype.values());
			
		
			++count_incompats;
			
			if(!incompatibilities.isEmpty()) {
				vcb.attribute(attributeName, incompatibilities.toArray());
				if( filterName!=null && !this.inverseFilter) vcb.filter(filterName);
				}
			else
				{
				if(this.discard_variants_without_mendelian_incompat) continue;
				if( filterName!=null && this.inverseFilter) vcb.filter(filterName);
				}
			w.add(vcb.make());
			
			if(w.checkError()) break;
			}	
		progress.finish();
		LOG.info("incompatibilities N="+count_incompats);
		return 0;
		}
	
	

	@Override
	public int doWork(final List<String> args) {
		if(this.pedigreeFile==null)
			{
			LOG.error("Pedigree undefined.");
			return -1;
			}
		if(this.discard_variants_without_mendelian_incompat && this.inverseFilter)
			{
			LOG.error("Cannot inverse filter and discard variants without problem at the same time");
			return -1;
			}
		
		BufferedReader in = null;
		try {
			LOG.info("reading pedigree "+this.pedigreeFile);
			 in=IOUtils.openFileForBufferedReading(this.pedigreeFile);
			this.pedigree=Pedigree.newParser().parse(in);
			in.close();
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		return doVcfToVcf(args, this.outputFile);
		}
	
	
	public static void main(String[] args)
		{
		new VCFTrios().instanceMainWithExit(args);
		}

	}
