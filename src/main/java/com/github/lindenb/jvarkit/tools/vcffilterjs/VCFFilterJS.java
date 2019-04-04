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


History:

* 2015 moving to knime and adding predictions snpeff and VEP
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcffilterjs;


import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import javax.script.Bindings;
import javax.script.CompiledScript;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC

## About the script

The user script is a javascript nashorn script [https://docs.oracle.com/javase/8/docs/technotes/guides/scripting/nashorn/api.html](https://docs.oracle.com/javase/8/docs/technotes/guides/scripting/nashorn/api.html).
The return value should be either:


* a boolean : true accept the variant, false reject the variant
* a (VariantContext)[https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html] to replace the current variant
* a (java.util.List)[https://docs.oracle.com/javase/8/docs/api/java/util/List.html]<(VariantContext)[https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html]> to replace the current variant with a list of variants.


## About Galaxy

At first, this tool is not safe for a public Galaxy server, because the javascript code can access the filesystem.
But you can use the JVM parameter

```
-J-Djava.security.manager
```

to prevent it to access the filesystem. See [http://stackoverflow.com/questions/40177810](http://stackoverflow.com/questions/40177810)


## See also

* VcfFilterJdk a faster version that doesn't use javascript but in-memory java compilation.

## History

 * 2017-05 : removed the variables 'sneff', 'vep','casecontrol'... 

## Variables

the script binds the following variables:


 *   **variant** : the current variation;  a org.broadinstitute.variant.variantcontext.VariantContext ( [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) )
 *   **header** : the VCF header org.broadinstitute.variant.vcf.VCFHeader ( [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html) ).
 *   **tools** : an instance of com.github.lindenb.jvarkit.util.vcf.VcfTools ( [https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/VcfTools.java](https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/VcfTools.java) ).
 *   **pedigree** a Pedigree  [https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/Pedigree.java](https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/Pedigree.java):
 *   **individuals** a List<Pedigree.Person> of the persons affected or non-affected and present in the VCF header.



###  Examples


####  Example 

the file filter.js


```
// prints a VARIATION if two samples at least have a DP<200 
function myfilterFunction()
	{
	var samples=header.genotypeSamples;
	var countOkDp=0;


	for(var i=0; i< samples.size();++i)
		{
		var sampleName=samples.get(i);
		if(! variant.hasGenotype(sampleName)) continue;
		var genotype = variant.genotypes.get(sampleName);
		if( ! genotype.hasDP()) continue;
		var dp= genotype.getDP();
		if(dp < 200 ) countOkDp++;
		}
	return (countOkDp>2)
	}
myfilterFunction();

```

```
 $ curl -sL "https://raw.github.com/jamescasbon/PyVCF/master/vcf/test/gatk.vcf" |\
   java -jar  dist/vcffilterjs.jar  -f filter.js
   
##fileformat=VCFv4.1
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BLANK	NA12878	NA12891	NA12892	NA19238NA19239	NA19240
chr22	42526449	.	T	A	151.47	.	.	GT:AD:DP:GQ:PL	0/1:23,8:31:99:190,0,694	0/0:188,0:190:99:0,478,5376	0/0:187,0:187:99:0,493,5322	0/0:247,0:249:99:0,634,6728	0/0:185,0:185:99:0,487,5515	0/0:202,0:202:99:0,520,5857	0/0:181,1:182:99:0,440,5362
chr22	42526634	.	T	C	32.60	.	.	GT:AD:DP:GQ:PL	0/1:21,4:25:71:71,0,702	0/0:187,2:189:99:0,481,6080	0/0:233,0:233:99:0,667,7351	0/0:230,0:230:99:0,667,7394	0/0:174,1:175:99:0,446,5469	0/0:194,2:196:99:0,498,6239	0/0:174,0:175:99:0,511,5894
chr22	42527793	rs1080989	C	T	3454.66	.	.	GT:AD:DP:GQ:PL	./.	0/1:72,90:162:99:1699,0,1767	0/1:103,96:202:99:1756,0,2532	0/0:188,0:188:99:0,526,5889	0/0:160,0:160:99:0,457,4983	0/0:197,0:198:99:0,544,6100	0/0:156,0:156:99:0,439,5041
```


#### Example 


add a new INFO tag 'MYKEY'

the script:

```
var ArrayList = Java.type("java.util.ArrayList");
var VariantContextBuilder = Java.type("htsjdk.variant.variantcontext.VariantContextBuilder");


function addInfo(v)
	{
	var vcb = new VariantContextBuilder(v);
	var atts = new ArrayList();
	atts.add(v.getType().name()+ (variant.isFiltered()?"_FILTERED":"_UNFILTERED"));
	atts.add(v.getType().name()+ (variant.hasID()?"_ID":"_NOID"));
	vcb.attribute("MYKEY",atts);
	return vcb.make();
	}


addInfo(variant);
```

run the program, but first use awk to insert the new INFO definition for 'MYKEY'

```
cat input.vcf |\
	awk '/^#CHROM/ {printf("##INFO=<ID=MYKEY,Number=.,Type=String,Description=\"My key\">\n");} {print;}' |\
	java -jar dist/vcffilterjs.jar -f script.js 
```




####  Example 


Script used for http://plindenbaum.blogspot.fr/2013/10/yes-choice-of-transcripts-and-software.html


```
 function has_missense(v)
	{
	if(!v.getClass().getName().equals("java.lang.String"))
		{
		var i;
		for(i=0;i< v.size();++i)
			{
			if(has_missense(v.get(i))) return true;
			}
		return false;
		}
	if(v.indexOf("non_coding_exon_variant")!=-1) return 0;
	return v.indexOf("missense")!=-1;
	}

function accept(v)
	{
	if(v.isIndel()) return 0;
	var vep=v.getAttribute("CSQ");
	if(vep==null ) return 0;

	var pred=v.getAttribute("PRED");
	if(pred==null ) return 0;	
	if(!has_missense(vep) && has_missense(pred)) return 1;
	return 0;
	}
	
accept(variant);

```





####  Example

Sample having a unique genotype:


```
 function accept(ctx)
	{
	var x,y,g1,g2,count_same=0;
	var sampleList=header.getSampleNamesInOrder();
	// loop over one sample 
	for(x=0;x < sampleList.size();++x)
		{
		g1=ctx.getGenotype( sampleList.get(x) );
		// ignore non-called
		if(! g1.isCalled() ) continue;
		count_same=0;
		// loop over the other samples 
		for(y=0;y< sampleList.size() && count_same==0 ;++y)
			{
			if(x==y) continue;// same sample ?
			g2=ctx.getGenotype( sampleList.get(y) );
			// ignore non-called 
			if(! g2.isCalled() ) continue;
			// is g1==g2 ? 
			if( g1.sameGenotype( g2 ) )
				{
				count_same++;
				}
			}
		// found no other same genotype
		if(count_same==0) return true;
		}
	return false;
	}
accept(variant);

```




```
  curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/vcffilterjs.jar  -f select.js |\
grep -v "#" | cut -f 1-5 

1	13957	rs201747181	TC	T
1	51914	rs190452223	T	G
1	54753	rs143174675	T	G
1	55313	rs182462964	A	T
1	55326	rs3107975	T	C
1	55330	rs185215913	G	A
1	55388	rs182711216	C	T
1	55416	rs193242050	G	A
1	55427	rs183189405	T	C
1	62156	rs181864839	C	T
1	63276	rs185977555	G	A
1	66457	rs13328655	T	A
1	69534	rs190717287	T	C
1	72148	rs182862337	C	T
1	77470	rs192898053	T	C
1	79137	rs143777184	A	T
1	81949	rs181567186	T	C
1	83088	rs186081601	G	C
1	83977	rs180759811	A	G
1	84346	rs187855973	T	C

```


verify:


```
 curl  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" |\
gunzip -c  |\
grep rs201747181 |\
cut -f 10- |\
tr "   " "\n" |\
cut -d ':' -f 1 |\
sort |\
uniq -c

   1013 0|0
     26 0|1
      7 1|0
      1 1|1

```




####  Example

filter homozygotes for sample NA12878


```
 curl -sL "https://raw.github.com/jamescasbon/PyVCF/master/vcf/test/gatk.vcf" |\
java -jar dist/vcffilterjs.jar -e 'variant.getGenotype("NA12878").isHom()'

```


END_DOC
 */
@Program(
		name="vcffilterjs",
		description="Filtering VCF with javascript expressions",
		keywords={"vcf","filter","javascript","json","nashorn"},
		biostars={88921,233587,104021,213032,215885,243972,111924,196057,142215,229935,
				181358,184966,245802,7403,245181,117974,242281,252580},
		references="\"bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files\" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734)."
		)
public class VCFFilterJS
	extends Launcher
	{	
	private static final Logger LOG = Logger.build(VCFFilterJS.class).make();

	private CompiledScript compiledScript = null;
	
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-F","--filter"},description="If not empty, variants won't be discarded and this name will be used in the FILTER column")
	private String filteredTag = "";
	@Parameter(names={"-ped","--pedigree"},description="inject a Pedigree. If it's not specified , the tool will try to extract a pedigree from the VCF header")
	private File pedigreeFile=null;
	@Parameter(names={"-e","--expression"},description=" (js expression). Optional.")
	private String scriptExpr=null;
	@Parameter(names={"-f","--script"},description=" (js file). Optional.")
	private File scriptFile=null;
	@Parameter(names={"-json","--json"},description="json files. syntax key=path/to/file.json . Inject the json object parsed with google gson into the javascript context as 'key'")
	private List<String> jsonFiles=new  ArrayList<>();
	@Parameter(names={"-casecontrol","--casecontrol"},description="deprecated",hidden=true)
	private boolean deprecated_use_casecontrol =false;
	
	public VCFFilterJS()
		{
		
		}
	
	@Override
	protected int doVcfToVcf(final String inputName,final VCFIterator r,final VariantContextWriter w) {
		try
			{
			
			final VCFHeader header = r.getHeader();
			final VcfTools vcfTools = new VcfTools(header);
			
			
			

			final  VCFHeader h2 = new VCFHeader(header);
			addMetaData(h2);
			
			final Pedigree pedigree;
			if(pedigreeFile!=null) {
				pedigree = Pedigree.newParser().parse(this.pedigreeFile);
				}
			else // try to read from VCF header
				{
				Pedigree p=null;
				try {
					p = Pedigree.newParser().parse(header);
					}
				catch(final Exception err) 
					{
					LOG.warn("cannot decode pedigree from vcf header");
					p = Pedigree.createEmptyPedigree();
					}
				pedigree = p;
				}
			
			
			final VCFFilterHeaderLine filterHeaderLine = (filteredTag.trim().isEmpty()?null:
				new VCFFilterHeaderLine(this.filteredTag.trim(),"Filtered with "+getProgramName())
				);
			
			
			if(filterHeaderLine!=null) h2.addMetaDataLine(filterHeaderLine);
			
			final  SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);

			final  Bindings bindings = this.compiledScript.getEngine()
					.createBindings();
			bindings.put("header", header);
			bindings.put("tools", vcfTools);
			bindings.put("pedigree", pedigree);
			

			bindings.put("individuals", Collections.unmodifiableList( pedigree.getPersons().stream().
				filter(P->(P.isAffected() || P.isUnaffected())).
				filter(P->P.hasUniqId()).
				filter(P->header.getSampleNamesInOrder().contains(P.getId())).
				collect(Collectors.toList())))
				;		
				
			
			for(final String jsonkv :this.jsonFiles)
				{
				int eq=jsonkv.indexOf("=");
				if(eq<=0) throw new JvarkitException.UserError("Bad format for json . expected key=/path/to/file.json but got '"+jsonkv+"'");
				final String key=jsonkv.substring(0,eq);
				final FileReader jsonFile = new FileReader(jsonkv.substring(eq+1));
				JsonParser jsonParser=new JsonParser();
				final JsonElement root=jsonParser.parse(jsonFile);
				jsonFile.close();
				bindings.put(key, root);
				}

			w.writeHeader(h2);
			while (r.hasNext() && !w.checkError())
				{
				final  VariantContext variation = progress.watch(r.next());
				bindings.put("variant", variation);

				final Object result = compiledScript.eval(bindings);
				// result is an array of a collection of variants
				if(result!=null && (result.getClass().isArray() || (result instanceof Collection)))
					{
					final  Collection<?> col;
					if(result.getClass().isArray())
						{
						final Object array[]=(Object[])result;
						col= Arrays.asList(array);
						}
					else
						{
						col =( Collection<?>)result;
						}
					// write all of variants
					for(final Object item:col)
						{
						if(item==null) throw new JvarkitException.UserError("item in array is null");
						if(!(item instanceof VariantContext)) throw new JvarkitException.UserError("item in array is not a VariantContext "+item.getClass());
						w.add(VariantContext.class.cast(item));
						}
					}
				// result is a VariantContext
				else if(result!=null && (result instanceof VariantContext)) {
					w.add(VariantContext.class.cast(result));
					}
				else
					{
					boolean accept=true;
					if(result==null)
						{
						accept=false;
						}
					else if(result instanceof Boolean)
						{
						if(Boolean.FALSE.equals(result)) accept = false;
						}
					else if(result instanceof Number)
						{
						if(((Number)result).intValue()!=1) accept = false;
						}
					else
						{
						LOG.warn("Script returned something that is not a boolean or a number:"+result.getClass());
						accept = false;
						}
					if (!accept)
						{
						if(filterHeaderLine!=null)
							{
							final VariantContextBuilder vcb = new VariantContextBuilder(variation);
							vcb.filter(filterHeaderLine.getID());
							w.add(vcb.make());
							}
						continue;
						}
					
					// set PASS filter if needed
					if(filterHeaderLine!=null && !variation.isFiltered())
						{
						w.add( new VariantContextBuilder(variation).passFilters().make());
						continue;
						}
					
					w.add(variation);
					}
				}
			return RETURN_OK;
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
	
	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			this.compiledScript = super.compileJavascript(this.scriptExpr,this.scriptFile);
			
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			this.compiledScript=null;
			}
		}
	
	
	public static void main(String[] args) throws Exception
		{
		new VCFFilterJS().instanceMainWithExit(args);
		}

	}
