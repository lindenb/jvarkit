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
package com.github.lindenb.jvarkit.tools.vcfcmp;

import htsjdk.samtools.SAMSequenceDictionary;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.Interval;

/**

BEGIN_DOC





### Synopsis




```
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) stdin 
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) file2.vcf(.gz) 

```



both vcf **must** be sorted on CHROM/POS/REF



### Column headers


 *  off_target_only_1 : variant for this sample is found in vcf1 but not in vcf2, but variant is off target
 *  off_target_only_2 : variant for this sample is found in vcf2 but not in vcf1, but variant is off target
 *  off_target_both : variant for this sample is found in both vcf, but it variant off target
 *  unique_to_file_1 :  variant for this sample is found in vcf1 but not in vcf2 (in target)
 *  unique_to_file_1_snp : variant for this sample is found in vcf1 but not in vcf2 (in target) and it is a snp
 *  unique_to_file_1_indel : variant for this sample is found in vcf1 but not in vcf2 (in target) and it is an indel
 *  unique_to_file_2 :  variant for this sample is found in vcf2 but not in vcf1 (in target)
 *  unique_to_file_2_snp: variant for this sample is found in vcf2 but not in vcf1 (in target) and it is a snp
 *  unique_to_file_2_indel: variant for this sample is found in vcf2 but not in vcf1 (in target) and it is an indel
 *  both_missing: there is no genotype for this variant (chrom/pos/ref), while some other samples can have a called genotype.
 *  common_context: a genotype is available for this sample for this context in vcf1 and vcf2.
 *  common_context_snp:  a genotype is available for this sample for this context in vcf1 and vcf2 and it's a snp.
 *  common_context_indel:  a genotype is available for this sample for this context in vcf1 and vcf2 and it's an indel.
 *  common_context_discordant_id:  a genotype is available for this sample for this context in vcf1 and vcf2  but the column ID was not the same.
 *  common_context_discordant_filter:  a genotype is available for this sample for this context in vcf1 and vcf2  but the column FILTER is set in one and not in the other
 *  called_and_same : vcf1 and vcf2 have the same genotype for this sample and variant (chrom/pos/ref).
 *  called_and_same_hom_ref :  vcf1 and vcf2 have the same hom-ref genotype for this sample and variant (chrom/pos/ref).
 *  called_and_same_hom_var :  vcf1 and vcf2 have the same hom-var genotype for this sample and variant (chrom/pos/ref).
 *  called_and_same_het : vcf1 and vcf2 have the same het genotype for this sample and variant (chrom/pos/ref).
 *  called_but_discordant : vcf1 and vcf2 don't have the same genotype for this variant and sample.
 *  called_but_discordant_snp : vcf1 and vcf2 don't have the same genotype for this variant and sample and the variant is a SNP.
 *  called_but_discordant_hom1_het2 :  vcf1 and vcf2 don't have the same genotype for this variant and sample. vcf1 is hom and vcf2 is het
 *  called_but_discordant_het1_hom2 :  vcf1 and vcf2 don't have the same genotype for this variant and sample. vcf1 is het and vcf2 is hom
 *  called_but_discordant_hom1_hom2 :  vcf1 and vcf2 don't have the same genotype for this variant and sample. vcf1 is hom and vcf2 is hom
 *  called_but_discordant_het1_het2 :  vcf1 and vcf2 don't have the same genotype for this variant and sample. vcf1 is het and vcf2 is het
 *  called_but_discordant_others :  vcf1 and vcf2 don't have the same genotype for this variant and sample. other cases.







### Example



```
$ java -jar dist-1.128/vcfcomparecallers.jar  Proj1.samtools.vcf.gz  Proj1.varscan.vcf.gz
#Sample	unique_to_file_1	unique_to_file_1_snp	unique_to_file_1_indel	unique_to_file_2	unique_to_file_2_snp	unique_to_file_2_indel	both_missing	common_context	common_context_snp	common_context_indel	common_context_discordant_id	called_and_same	called_and_same_hom_ref	called_and_same_hom_varcalled_and_same_het	called_but_discordant	called_but_discordant_hom1_het2called_but_discordant_het1_hom2	called_but_discordant_hom1_hom2	called_but_discordant_het1_het2	called_but_discordant_others
B00G5XG	43739	15531	27518	0	10773	11730	2182	558753	535010	22508	55052	1043356	0	26920	41136	3047	698	1993	152	204	0
B00G74M	43629	15445	27503	0	10739	11747	2346	558716	534939	22526	55092	1043355	0	27962	40295	2910	742	1823	164	181	0
B00G5XF	43542	15344	27515	0	10742	11691	2185	559017	535236	22533	55089	1044311	0	26842	40961	2960	809	1821	157	173	0
B00G74L	43705	15461	27543	0	10765	11745	2356	558606	534872	22509	55053	1041955	0	26849	42430	2989	725	1904	175	185	0
B00G5XE	43589	15393	27515	0	10764	11708	2425	558691	534970	22481	55052	1042648	0	27088	41698	2974	746	1906	152	170	0

```





```

$ java -jar dist-1.128/vcfcomparecallers.jar  Proj1.samtools.vcf.gz  Proj1.varscan.vcf.gz | verticalize

>>> 2
$1	#Sample	B00G5XG
$2	unique_to_file_1	43739
$3	unique_to_file_1_snp	15531
$4	unique_to_file_1_indel	27518
$5	unique_to_file_2	0
$6	unique_to_file_2_snp	10773
$7	unique_to_file_2_indel	11730
$8	both_missing	2182
$9	common_context	558753
$10	common_context_snp	535010
$11	common_context_indel	22508
$12	common_context_discordant_id	55052
$13	called_and_same	1043356
$14	called_and_same_hom_ref	0
$15	called_and_same_hom_var	26920
$16	called_and_same_het	41136
$17	called_but_discordant	3047
$18	called_but_discordant_hom1_het2	698
$19	called_but_discordant_het1_hom2	1993
$20	called_but_discordant_hom1_hom2	152
$21	called_but_discordant_het1_het2	204
$22	called_but_discordant_others	0
<<< 2

>>> 3
$1	#Sample	B00G74M
$2	unique_to_file_1	43629
$3	unique_to_file_1_snp	15445
$4	unique_to_file_1_indel	27503
$5	unique_to_file_2	0
$6	unique_to_file_2_snp	10739
$7	unique_to_file_2_indel	11747
$8	both_missing	2346
$9	common_context	558716
$10	common_context_snp	534939
$11	common_context_indel	22526
$12	common_context_discordant_id	55092
$13	called_and_same	1043355
$14	called_and_same_hom_ref	0
$15	called_and_same_hom_var	27962
$16	called_and_same_het	40295
$17	called_but_discordant	2910
$18	called_but_discordant_hom1_het2	742
$19	called_but_discordant_het1_hom2	1823
$20	called_but_discordant_hom1_hom2	164
$21	called_but_discordant_het1_het2	181
$22	called_but_discordant_others	0
<<< 3

>>> 4
$1	#Sample	B00G5XF
$2	unique_to_file_1	43542
$3	unique_to_file_1_snp	15344
$4	unique_to_file_1_indel	27515
$5	unique_to_file_2	0
$6	unique_to_file_2_snp	10742
$7	unique_to_file_2_indel	11691
$8	both_missing	2185
$9	common_context	559017
$10	common_context_snp	535236
$11	common_context_indel	22533
$12	common_context_discordant_id	55089
$13	called_and_same	1044311
$14	called_and_same_hom_ref	0
$15	called_and_same_hom_var	26842
$16	called_and_same_het	40961
$17	called_but_discordant	2960
$18	called_but_discordant_hom1_het2	809
$19	called_but_discordant_het1_hom2	1821
$20	called_but_discordant_hom1_hom2	157
$21	called_but_discordant_het1_het2	173
$22	called_but_discordant_others	0
<<< 4


<<< 4


```





### Example


the following XSLT stylesheet can be used to produce a HTML table for a few differences:



```

<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml" version="1.0">
  <xsl:output method="xml"/>
  <xsl:template match="/">
  	<table>
  	<thead>
  		<tr><th>Category</th><th>Sample</th><th>Variant 1</th><th>Genotype 1</th><th>Variant 2</th><th>Genotype 2</th></tr>
  	</thead>
  	<tbody>
    <xsl:apply-templates select="compare-callers/diff">
    	 <xsl:sort select="@sample" />
    	 <xsl:sort select="@type" /> 
    </xsl:apply-templates>
    </tbody>
    </table>
  </xsl:template>
  
  <xsl:template match="diff">
  	<tr>
  		<td><xsl:value-of select="@type"/></td>
  		<td><xsl:value-of select="@sample"/></td>
  		<td><xsl:apply-templates select="variant[@file='1']"/></td>
  		<td><xsl:apply-templates select="variant[@file='1']/genotype"/></td>
  		<td><xsl:apply-templates select="variant[@file='2']"/></td>
  		<td><xsl:apply-templates select="variant[@file='2']/genotype"/></td>
  	</tr>
  	<xsl:text>
</xsl:text>
  </xsl:template>
 
   <xsl:template match="variant">
    <xsl:value-of select="concat('(',@type,') ',chrom,':',pos,' ',id,' ',ref,'/',alts)"/>
   </xsl:template>
 
  <xsl:template match="genotype">
  	<xsl:value-of select="concat('(',@type,')')"/>
  	<xsl:for-each select="allele">
  	<xsl:if test="position()>1">/</xsl:if>
  	<xsl:value-of select="."/>
  	</xsl:for-each>
  	<xsl:if test="dp"> DP:<xsl:value-of select="dp"/></xsl:if>
   </xsl:template>
 

 
</xsl:stylesheet>

```



This other stylesheet can be used to print the context of the variant for a given BAM (pipe the output to bash )



```

<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml" version="1.0">

  <xsl:output method="text"/>
  <xsl:template match="/">

    <xsl:apply-templates select="compare-callers/diff">
    	 <xsl:sort select="@sample" />
    	 <xsl:sort select="@type" /> 
    </xsl:apply-templates>
  </xsl:template>
  
  <xsl:template match="diff">

  	
  	
  	<xsl:variable name="type" select="@type"/>  	
  	<xsl:variable name="sample" select="@sample"/>  	

  	
  	echo "Sample **<xsl:value-of select="@sample"/>** and category **<xsl:value-of select="@type"/>** at **<xsl:choose>
  		<xsl:when test="variant[@file='1']">
  			<xsl:apply-templates select="variant[@file='1']"/>
  		</xsl:when>
  		<xsl:when test="variant[@file='2']">
  			<xsl:apply-templates select="variant[@file='2']"/>
  		</xsl:when>
  	</xsl:choose>**"
  	
  	echo '```'
  	
  	COLUMNS=30 samtools tview -d T -p "<xsl:choose>
  		<xsl:when test="variant[@file='1']">
  			<xsl:apply-templates select="variant[@file='1']" mode="upstream"/>
  		</xsl:when>
  		<xsl:when test="variant[@file='2']">
  			<xsl:apply-templates select="variant[@file='2']" mode="upstream"/>
  		</xsl:when>
  	</xsl:choose>" path/to/file.bam ref.fasta
  	echo '```'
  	

  </xsl:template>
 
   <xsl:template match="variant">
    <xsl:value-of select="concat(chrom,':',pos)"/>
   </xsl:template>
 
    <xsl:template match="variant" mode="upstream">
    <xsl:value-of select="concat(chrom,':',number(pos)-10)"/>
   </xsl:template>
 
</xsl:stylesheet>

```

END_DOC
*/

@Program(name="vcfcomparecallers",description="Compare two VCFs and print common/exclusive information for each sample/genotype")
public class VcfCompareCallers
	extends Launcher
	{

	private static final Logger LOG = Logger.build(VcfCompareCallers.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-n","--num"},description="number of variants to dump in the example file")
	private int numberOfExampleVariants = 10 ;

	@Parameter(names={"-e","--examplefile"},description="Write a few Variants in this XML file. Optional")
	private File exampleFile = null;

	@Parameter(names={"-B","--bed"},description="Limit to variants in that BED region")
	private File captureFile = null;

	@Parameter(names={"-c","--homref2nocall"},description="Treat HomRef as No Call (created when comparing merged vcf with GATK: there is no homref, everything is nocall)")
	private boolean homRefIsNoCall = false;

	@Parameter(names={"-d","--noindel"},description="Ignore indels (SNP only)")
	private boolean ignoreIndels = false;

	@Parameter(names={"--vertical"},description="Vertical layout, to use with group-by tools")
	private boolean verticalLayout = false;

	
	private enum Category
		{
		off_target_only_1,
		off_target_only_2,
		off_target_both,
		unique_to_file_1,
		unique_to_file_1_snp,
		unique_to_file_1_indel,
		unique_to_file_2,
		unique_to_file_2_snp,
		unique_to_file_2_indel,
		both_missing,
		common_context,
		common_context_snp,
		common_context_indel,
		common_context_discordant_id,
		common_context_discordant_filter,
		called_and_same,
		called_and_same_hom_ref,
		called_and_same_hom_var,
		called_and_same_het,
		called_but_discordant,
		called_but_discordant_snp,
		called_but_discordant_hom1_het2,
		called_but_discordant_het1_hom2,
		called_but_discordant_hom1_hom2,
		called_but_discordant_het1_het2,
		called_but_discordant_others,
		
		}
	
	public VcfCompareCallers()
		{
		}
	
	
	
	private void watch(
			final XMLStreamWriter out,
			final VariantContext ctx0,
			final VariantContext ctx1,
			final Genotype g0,
			final Genotype g1,
			final String sampleName,
			final Counter<Category> count,
			final Category cat
			) throws XMLStreamException
		{
		long n=count.incr(cat);
		if(out==null || n> this.numberOfExampleVariants) return;
		final VariantContext variants[]=new VariantContext[]{ctx0,ctx1};
		final Genotype gts[]=new Genotype[]{g0,g1};
		out.writeStartElement("diff");
		out.writeAttribute("type", cat.name());
		out.writeAttribute("sample",sampleName);
		for(int i=0;i< 2;++i)
			{
			if(variants[i]==null) continue;
			out.writeStartElement("variant");
			out.writeAttribute("file",String.valueOf(i+1));
			out.writeAttribute("type",String.valueOf(variants[i].getType()));
			if(variants[i].isFiltered()) out.writeAttribute("filtered","true");
			
			
			out.writeStartElement("chrom");
			out.writeCharacters(variants[i].getContig());
			out.writeEndElement();

			out.writeStartElement("pos");
			out.writeCharacters(String.valueOf(variants[i].getStart()));
			out.writeEndElement();

			
			
			out.writeStartElement("id");
			out.writeCharacters(variants[i].hasID()?variants[i].getID():".");
			out.writeEndElement();
				
			
			out.writeStartElement("ref");
			out.writeCharacters(String.valueOf(variants[i].getReference()));
			out.writeEndElement();
			
			out.writeStartElement("alts");
			out.writeCharacters(String.valueOf(variants[i].getAlternateAlleles()));
			out.writeEndElement();

			
			
			if(gts[i]!=null)
				{
				out.writeStartElement("genotype");
				out.writeAttribute("type",String.valueOf(gts[i].getType()));
				if(gts[i].isFiltered()) out.writeAttribute("filtered","true");

				for(Allele a:gts[i].getAlleles())
					{
					out.writeStartElement("allele");
					out.writeCharacters(a.toString());
					out.writeEndElement();
					}
				if(gts[i].hasDP())
					{
					out.writeStartElement("dp");
					out.writeCharacters(String.valueOf(gts[i].getDP()));
					out.writeEndElement();
					}
				out.writeEndElement();
				}
			
			out.writeEndElement();
			}
		
	
		
		out.writeEndElement();
		out.writeCharacters("\n");
		}
	
	private void internalTests()
		{
		Allele ref = Allele.create("A",true);
		List<Allele> alleles = new java.util.ArrayList<Allele>();
		alleles.add(ref);
		alleles.add(ref);
		Genotype g=new GenotypeBuilder("sample",alleles).make();
		if(g.isNoCall()) throw new IllegalStateException();
		if(!g.isCalled()) throw new IllegalStateException();
		if(!g.isAvailable()) throw new IllegalStateException();
		if(!g.isHomRef()) throw new IllegalStateException();
		}
	
	@Override
	public int doWork(final List<String> args) {
		final VcfIterator vcfInputs[]=new VcfIterator[]{null,null};
		
		internalTests();
		
		try {
			if(args.size()==1)
				{
				LOG.info("Reading from stdin and "+ args.get(0));
				vcfInputs[0] = VCFUtils.createVcfIteratorStdin();
				vcfInputs[1] = VCFUtils.createVcfIterator( args.get(0));
				}
			else if(args.size()==2)
				{
				LOG.info("Reading from stdin and "+ args.get(0)+" and "+ args.get(1));
				vcfInputs[0] = VCFUtils.createVcfIterator( args.get(0));
				vcfInputs[1] = VCFUtils.createVcfIterator( args.get(1));
				}
			else
				{
				LOG.error("illegal.number.of.arguments");
				return -1;
				}
			return compare(vcfInputs[0],vcfInputs[1]);
			}
		catch (final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfInputs[0]);
			CloserUtil.close(vcfInputs[1]);
			}
		}
		
	public  int compare(VcfIterator vcfIterator1,VcfIterator vcfIterator2)
		{
		PrintWriter exampleWriter=null;
		XMLStreamWriter  exampleOut=null;
		PrintStream pw=null;
		final VCFHeader headers[]=new VCFHeader[]{null,null};
		htsjdk.samtools.util.IntervalTreeMap<Boolean> capture = null;
		final VcfIterator vcfInputs[]=new VcfIterator[]{vcfIterator1,vcfIterator2};
		try {	
			if( this.captureFile !=null )
				{
				LOG.info("Reading "+this.captureFile);
				capture = super.readBedFileAsBooleanIntervalTreeMap(this.captureFile);
				}
			
			for(int i=0;i< vcfInputs.length;++i)
				{
				headers[i] = vcfInputs[i].getHeader();
				}
			/* dicts */
			final SAMSequenceDictionary dict0 = headers[0].getSequenceDictionary();
			final SAMSequenceDictionary dict1 = headers[1].getSequenceDictionary();
			final Comparator<VariantContext> ctxComparator;
			if(dict0==null && dict1==null)
				{
				LOG.info("using createChromPosRefComparator");
				ctxComparator = VCFUtils.createChromPosRefComparator();
				}
			else if(dict0!=null && dict1!=null)
				{
				if( !SequenceUtil.areSequenceDictionariesEqual(dict0, dict1))
					{
					LOG.error("not.the.same.sequence.dictionaries");
					return -1;
					}
				LOG.info("using createTidPosRefComparator");
				ctxComparator = VCFUtils.createTidPosRefComparator(dict0);
				}
			else
				{
				LOG.error("not.the.same.sequence.dictionaries");
				return -1;
				}
			/* samples */
			final Set<String> samples0=new HashSet<>(headers[0].getSampleNamesInOrder());
			final Set<String> samples1=new HashSet<>(headers[1].getSampleNamesInOrder());
			final Set<String> commonSamples= new TreeSet<>(samples0);
			commonSamples.retainAll(samples1);
			
			if(commonSamples.size()!=samples0.size() || commonSamples.size()!=samples1.size())
				{
				LOG.warn("Warning: Not the same samples set. Using intersection of both lists.");
				}
			if(commonSamples.isEmpty())
				{	
				LOG.error("No common samples");
				return -1;
				}
			
			final Map<String, Counter<Category>> sample2info=new HashMap<String, Counter<Category>>(1+commonSamples.size());
			for(final String sampleName:commonSamples)
				{
				sample2info.put(sampleName, new  Counter<Category>());
				}
			
			if(this.exampleFile!=null)
				{
				exampleWriter=new PrintWriter(exampleFile,"UTF-8");
				final XMLOutputFactory xof=XMLOutputFactory.newFactory();
				exampleOut=xof.createXMLStreamWriter(exampleWriter);
				exampleOut.writeStartDocument("UTF-8", "1.0");
				exampleOut.writeStartElement("compare-callers");
				}
			
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(dict0);
			final VariantContext buffer[]=new VariantContext[vcfInputs.length];
			final VariantContext prev[]=new VariantContext[vcfInputs.length];
			for(;;)
				{
				VariantContext smallest=null;
				//refill buffer
				for(int i=0;i< vcfInputs.length;++i)
					{
					if(buffer[i]==null && vcfInputs[i]!=null )
						{
						if(vcfInputs[i].hasNext())
							{
							buffer[i]=vcfInputs[i].peek();
							/* check data are sorted */
							if(prev[i]!=null && ctxComparator.compare(prev[i], buffer[i])>0)
								{
								LOG.error("Input "+(i+1)+"/2 is not sorted"+(
									((i==0 && dict0==null) ||(i==1 && dict1==null))?
									"on chrom/pos/ref":
									"on sequence dictionary"
									)+". got\n"+buffer[i]+"\nafter\n"+prev[i]);
								return -1;
								}
							}
						else
							{
							vcfInputs[i].close();
							vcfInputs[i]=null;
							}							
						}
					
					if(buffer[i]!=null )
						{
						if(smallest==null || ctxComparator.compare(buffer[i],smallest)<0)
							{
							smallest=buffer[i];
							}
						}
					}
				
				if(smallest==null) break;
				
				VariantContext ctx0=null;
				VariantContext ctx1=null;
				Interval interval= null;

				if(buffer[0]!=null && ctxComparator.compare(buffer[0],smallest)==0)
					{
					prev[0] = progress.watch(vcfInputs[0].next());
					ctx0= prev[0];
					buffer[0]=null;
					interval = new Interval(ctx0.getContig(),ctx0.getStart(),ctx0.getEnd());
					}
				if(buffer[1]!=null && ctxComparator.compare(buffer[1],smallest)==0)
					{
					prev[1]= progress.watch(vcfInputs[1].next());
					ctx1= prev[1];
					buffer[1]=null;
					interval = new Interval(ctx1.getContig(),ctx1.getStart(),ctx1.getEnd());
					}
				
				if(this.ignoreIndels)
					{
					if(ctx0!=null && !ctx0.isSNP()) continue;
					if(ctx1!=null && !ctx1.isSNP()) continue;
					}
				
				if( buffer[0]!=null && buffer[1]!=null)
					{
					
					//paranoid check
					if( ctxComparator.compare(buffer[0],buffer[1])!=0) {
						exampleWriter.close();
						throw new IllegalStateException();
						}
					
					}
				
				boolean in_capture=true;
				if(capture!=null && interval!=null)
					{
					in_capture = capture.containsOverlapping(interval);
					}
				
				for(final String sampleName: sample2info.keySet())
					{
					final Counter<Category> sampleInfo=sample2info.get(sampleName);
					
					Genotype g0,g1;
					
					
				
					 g0=(ctx0==null?null:ctx0.getGenotype(sampleName));
					 g1=(ctx1==null?null:ctx1.getGenotype(sampleName));
					
					
					
					
					if(g0!=null && (g0.isNoCall() || !g0.isAvailable())) g0=null;
					if(g1!=null && (g1.isNoCall() || !g1.isAvailable())) g1=null;
					
					/* when merging multiple vcf with GATK. there is no homref
					 * because everything is no call. Problem when comparing 
					 * with multiple called VCF.
					 */
					if(g0!=null && this.homRefIsNoCall && g0.isHomRef()) g0=null;
					if(g1!=null && this.homRefIsNoCall && g1.isHomRef()) g1=null;
					

					if(g0==null && g1==null)
						{
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.both_missing);
						continue;
						}
					else if(g0!=null && g1==null)
						{
						if(!in_capture)  {
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.off_target_only_1);
							continue;
							}
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_1);
						
						if(ctx0.isIndel())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_1_indel);
							}
						else if(ctx0.isSNP())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_1_snp);
							}
						continue;
						}
					else if(g0==null && g1!=null)
						{
						if(!in_capture) {
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.off_target_only_2);
							continue;
							}
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_2);
						if(ctx1.isIndel())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_2_indel);
							}
						else if(ctx1.isSNP())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_2_snp);
							}
						continue;
						}
					else
						{	
						if(!in_capture)
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.off_target_both);
							continue;
							}
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context);
						if(ctx0.isIndel() && ctx1.isIndel())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context_indel);
							}
						else if(ctx0.isSNP() && ctx1.isSNP())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context_snp);
							}
						
						if( (ctx0.hasID() && !ctx1.hasID()) ||
							(!ctx0.hasID() && ctx1.hasID()) ||
							(ctx0.hasID() && ctx1.hasID() && !ctx0.getID().equals(ctx1.getID()))
							)
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context_discordant_id);
							}
						
						if( (ctx0.isFiltered() != ctx1.isFiltered())
							)
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context_discordant_filter);
							}
						
						
						if(g0.sameGenotype(g1))
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same);

							if(g0.isHomRef())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same_hom_ref);
								}
							if(g0.isHomVar())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same_hom_var);
								}
							else if(g0.isHet())
								{	
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same_het);
								}
							}
						else
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant);
							if(ctx0.isSNP() && ctx1.isSNP())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_snp);
								}
							
							if(g0.isHom() && g1.isHet())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_hom1_het2);
								}
							else if(g0.isHet() && g1.isHom())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_het1_hom2);
								}
							else if(g0.isHom() && g1.isHom())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_hom1_hom2);
								}
							else if(g0.isHet() && g1.isHet())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_het1_het2);
								}
							else 
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_others);
								}
							}
						
						}
					}
				}
			progress.finish();
		
			pw  = openFileOrStdoutAsPrintStream(this.outputFile);
			
			if(verticalLayout)
				{
				pw.print("#Sample\tCategory\tCount");
				
				for(final String sample: sample2info.keySet())
					{
					final Counter<Category> count=sample2info.get(sample);
					for(final Category c:Category.values())
						{
						pw.print(sample);
						pw.print('\t');
						pw.print(c.name());
						pw.print('\t');
						pw.println(count.count(c));
						}
					if(pw.checkError()) break;
					}
				}
			else
				{
				pw.print("#Sample");
				for(Category c:Category.values())
					{
					pw.print('\t');
					pw.print(c.name());
					}
				pw.println();
				for(final String sample: sample2info.keySet())
					{
					final Counter<Category> count=sample2info.get(sample);
					pw.print(sample);
					for(Category c:Category.values())
						{
						pw.print('\t');
						pw.print(count.count(c));
						}
					pw.println();
					if(pw.checkError()) break;
					}
				}
			pw.flush();
			
			if(exampleOut!=null)
				{
				exampleOut.writeEndElement();
				exampleOut.writeEndDocument();
				exampleOut.flush();
				exampleOut.close();
				}
			return RETURN_OK;
			} 
		catch (Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(pw);
			CloserUtil.close(exampleWriter);
			}
		}

	
	
	public static void main(String[] args) {
		new VcfCompareCallers().instanceMainWithExit(args);
	}
}
