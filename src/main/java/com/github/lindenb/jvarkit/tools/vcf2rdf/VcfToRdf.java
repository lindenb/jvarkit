/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcf2rdf;

import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.rdf.ns.DC;
import com.github.lindenb.jvarkit.rdf.ns.RDFS;
import com.github.lindenb.jvarkit.util.Pedigree;

import htsjdk.variant.vcf.VCFIterator;

import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;


/**

BEGIN_DOC


### Example


```
$  java -jar dist/vcf2rdf.jar < in.vcf | xmllint --format -
```


```xml
<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:vcf="http://github.com/lindenb/jvarkit/" xmlns:xsd="http://www.w3.org/2001/XMLSchema#">
  <vcf:Chromosome rdf:about="urn:chromosome/1">
    <dc:title>1</dc:title>
    <vcf:length rdf:datatype="xsd:int">249250621</vcf:length>
    <vcf:index rdf:datatype="xsd:int">0</vcf:index>
  </vcf:Chromosome>
  <vcf:Chromosome rdf:about="urn:chromosome/2">
    <dc:title>2</dc:title>
    <vcf:length rdf:datatype="xsd:int">243199373</vcf:length>
(...)
  <vcf:Chromosome rdf:about="urn:chromosome/GL000192.1">
    <dc:title>GL000192.1</dc:title>
    <vcf:length rdf:datatype="xsd:int">547496</vcf:length>
    <vcf:index rdf:datatype="xsd:int">83</vcf:index>
  </vcf:Chromosome>
  <vcf:Filter rdf:about="urn:filter/FILTER">
    <dc:title>FILTER</dc:title>
    <dc:description/>
  </vcf:Filter>
  <vcf:Sample rdf:about="urn:sample/V2528">
    <dc:title>V2528</dc:title>
  </vcf:Sample>
  <vcf:Variant rdf:about="urn:variant/1">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467571</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467571</vcf:end>
    <vcf:ref>C</vcf:ref>
    <vcf:alt>C</vcf:alt>
    <vcf:alt>T</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">1311.77</vcf:qual>
    <vcf:BaseQRankSum>-2.441</vcf:BaseQRankSum>
    <vcf:HaplotypeScore>365.0758</vcf:HaplotypeScore>
    <vcf:QD>5.33</vcf:QD>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>49.60</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AC>1</vcf:AC>
    <vcf:FS>88.037</vcf:FS>
    <vcf:MQRankSum>-10.623</vcf:MQRankSum>
    <vcf:ReadPosRankSum>-5.556</vcf:ReadPosRankSum>
    <vcf:Dels>0.02</vcf:Dels>
    <vcf:DP>246</vcf:DP>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/2">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/1"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>C</vcf:allele>
    <vcf:allele>T</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">229</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
    <vcf:pl>1340,0,5602</vcf:pl>
  </vcf:Genotype>
  <vcf:Variant rdf:about="urn:variant/3">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467576</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467576</vcf:end>
    <vcf:ref>T</vcf:ref>
    <vcf:alt>T</vcf:alt>
    <vcf:alt>C</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">624.77</vcf:qual>
    <vcf:BaseQRankSum>-4.704</vcf:BaseQRankSum>
    <vcf:HaplotypeScore>415.3829</vcf:HaplotypeScore>
    <vcf:QD>2.58</vcf:QD>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>48.90</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AC>1</vcf:AC>
    <vcf:FS>26.182</vcf:FS>
    <vcf:MQRankSum>4.246</vcf:MQRankSum>
    <vcf:ReadPosRankSum>0.781</vcf:ReadPosRankSum>
    <vcf:Dels>0.00</vcf:Dels>
    <vcf:DP>242</vcf:DP>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/4">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/3"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>T</vcf:allele>
    <vcf:allele>C</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">230</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
    <vcf:pl>653,0,6159</vcf:pl>
  </vcf:Genotype>
  <vcf:Variant rdf:about="urn:variant/5">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467578</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467581</vcf:end>
    <vcf:ref>TTTC</vcf:ref>
    <vcf:alt>TTTC</vcf:alt>
    <vcf:alt>T</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">2.14748360973E9</vcf:qual>
    <vcf:FS>7.426</vcf:FS>
    <vcf:AC>1</vcf:AC>
    <vcf:BaseQRankSum>-3.400</vcf:BaseQRankSum>
    <vcf:MQRankSum>-7.221</vcf:MQRankSum>
    <vcf:ReadPosRankSum>-3.139</vcf:ReadPosRankSum>
    <vcf:DP>243</vcf:DP>
    <vcf:QD>29.69</vcf:QD>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>49.83</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/6">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/5"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>TTTC</vcf:allele>
    <vcf:allele>T</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">243</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
   <h:pre><![CDATA[ <vcf:pl>5687,0,7242</vcf:pl>
  </vcf:Genotype>
  <vcf:Variant rdf:about="urn:variant/7">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467588</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467588</vcf:end>
    <vcf:ref>T</vcf:ref>
    <vcf:alt>T</vcf:alt>
    <vcf:alt>C</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">610.77</vcf:qual>
    <vcf:BaseQRankSum>5.163</vcf:BaseQRankSum>
    <vcf:HaplotypeScore>106.1491</vcf:HaplotypeScore>
    <vcf:QD>2.51</vcf:QD>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>54.46</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AC>1</vcf:AC>
    <vcf:FS>18.558</vcf:FS>
    <vcf:MQRankSum>3.886</vcf:MQRankSum>
    <vcf:ReadPosRankSum>0.762</vcf:ReadPosRankSum>
    <vcf:Dels>0.00</vcf:Dels>
    <vcf:DP>243</vcf:DP>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/8">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/7"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>T</vcf:allele>
    <vcf:allele>C</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">231</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
    <vcf:pl>639,0,6257</vcf:pl>
  </vcf:Genotype>
  <vcf:Variant rdf:about="urn:variant/9">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467589</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467589</vcf:end>
    <vcf:ref>T</vcf:ref>
    <vcf:alt>T</vcf:alt>
    <vcf:alt>C</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">3705.77</vcf:qual>
    <vcf:BaseQRankSum>6.528</vcf:BaseQRankSum>
    <vcf:HaplotypeScore>79.8219</vcf:HaplotypeScore>
    <vcf:QD>15.19</vcf:QD>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>55.35</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AC>1</vcf:AC>
    <vcf:FS>19.873</vcf:FS>
    <vcf:MQRankSum>2.222</vcf:MQRankSum>
    <vcf:ReadPosRankSum>4.388</vcf:ReadPosRankSum>
    <vcf:Dels>0.00</vcf:Dels>
    <vcf:DP>244</vcf:DP>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/10">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/9"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>T</vcf:allele>
    <vcf:allele>C</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">232</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
    <vcf:pl>3734,0,3449</vcf:pl>
  </vcf:Genotype>
</rdf:RDF>
```


### Example

load with jena TDB


```

(cd /home/lindenb/src/jvarkit-git/ && make vcf2rdf)
export JENAROOT=/home/lindenb/package/apache-jena-2.11.0
 
rm -rf TMPTDB
for F in *.vcf.gz
do
gunzip -c ${F} |\
java -jar ${HOME}/src/jvarkit-git/dist/vcf2rdf.jar > tmp.rdf
${JENAROOT}/bin/tdbloader2 --loc=TMPTDB tmp.rdf
du -hs TMPTDB
done
rm tmp.rdf
${JENAROOT}/bin/tdbdump --loc=TMPTDB | head -n 100

loaded 78,258 variants / 20 Mbytes in Jena/TDB 11,193,250 triples / 892Mbytes 

```





END_DOC
*/


@Program(
	name="vcf2rdf",
	description="convert VCF to RDF (N3 notation)",
	keywords={"vcf","rdf"},
	creationDate = "20191213",
	modificationDate = "20250901",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfToRdf extends Launcher
	{
	private static final Logger LOG = Logger.of(VcfToRdf.class);


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--hide","--remove"},description="hide those features. Comma separated")
	private  String hideStr="FILTER,ALT,VEP,GT";

	
	
	
	
	
	//private long id_generator=0L;
	private static final String XSD=com.github.lindenb.jvarkit.rdf.ns.XSD.NS; ;
	private static final String RDF=com.github.lindenb.jvarkit.rdf.ns.RDF.NS;
	private static final String NS="http://github.com/lindenb/jvarkit/";
	private final Set<String> suffixes=new HashSet<String>();
	
	public VcfToRdf()
		{
		}
	
	private void prefix(final PrintWriter w, String pfx,final String uri)
		{
		w.print("@prefix ");
		w.print(pfx);
		w.print(": <");
		w.print(uri);
		w.println("> .");
		this.suffixes.add(pfx);
		}
	
	private void emitResource(final PrintWriter w,final Object r)
		{
		if(r.getClass()==String.class)
			{
			w.print(String.class.cast(r));
			}
		else if(r.getClass()==URI.class)
			{
			w.print("<");
			w.print(URI.class.cast(r).toString());
			w.print(">");
			}
		else
			{
			throw new IllegalStateException(r.getClass().toString());
			}
		}
	
	private void emitObject(final PrintWriter w,final Object r)
		{
		if(r.getClass()==String.class)
			{
			String s=String.class.cast(r);
			int colon = s.indexOf(':');
			if(colon!=-1 && this.suffixes.contains(s.substring(0, colon)))
				{
				emitResource(w,URI.create(s));
				}
			else
				{
				w.print("\"");
				w.print(s);
				w.print("\"");
				}
			}
		else if(r.getClass()==Integer.class)
			{
			w.print("\"");
			w.print(Integer.class.cast(r));
			w.print("\"^^xsd:int");
			}
		else if(r.getClass()==Double.class)
			{
			w.print("\"");
			w.print(Double.class.cast(r));
			w.print("\"^^xsd:double");
			}
		else if(r.getClass()==URI.class)
			{
			emitResource(w,r);
			}
		else
			{
			throw new IllegalStateException(r.getClass().toString());
			}
		}

	
	private void emit(final PrintWriter w,URI r,Object...array)
		{
		if(array.length==0 || array.length%2!=0)
			throw new IllegalStateException(String.valueOf(array.length));
		boolean subject_printed=false;
		
		for(int i=0;i< array.length;i+=2)
			{
			if(array[i+0]==null || array[i+1]==null) continue;
			if(!subject_printed)
				{
				emitResource(w,r);
				subject_printed=true;
				}
			w.print(" ");
			emitResource(w,array[i+0]);
			w.print(" ");
			emitObject(w,array[i+1]);
			if( i+2 == array.length)
				{
				w.println(" . ");
				}
			else
				{
				w.print(" ; ");
				}
			}
		}
	
	private void writeHeader(final PrintWriter w,final VCFHeader header,final URI source)
		{
		if(source!=null)
			{
			emit(w,source,"rdf:type","vcf:File");
			}
		
		final SAMSequenceDictionary dict=header.getSequenceDictionary();
		if(dict!=null)
			{
			for(final SAMSequenceRecord ssr:dict.getSequences())
				{
				emit(w,URI.create("urn:chrom/"+ssr.getSequenceName()),
					"rdf:type","vcf:Chromosome",
					"dc:title",ssr.getSequenceName(),
					"vcf:length",ssr.getSequenceLength(),
					"vcf:sequenceIndex",ssr.getSequenceIndex()
					);
				}
			}
		
		for(final VCFFilterHeaderLine h:header.getFilterLines())
			{
			emit(w,URI.create("urn:filter/"+h.getID()),
					"rdf:type","vcf:Filter",
					"dc:title",h.getID(),
					"dc:description",h.getValue()
					);
			}
		final Pedigree pedigree = Pedigree.newParser().parse(header);
		for(final Pedigree.Person pe:pedigree.getPersons()) {
			final URI sample = URI.create("urn:sample/"+pe.getId());
			emit(w,sample,
				"rdf:type","foaf:Person",
				"foaf:name",pe.getId(),
				"foaf:family",pe.getFamily().getId()
				);
			if(pe.isMale())  emit(w,sample,"idt:gender","male");
			if(pe.isFemale())  emit(w,sample,"idt:gender","female");
			if(pe.isAffected())  emit(w,sample,"idt:status","affected");
			if(pe.isUnaffected())  emit(w,sample,"idt:status","unaffected");
			
		}

		//Sample
		for(final String sample:header.getSampleNamesInOrder())
			{
			emit(w,URI.create("urn:sample/"+sample),
					"rdf:type","vcf:Sample",
					"dc:title",sample 
					);
			
			}

		}
	
	
	private void scanVCF(final PrintWriter w,URI source,VCFIterator in) throws IOException
		{
		final Set<String> hides = new HashSet<>(Arrays.asList(CharSplitter.COMMA.split(this.hideStr.toUpperCase())));
		final boolean print_ALT_alleles = !hides.contains("ALT");
		final boolean print_filters = !(hides.contains("FILTER") || hides.contains("FILTERS"));
		final boolean print_vep = !(hides.contains("CSQ") || hides.contains("VEP"));
		final boolean print_gt = !(hides.contains("GT") || hides.contains("GENOTYPES")|| hides.contains("GENOTYPE"));
		final Map<String,URI> vepCol2uri=new HashMap<>();
		
			final VCFHeader header = in.getHeader();
			
			
			final VepPredictionParser vepPredictionParser=new VepPredictionParserFactory(header).get();
			writeHeader(w,header,source);
			
			if(print_vep && vepPredictionParser.isValid()) {
				for(String cat:vepPredictionParser.getCategories()) {
					URI uri=URI.create("urn:vep/"+cat);
					vepCol2uri.put(cat, uri);
					emit(w,uri,"rdfs:label",cat);
				}
			}
			
			while(in.hasNext())
				{
				
				final VariantContext ctx = in.next(); 
				
				/* Variant */
				final URI variant = URI.create("urn:variant/"+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getBaseString());
				
				
				
				emit(w,variant,
					"rdf:type","vcf:Variant",
					"vcf:chrom",URI.create("urn:chrom/"+ctx.getContig()),
					"vcf:position",ctx.getStart(),
					"vcf:ref",ctx.getReference().getBaseString(),
					"vcf:id",(ctx.hasID()?ctx.getID():null),
					"vcf:qual",(ctx.hasLog10PError()?ctx.getPhredScaledQual():null)
					);
				
				if(print_ALT_alleles) {
					for(final Allele alt: ctx.getAlternateAlleles())
						{
						emit(w,variant,"vcf:alt",alt.getBaseString());
						}
					}
				
				if(print_filters) {
				for(final String f:ctx.getFilters())
					{
					emit(w,variant,"vcf:filter",URI.create("urn:filter/"+f));
					}
				}
				
				if(print_vep && vepPredictionParser.isValid()) {
					for(final VepPrediction prediction : vepPredictionParser.getPredictions(ctx))
							{
							URI preduri=URI.create("urn:vep/"+StringUtils.md5(ctx.getContig()+":"+ctx.getStart()+":"+prediction.getOriginalAttributeAsString()));
							
							emit(w,preduri,
									"rdf:type","vep:Prediction",
									"vcf:variant",variant
									);
							for(String cat:vepPredictionParser.getCategories()) {
								String s = prediction.get(cat);
								if(StringUtils.isBlank(s)) continue;
								emit(w,preduri,vepCol2uri.get(cat),s);
	 							}
						}
					}
				
				
			if(print_gt) {
				for(final String sample: ctx.getSampleNames())
					{
					final Genotype g = ctx.getGenotype(sample);
					final List<Object> L=new ArrayList<>();
					L.add("vcf:sample");  L.add(URI.create("urn:sample/"+sample));
					L.add("vcf:variant"); L.add(variant);
					L.add("rdf:type");L.add("vcf:Genotype");
					if(g.hasDP())
						{
						L.add("vcf:dp");
						L.add(g.getDP());
						}
					if(g.hasGQ())
						{
						L.add("vcf:gq");
						L.add(g.getGQ());
						}
				
					if(g.isCalled())
						{
						if(g.isHet())
							{
							if(g.isHetNonRef())
								{
								L.add("rdf:type");L.add("vcf:HetNonRefGenotype");
								}
							else
								{
								L.add("rdf:type");L.add("vcf:HetGenotype");
								}
							}
						else if(g.isHom())
							{
							if(g.isHomRef())
								{
								L.add("rdf:type");L.add("vcf:HomRefGenotype");
								}
							else
								{
								L.add("rdf:type");L.add("vcf:HomVarGenotype");
								}
							}
						for(final Allele a:g.getAlleles())
							{
							L.add("vcf:allele");L.add(a.getBaseString());
							}
						}
					
					emit(
						w,
						URI.create("urn:gt/"+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getBaseString()+":"+sample),
						L.toArray()
						);
					}
			}
				
				}
		
		}
	
	@Override
	public int doWork(List<String> args) {
		try
			{
				try(PrintWriter w= super.openPathOrStdoutAsPrintWriter(this.outputFile) ) {
				prefix(w,"rdf",RDF);
				prefix(w,"rdfs",RDFS.NS);
				prefix(w,"dc", DC.NS);
				prefix(w,"vcf", NS);
				prefix(w,"xsd", XSD);
				prefix(w,"vep", "http://www.ensembl.org/info/docs/tools/vep/");
				prefix(w,"uniprot", "http://purl.uniprot.org/core/");//used in RDF dump of uniprot
				
				
				//prefix("so", "<http://www.sequenceontology.org/browser/current_svn/term/");
	
				
				if(args.isEmpty())
					{
					try(VCFIterator iter= new BcfIteratorBuilder().open(stdin())) {
						scanVCF(w,null,iter);
						}
					
					}
				else
					{
					for(final String file: IOUtils.unrollStrings(args))
						{
						try(VCFIterator iter= new BcfIteratorBuilder().open(file)) {
							scanVCF(w,
									IOUtils.isRemoteURI(file)?
									IOUtils.toURL(file).toURI():
									Paths.get(file).toUri()
									,iter);
							}
						
						}
					}
				
				w.flush();
				}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	public static void main(final String[] args)
		{
		new VcfToRdf().instanceMainWithExit(args);
		}
	}
