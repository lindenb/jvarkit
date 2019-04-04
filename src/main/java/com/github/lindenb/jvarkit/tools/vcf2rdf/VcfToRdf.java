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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcf2rdf;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;


import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;


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


@Program(name="vcf2rdf",description="convert VCF to RDF (N3 notation)",
keywords={"vcf","rdf"})
public class VcfToRdf extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfToRdf.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-a","--alleles"},description="print ALT alleles")
	private boolean printAlleles = false;

	@Parameter(names={"-f","--filters"},description="print FILTERs")
	private boolean printFilters = false;

	@Parameter(names={"-vep","--vep"},description="print VEP informations")
	private boolean printVep = false;

	@Parameter(names={"-g","--genotypes"},description="print Genotypes informations")
	private boolean printGenotypes = false;

	
	
	
	
	
	//private long id_generator=0L;
	private static final String XSD="http://www.w3.org/2001/XMLSchema#";
	private static final String RDF=com.github.lindenb.jvarkit.util.ns.RDF.NS;
	private static final String DC="http://purl.org/dc/elements/1.1/";
	private static final String NS="http://github.com/lindenb/jvarkit/";
	private final Set<String> suffixes=new HashSet<String>();
	private PrintWriter w = null;
	
	public VcfToRdf()
		{
		}
	
	private void prefix(final String pfx,final String uri)
		{
		this.w.print("@prefix ");
		this.w.print(pfx);
		this.w.print(": <");
		this.w.print(uri);
		this.w.println("> .");
		this.suffixes.add(pfx);
		}
	
	private void emitResource(final Object r)
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
	
	private void emitObject(final Object r)
		{
		if(r.getClass()==String.class)
			{
			String s=String.class.cast(r);
			int colon = s.indexOf(':');
			if(colon!=-1 && this.suffixes.contains(s.substring(0, colon)))
				{
				emitResource(URI.create(s));
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
			emitResource(r);
			}
		else
			{
			throw new IllegalStateException(r.getClass().toString());
			}
		}

	
	private void emit(URI r,Object...array)
		{
		if(array.length==0 || array.length%2!=0)
			throw new IllegalStateException(String.valueOf(array.length));
		boolean subject_printed=false;
		
		for(int i=0;i< array.length;i+=2)
			{
			if(array[i+0]==null || array[i+1]==null) continue;
			if(!subject_printed)
				{
				emitResource(r);
				subject_printed=true;
				}
			w.print(" ");
			emitResource(array[i+0]);
			w.print(" ");
			emitObject(array[i+1]);
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
	
	private void writeHeader(final VCFHeader header,final URI source)
		{
		if(source!=null)
			{
			emit(source,"rdf:type","vcf:File");
			}
		
		final SAMSequenceDictionary dict=header.getSequenceDictionary();
		if(dict!=null)
			{
			for(final SAMSequenceRecord ssr:dict.getSequences())
				{
				emit(URI.create("urn:chrom/"+ssr.getSequenceName()),
					"rdf:type","vcf:Chromosome",
					"dc:title",ssr.getSequenceName(),
					"vcf:length",ssr.getSequenceLength(),
					"vcf:sequenceIndex",ssr.getSequenceIndex()
					);
				}
			}
		
		for(final VCFFilterHeaderLine h:header.getFilterLines())
			{
			emit(URI.create("urn:filter/"+h.getID()),
					"rdf:type","vcf:Filter",
					"dc:title",h.getID(),
					"dc:description",h.getValue()
					);
			}
		final Pedigree pedigree = Pedigree.newParser().parse(header);
		for(final Pedigree.Person pe:pedigree.getPersons()) {
			final URI sample = URI.create("urn:sample/"+pe.getId());
			emit(sample,
				"rdf:type","foaf:Person",
				"foaf:name",pe.getId(),
				"foaf:family",pe.getFamily().getId()
				);
			if(pe.isMale())  emit(sample,"idt:gender","male");
			if(pe.isFemale())  emit(sample,"idt:gender","female");
			if(pe.isAffected())  emit(sample,"idt:status","affected");
			if(pe.isUnaffected())  emit(sample,"idt:status","unaffected");
			
		}

		//Sample
		for(final String sample:header.getSampleNamesInOrder())
			{
			emit(URI.create("urn:sample/"+sample),
					"rdf:type","vcf:Sample",
					"dc:title",sample 
					);
			
			}

		}
	
	
	private void scanVCF(final File filein) throws IOException
		{
		VCFIterator in=null;
		URI source=null;
		
		try {
			if(filein!=null) source= filein.toURI();
			in=(filein==null?VCFUtils.createVCFIteratorStdin():VCFUtils.createVCFIteratorFromFile(filein));
			final VCFHeader header = in.getHeader();
			
			
			final VepPredictionParser vepPredictionParser=new VepPredictionParserFactory(header).get();
			writeHeader(header,source);
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			while(in.hasNext())
				{
				if(this.w.checkError())
					{
					LOG.warn("I/O interruption");
					break;
					}
				final VariantContext ctx = progress.watch(in.next()); 
				
				/* Variant */
				final URI variant = URI.create("urn:variant/"+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getBaseString());
				
				
				
				emit(variant,
					"rdf:type","vcf:Variant",
					"vcf:chrom",URI.create("urn:chrom/"+ctx.getContig()),
					"vcf:position",ctx.getStart(),
					"vcf:ref",ctx.getReference().getBaseString(),
					"vcf:id",(ctx.hasID()?ctx.getID():null),
					"vcf:qual",(ctx.hasLog10PError()?ctx.getPhredScaledQual():null)
					);
				
				if(this.printAlleles) {
					for(final Allele alt: ctx.getAlternateAlleles())
						{
						emit(variant,"vcf:alt",alt.getBaseString());
						}
				}
				
				if(this.printFilters) {
				for(final String f:ctx.getFilters())
					{
					emit(variant,"vcf:filter",URI.create("urn:filter/"+f));
					}
				}
				
				if(this.printVep) {
				for(final VepPrediction prediction : vepPredictionParser.getPredictions(ctx))
					{
					/* 
					final List<Object> L=new ArrayList<>();
					L.add("rdf:type");L.add("vep:Prediction");
					L.add("vcf:variant"); L.add(variant);
					L.add("vcf:allele");L.add(prediction.getAllele().getBaseString());
					for(final SequenceOntologyTree.Term term:prediction.getSOTerms())
						{
						L.add("vcf:so");
						L.add(URI.create(term.getUri()));
						}
					if(prediction.getEnsemblTranscript()!=null)
						{
						final  URI transcriptid=URI.create("http://www.ensembl.org/id/"+prediction.getEnsemblTranscript());
						L.add("vep:transcript");
						L.add(transcriptid);

						
						if(prediction.getEnsemblGene()!=null)
							{
							emit(transcriptid,
								"uniprot:transcribedFrom",//used  in uniprot dump
								URI.create("http://www.ensembl.org/id/"+prediction.getEnsemblGene())
								);
							}
						
						if(prediction.getEnsemblProtein()!=null)
							{
							emit(
								transcriptid,
								"uniprot:translatedTo",//used  in uniprot dump
								URI.create("http://www.ensembl.org/id/"+prediction.getEnsemblProtein())
								);
							}
						}
					
					
					
					emit(
						URI.create("urn:vep/"+(++id_generator)),
						L.toArray()
						);
					*/
					}
				}
				
				
			if(this.printGenotypes) {
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
						URI.create("urn:gt/"+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getBaseString()+":"+sample),
						L.toArray()
						);
					}
			}
				
				}
			in.close(); in = null;
			progress.finish();
			} 
		catch (final Exception e) {
			throw new IOException(e);
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try
			{
			this.w= super.openFileOrStdoutAsPrintWriter(this.outputFile);
			prefix("rdf",RDF);
			prefix("dc", DC);
			prefix("vcf", NS);
			prefix("xsd", XSD);
			prefix("vep", "http://www.ensembl.org/info/docs/tools/vep/");
			prefix("uniprot", "http://purl.uniprot.org/core/");//used in RDF dump of uniprot
			
			
			//prefix("so", "<http://www.sequenceontology.org/browser/current_svn/term/");

			
			if(args.isEmpty())
				{
				scanVCF(null);
				}
			else
				{
				for(final String file: IOUtils.unrollFiles(args))
					{
					scanVCF(new File(file));
					}
				}
			
			w.flush();
			w.close();
			w=null;
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			}
		
		}
	public static void main(String[] args)
		{
		new VcfToRdf().instanceMainWithExit(args);
		}
	}
