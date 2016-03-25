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
package com.github.lindenb.jvarkit.tools.vcf2rdf;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;


import java.util.ArrayList;
import java.util.Collection;
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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;

public class VcfToRdf extends AbstractVcfToRdf
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(AbstractVcfToRdf.class);

	private long id_generator=0L;
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
		VcfIterator in=null;
		URI source=null;
		
		try {
			if(filein!=null) source= filein.toURI();
			in=(filein==null?VCFUtils.createVcfIteratorStdin():VCFUtils.createVcfIteratorFromFile(filein));
			final VCFHeader header = in.getHeader();
			final VepPredictionParser vepPredictionParser=new VepPredictionParser(header);
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
				
				for(final Allele alt: ctx.getAlternateAlleles())
					{
					emit(variant,"vcf:alt",alt.getBaseString());
					}
				
				for(final String f:ctx.getFilters())
					{
					emit(variant,"vcf:filter",URI.create("urn:filter/"+f));
					}
				
				
				for(final VepPrediction prediction : vepPredictionParser.getPredictions(ctx))
					{
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
					}
				
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
	public Collection<Throwable> call() throws Exception {
		final List<String > args = super.getInputFiles();
		try
			{
			this.w= super.openFileOrStdoutAsPrintWriter();
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
			return wrapException(err);
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
