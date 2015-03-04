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

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;

public class VcfToRdf extends AbstractCommandLineProgram
	{
	private long id_generator=0L;
	private static final String XSD="http://www.w3.org/2001/XMLSchema#";
	private static final String RDF=com.github.lindenb.jvarkit.util.ns.RDF.NS;
	private static final String DC="http://purl.org/dc/elements/1.1/";
	private static final String NS="http://github.com/lindenb/jvarkit/";
	private Set<String> suffixes=new HashSet<String>();
	private PrintWriter w = null;
	
	
	private VcfToRdf()
		{
		}
	
	@Override
	public String getProgramDescription() {
		return "convert VCF to RDF (N3 notation)";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfToRdf";
		}
	
	private void prefix(String pfx,String uri)
		{
		this.w.print("@prefix ");
		this.w.print(pfx);
		this.w.print(" <");
		this.w.print(uri);
		this.w.println("> .");
		this.suffixes.add(pfx);
		}
	
	private void emitResource(Object r)
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
	
	private void emitObject(Object r)
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
				w.println(". ");
				}
			else
				{
				w.print("; ");
				}
			}
		}
	
	private void writeHeader(VCFHeader header,URI source)
		{
		if(source!=null)
			{
			emit(source,"rdf:type","vcf:File");
			}
		
		SAMSequenceDictionary dict=header.getSequenceDictionary();
		if(dict!=null)
			{
			for(SAMSequenceRecord ssr:dict.getSequences())
				{
				emit(URI.create("urn:chrom/"+ssr.getSequenceName()),
					"rdf:type","vcf:Chromosome",
					"dc:title",ssr.getSequenceName(),
					"vcf:length",ssr.getSequenceLength(),
					"vcf:sequenceIndex",ssr.getSequenceIndex()
					);
				}
			}
		
		for(VCFFilterHeaderLine h:header.getFilterLines())
			{
			emit(URI.create("urn:filter/"+h.getID()),
					"rdf:type","vcf:Filter",
					"dc:title",h.getID(),
					"dc:description",h.getValue()
					);
			}
		
		//Sample
		for(String sample:header.getSampleNamesInOrder())
			{
			emit(URI.create("urn:sample/"+sample),
					"rdf:type","vcf:Sample",
					"dc:title",sample 
					);
			
			}

		}
	
	
	private void scanVCF(File filein) throws IOException
		{
		VcfIterator in=null;
		URI source=null;
		
		try {
			if(filein!=null) source= filein.toURI();
			in=VCFUtils.createVcfIteratorFromFile(filein);
			VCFHeader header = in.getHeader();
			VepPredictionParser vepPredictionParser=new VepPredictionParser(header);
			writeHeader(header,source);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			while(in.hasNext())
				{
				if(this.w.checkError()) break;
				VariantContext ctx = progress.watch(in.next()); 
				
				/* Variant */
				URI variant = URI.create("urn:variant/"+ctx.getChr()+":"+ctx.getStart()+":"+ctx.getReference().getBaseString());
				
				
				
				emit(variant,
					"rdf:type","vcf:Variant",
					"vcf:chrom",URI.create("urn:chrom/"+ctx.getChr()),
					"vcf:position",ctx.getStart(),
					"vcf:ref",ctx.getReference().getBaseString(),
					"vcf:id",(ctx.hasID()?ctx.getID():null),
					"vcf.qual",(ctx.hasLog10PError()?ctx.getPhredScaledQual():null)
					);
				
				for(Allele alt: ctx.getAlternateAlleles())
					{
					emit(variant,"vcf:alt",alt.getBaseString());
					}
				
				for(String f:ctx.getFilters())
					{
					emit(variant,"vcf:filter",URI.create("urn:filter/"+f));
					}
				
				
				for(VepPrediction prediction : vepPredictionParser.getPredictions(ctx))
					{
					List<Object> L=new ArrayList<>();
					L.add("rdf:type");L.add("vep:Prediction");
					L.add("vcf:variant"); L.add(variant);
					L.add("vcf:allele");L.add(prediction.getAllele().getBaseString());
					for(SequenceOntologyTree.Term term:prediction.getSOTerms())
						{
						L.add("vcf:so");
						L.add(URI.create(term.getUri()));
						}
					if(prediction.getEnsemblTranscript()!=null)
						{
						URI transcriptid=URI.create("http://www.ensembl.org/id/"+prediction.getEnsemblTranscript());
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
								URI.create("http://purl.uniprot.org/ensembl/"+prediction.getEnsemblProtein())
								);
							}
						}
					
					
					
					emit(
						URI.create("urn:vep/"+(++id_generator)),
						L.toArray()
						);
					}
				
				for(String sample: ctx.getSampleNames())
					{
					Genotype g = ctx.getGenotype(sample);
					
					
					List<Object> L=new ArrayList<>();
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
						for(Allele a:g.getAlleles())
							{
							L.add("vcf:allele");L.add(a.getBaseString());
							}
						}
					
					emit(
						URI.create("urn:gt/"+ctx.getChr()+":"+ctx.getStart()+":"+ctx.getReference().getBaseString()+":"+sample),
						L.toArray()
						);
					}
				
				
				}
			in.close(); in = null;
			progress.finish();
			} 
		catch (Exception e) {
			
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		try
			{
			this.w= new PrintWriter(System.out);
			prefix("rdf",RDF);
			prefix("dc", DC);
			prefix("vcf", NS);
			prefix("xsd", XSD);
			prefix("vep", "http://www.ensembl.org/info/docs/tools/vep/");
			prefix("uniprot", "http://purl.uniprot.org/core/");//used in RDF dump of uniprot
			
			
			//prefix("so", "<http://www.sequenceontology.org/browser/current_svn/term/");

			
			if(opt.getOptInd()==args.length)
				{
				scanVCF(null);
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i )
					{
					scanVCF(new File(args[i]));
					}
				}
			
			w.flush();
			w.close();
			w=null;
			return 0;
			}
		catch(Exception err)
			{
			error(err);
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
