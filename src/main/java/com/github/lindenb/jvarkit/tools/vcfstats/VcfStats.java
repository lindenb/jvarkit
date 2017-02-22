/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.InputStream;
import java.io.OutputStream;
import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfIteratorImpl;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser.SnpEffPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

public class VcfStats extends AbstractVcfStats
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfStats.class);

	private static final int HISTOGRAM_STEP=5;
	private VepPredictionParser vepPredictionParser=null;
	private SnpEffPredictionParser snpEffPredictionParser=null;
	private final SequenceOntologyTree.Term coding_exon_variant=SequenceOntologyTree.getInstance().getTermByAcn("SO:0001791");
	private boolean DP_info_is_depth=false;
	
	private class Stats
		{
		Counter<String> counter=new Counter<String>();
		Counter<Integer> alternate_alleles=new Counter<Integer>();
		Counter<Integer> depth=new Counter<Integer>();
		Counter<Integer> qual=new Counter<Integer>();
		Counter<SequenceOntologyTree.Term> snpEffSo=new Counter<SequenceOntologyTree.Term>();
		Counter<SequenceOntologyTree.Term> vepSo=new Counter<SequenceOntologyTree.Term>();

		
		private void watch(boolean is_in_coding,String prefix,Allele allele1,Allele allele2)
			{
			for(int i=0;i< 2;++i)
				{
				Character a1= simpleATGC(allele1);
				Character a2= simpleATGC(allele2);
				if(a1==null || a2==null)
					{
					if(allele1.getDisplayString().length()!=allele2.getDisplayString().length())
						{
						counter.incr(prefix+".indels");
						}
					}
				else
					{
					counter.incr(prefix+".substitutions");
					if(isTransition(a1, a2))
						{
						counter.incr(prefix+".transitions");
						}
					else if(isTransversion(a1, a2))
						{
						counter.incr(prefix+".transversions");
						}
					}
				if(!is_in_coding) break;
				prefix=prefix+".coding";
				}
			}
		
		void watch(String sampleName,VariantContext ctx)
			{
			List<Allele> alleles=null;
			boolean is_in_coding=false;
			Set<SequenceOntologyTree.Term> vepSet=new HashSet<SequenceOntologyTree.Term >();
			Set<SequenceOntologyTree.Term> snpEffSet=new HashSet<SequenceOntologyTree.Term >();
			
			
			
			for(SnpEffPrediction eff :snpEffPredictionParser.getPredictions(ctx))
				{
				for(SequenceOntologyTree.Term t:eff.getSOTerms())
					{
					if(t.equals(coding_exon_variant)) is_in_coding=true;
					snpEffSet.add(t);
					}
				}
			for(VepPrediction eff :vepPredictionParser.getPredictions(ctx))
				{
				for(SequenceOntologyTree.Term t:eff.getSOTerms())
					{
					if(t.equals(coding_exon_variant)) is_in_coding=true;
					vepSet.add(t);
					}
				}
			
			for(SequenceOntologyTree.Term t:snpEffSet)
				{
				snpEffSo.incr(t);
				}
			for(SequenceOntologyTree.Term t:vepSet)
				{
				vepSo.incr(t);
				}
			
			if(sampleName!=null)
				{
				Genotype g=ctx.getGenotype(sampleName);
				
				if(!g.isAvailable())
					{
					counter.incr("genotype.unavailable");
					return;
					}
				if(g.isHomRef()) counter.incr("genotype.hom.ref");
				if(g.isHomVar()) counter.incr("genotype.hom.var");
				if(g.isHet()) counter.incr("genotype.het");
				if(g.isHom()) counter.incr("genotype.hom");
				if(g.isFiltered()) counter.incr("genotype.filtered");
				if(g.isMixed()) counter.incr("genotype.is.mixed");
				alleles=g.getAlleles();
				if(!g.isHomRef() && alleles.size()==2)
					{
					watch(is_in_coding,"genotype",alleles.get(0),alleles.get(1));
					}
				
				if(g.hasDP())
					{
					int dp=g.getDP();
					depth.incr(dp/HISTOGRAM_STEP);
					}
				if(g.hasGQ())
					{
					int q=g.getGQ();
					qual.incr(q/HISTOGRAM_STEP);
					}
				}
			else
				{
				
				
				alternate_alleles.incr(ctx.getAlternateAlleles().size());
				if(ctx.getReference()!=null && ctx.getAlternateAlleles().size()==1)
					{
					Allele allele1=ctx.getReference();
					Allele allele2=ctx.getAlternateAllele(0);
					watch(is_in_coding,"variant",allele1,allele2);
					}
				
				
				if(ctx.hasLog10PError())
					{
					int q=(int)ctx.getPhredScaledQual();
					qual.incr(q/HISTOGRAM_STEP);
					}
				if(DP_info_is_depth)
					{
					int dp=ctx.getAttributeAsInt("DP", -1);
					if(dp!=-1)
						{
						depth.incr(dp/HISTOGRAM_STEP);
						}
					}
				
				
				}
			
			
			
			
			counter.incr("num.variants");
			
			if(ctx.isBiallelic())
				{
				counter.incr("bi.allelic");
				}
			else
				{
				counter.incr("not.bi.allelic");
				}
			
			
			if(ctx.hasSymbolicAlleles())
				{
				counter.incr("symbolic.alleles");
				}
			
			if(ctx.hasID())
				{
				counter.incr("variant.having.ID");
				if(ctx.getID().matches("rs[0-9]+"))
					{
					counter.incr("variant.having.rsId");
					}
				}
			
			}
		
		
		void xml(XMLStreamWriter out) throws XMLStreamException
			{
			if(!counter.isEmpty())
				{
				out.writeStartElement("div");
				
				out.writeStartElement("h4");
				out.writeCharacters("Counts");
				out.writeEndElement();

				
				out.writeStartElement("table");
				out.writeAttribute("j:title", "counts");
				
				out.writeStartElement("thead");
				out.writeStartElement("caption");
				out.writeCharacters("Counts");
				out.writeEndElement();
				out.writeStartElement("tr");
				for(String th:new String[]{"Key","Count"})
					{
					out.writeStartElement("th");
					out.writeCharacters(th);
					out.writeEndElement();
					}
				out.writeEndElement();//tr
				out.writeEndElement();//thead

				
				for(String s:counter.keySet())
					{
					out.writeStartElement("tr");
					
					out.writeStartElement("th");
					out.writeCharacters(s);
					out.writeEndElement();
					
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(counter.count(s)));
					out.writeEndElement();
					
					out.writeEndElement();//tr
					}
				out.writeEndElement();//table
				
				out.writeEndElement();//
				}
			if(!alternate_alleles.isEmpty())
				{
				out.writeStartElement("div");
				
				out.writeStartElement("h4");
				out.writeCharacters("Alternate Alleles");
				out.writeEndElement();
				
				
				out.writeStartElement("table");
				out.writeAttribute("j:title", "alt-alleles");
				out.writeStartElement("thead");
				out.writeStartElement("caption");
				out.writeCharacters("Alternate Alleles");
				out.writeEndElement();
				out.writeStartElement("tr");
				for(String th:new String[]{"Num-Alt-Alleles","Count"})
					{
					out.writeStartElement("th");
					out.writeCharacters(th);
					out.writeEndElement();
					}
				out.writeEndElement();//tr
				out.writeEndElement();//thead

				
				out.writeStartElement("tbody");
				for(Integer i:(new TreeSet<Integer>(alternate_alleles.keySet())))
					{
					out.writeStartElement("tr");
					
					out.writeStartElement("th");
					out.writeCharacters(String.valueOf(i));
					out.writeEndElement();

					
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(alternate_alleles.count(i)));
					out.writeEndElement();

					out.writeEndElement();
					}
				out.writeEndElement();//tbody
				out.writeEndElement();//table
				out.writeEndElement();//div
				}
			if(!depth.isEmpty())
				{
				out.writeStartElement("div");
				
				out.writeStartElement("h4");
				out.writeCharacters("Depth");
				out.writeEndElement();

				
				out.writeStartElement("table");
				out.writeAttribute("j:title", "depth");
				out.writeStartElement("thead");
				out.writeStartElement("caption");
				out.writeCharacters("Depth");
				out.writeEndElement();
				out.writeStartElement("tr");
				for(String th:new String[]{"Lower","Upper","Count"})
					{
					out.writeStartElement("th");
					out.writeCharacters(th);
					out.writeEndElement();
					}
				out.writeEndElement();//tr
				out.writeEndElement();//thead

				out.writeStartElement("tbody");
				for(Integer i:(new TreeSet<Integer>(depth.keySet())))
					{
					out.writeStartElement("tr");
					
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(i*HISTOGRAM_STEP));
					out.writeEndElement();

					out.writeStartElement("td");
					out.writeCharacters(String.valueOf((i+1)*HISTOGRAM_STEP));
					out.writeEndElement();

					
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(depth.count(i)));
					out.writeEndElement();

					out.writeEndElement();
					}
				out.writeEndElement();//tbody
				out.writeEndElement();//table
				out.writeEndElement();//div
				}
			else
				{
				out.writeComment("No Depth available");
				}
			
			if(!qual.isEmpty())
				{
				out.writeStartElement("div");
				
				
				
				out.writeStartElement("h4");
				out.writeCharacters("Quality");
				out.writeEndElement();

				
				out.writeStartElement("table");
				out.writeAttribute("j:title", "qualities");
				out.writeStartElement("thead");
				out.writeStartElement("caption");
				out.writeCharacters("Qualities");
				out.writeEndElement();
				out.writeStartElement("tr");
				for(String th:new String[]{"Lower","Upper","Count"})
					{
					out.writeStartElement("th");
					out.writeCharacters(th);
					out.writeEndElement();
					}
				out.writeEndElement();//tr
				out.writeEndElement();//thead

				out.writeStartElement("tbody");
				for(Integer i:(new TreeSet<Integer>(qual.keySet())))
					{
					out.writeStartElement("tr");
					
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(i*HISTOGRAM_STEP));
					out.writeEndElement();

					out.writeStartElement("td");
					out.writeCharacters(String.valueOf((i+1)*HISTOGRAM_STEP));
					out.writeEndElement();

					
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(qual.count(i)));
					out.writeEndElement();

					out.writeEndElement();
					}
				out.writeEndElement();//tbody
				out.writeEndElement();//table
				out.writeEndElement();//div
				}
			else
				{
				out.writeComment("No QUAL available");
				}
			for(int i=0;i< 2;++i)
				{
				String predName=(i==0?"SnpEff":"Vep");
				Counter<SequenceOntologyTree.Term> set=(i==0?this.snpEffSo:this.vepSo);
				if(set.isEmpty())
					{
					out.writeComment("No prediction for "+predName+" available");
					continue;
					}
				
				
				out.writeStartElement("div");
				
				out.writeStartElement("h4");
				out.writeCharacters("Predictions for "+predName);
				out.writeEndElement();

				
				out.writeStartElement("table");
				out.writeAttribute("j:title", predName+":predictions");
				out.writeStartElement("thead");
				out.writeStartElement("caption");
				out.writeCharacters("Predictions for "+predName);
				out.writeEndElement();
				out.writeStartElement("tr");
				for(String th:new String[]{"ACN","label","Count"})
					{
					out.writeStartElement("th");
					out.writeCharacters(th);
					out.writeEndElement();
					}
				out.writeEndElement();//tr
				out.writeEndElement();//thead
				
				out.writeStartElement("tbody");
				for(SequenceOntologyTree.Term t:set.keySet())
					{
					out.writeStartElement("tr");
					
					out.writeStartElement("th");
					out.writeStartElement("a");
					out.writeAttribute("href", "http://www.sequenceontology.org/miso/current_release/term/"+t.getAcn());
					out.writeCharacters(t.getAcn());
					out.writeEndElement();
					out.writeEndElement();

					out.writeStartElement("td");
					out.writeCharacters(t.getLabel());
					out.writeEndElement();

					
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(set.count(t)));
					out.writeEndElement();

					out.writeEndElement();
					}
				out.writeEndElement();//tbody
				out.writeEndElement();//table
				out.writeEndElement();//div
				}
			
			}
		}
	
	// https://en.wikipedia.org/wiki/File:Transitions-transversions-v3.png
	private static boolean isTransversion(Character a1, Character a2)
		{
		if(a1==null || a2==null) return false;
		if(a1=='A' &&  a2=='C') return true;
		if(a1=='C' &&  a2=='A') return true;
		if(a1=='G' &&  a2=='T') return true;
		if(a1=='T' &&  a2=='G') return true;
		return false;
		}

	private static boolean isTransition(Character a1, Character a2)
		{
		if(a1==null || a2==null) return false;
		if(a1=='A' &&  a2=='G') return true;
		if(a1=='G' &&  a2=='A') return true;
		if(a1=='C' &&  a2=='T') return true;
		if(a1=='T' &&  a2=='C') return true;
		return false;
		}

	
	
	private static Character simpleATGC(Allele al)
		{
		if(al==null) return null;
		String s=al.getBaseString().toUpperCase();
		if(s==null || s.equals(".") || s.length()!=1 ) return null;
		switch(s.charAt(0))
			{
			case 'A': case 'T': case 'G': case 'C': return s.charAt(0);
			default: return null;
			}
		}
	
	@Override
	protected Collection<Throwable> call(String filename) throws Exception {
		VcfIterator iter=null;
		XMLStreamWriter xout=null;
		InputStream vcfInputStream=null;
		OutputStream ostream = null;
		try
			{			
			
			if(filename==null)
				{
				filename="stdin";
				vcfInputStream=stdin();
				}
			else
				{
				vcfInputStream=IOUtils.openURIForReading(filename);
				}
			
			LOG.info("Reading from "+filename);
			ostream= super.openFileOrStdoutAsStream();
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			xout=xof.createXMLStreamWriter(ostream,"UTF-8");
			
			
			iter=new VcfIteratorImpl(vcfInputStream);
			final VCFHeader header=iter.getHeader();
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			
			
			VCFInfoHeaderLine vihl=header.getInfoHeaderLine("DP");
			if		(vihl!=null &&
					String.valueOf(vihl.getDescription()).toLowerCase().contains("depth") &&
					vihl.getType()==VCFHeaderLineType.Integer &&
					vihl.getCountType()==VCFHeaderLineCount.INTEGER &&
					vihl.getCount()==1
					)
				{
				this.DP_info_is_depth=true;
				}
			
			this.vepPredictionParser=new VepPredictionParserFactory(header).get();
			this.snpEffPredictionParser=new SnpEffPredictionParserFactory(header).get();
			
			xout.writeStartElement("html");
			xout.writeAttribute("xmlns", "http://www.w3.org/1999/xhtml");
			xout.writeAttribute("xmlns:j", "https://github.com/lindenb/jvarkit/wiki/VcfStats");
			
			xout.writeStartElement("head");
			
			xout.writeEmptyElement("meta");
			xout.writeAttribute("name", "version");
			xout.writeAttribute("content", String.valueOf(getVersion()));
			
			xout.writeEmptyElement("meta");
			xout.writeAttribute("name", "input");
			xout.writeAttribute("content", filename);

			xout.writeEmptyElement("meta");
			xout.writeAttribute("name", "date");
			xout.writeAttribute("content", String.valueOf(new Date()));

			
			xout.writeStartElement("style");
			xout.writeCharacters(
					"table { border-collapse:collapse;width:100%;}\n"+
					"table,th, td { border: 1px solid black; }\n"
					);
			xout.writeEndElement();//style
			
			
			xout.writeEndElement();//head
			
			xout.writeStartElement("body");
			
						
			
			Stats g_stats=new Stats();
			Map<String,Stats> stats_per_samples=new TreeMap<String,Stats>();
			Map<String,Stats> stats_per_chr=new TreeMap<String,Stats>();
			for(String sample:header.getSampleNamesInOrder())
				{
				stats_per_samples.put(sample, new Stats());
				}
			
			g_stats.counter.incr("num.samples",header.getSampleNamesInOrder().size());
			
			g_stats.counter.incr("num.dictionary.chromosomes",
					header.getSequenceDictionary()==null?
					0
					:header.getSequenceDictionary().getSequences().size()
					);
			
			while(iter.hasNext())
				{
				VariantContext ctx=progress.watch(iter.next());
				
				
				
				Stats k_stats=stats_per_chr.get(ctx.getContig());
				if(k_stats==null)
					{
					k_stats=new Stats();
					stats_per_chr.put(ctx.getContig(), k_stats);
					}
				g_stats.watch(null,ctx);
				k_stats.watch(null,ctx);
				for(String sample:header.getSampleNamesInOrder())
					{
					stats_per_samples.get(sample).watch(sample,ctx);
					}
				}
			
			//hyperlink for chromosomes
			xout.writeStartElement("div");
			xout.writeCharacters("Chromosomes : ");
			for(String k:stats_per_chr.keySet())
				{
				xout.writeCharacters(" [");
				xout.writeStartElement("a");
				xout.writeAttribute("href","#chrom_"+k);
				xout.writeCharacters(k);
				xout.writeEndElement();//a
				xout.writeCharacters(" ]");
				}
			xout.writeEndElement();//div
			
			//hyperlink for samples
			xout.writeStartElement("div");
			xout.writeCharacters("Samples : ");
			for(String s:stats_per_samples.keySet())
				{
				xout.writeCharacters(" [");
				xout.writeStartElement("a");
				xout.writeAttribute("href","#sample_"+s);
				xout.writeCharacters(s);
				xout.writeEndElement();//a
				xout.writeCharacters(" ]");
				}
			xout.writeEndElement();//div

			
			
			g_stats.counter.incr("num.seen.chromosomes",stats_per_chr.size());
			
			xout.writeStartElement("div");
			xout.writeAttribute("j:title","General");
			xout.writeStartElement("h2");
			xout.writeCharacters("General");
			xout.writeEndElement();
			g_stats.xml(xout);
			xout.writeEndElement();//div
			
			xout.writeStartElement("div");
			xout.writeAttribute("j:title","Samples");
			xout.writeStartElement("h2");
			xout.writeCharacters("Samples");
			xout.writeEndElement();
			for(String sample:stats_per_samples.keySet())
				{
				xout.writeStartElement("div");
				
				
				
				xout.writeAttribute("j:name",sample);
				xout.writeAttribute("j:sample",sample);
				
				
				xout.writeEmptyElement("a");
				xout.writeAttribute("name","sample_"+sample);

				
				xout.writeStartElement("h3");
				xout.writeCharacters(sample);
				xout.writeEndElement();
				
				stats_per_samples.get(sample).xml(xout);
				xout.writeEndElement();//div
				}
			xout.writeEndElement();//div
			
			
			xout.writeStartElement("div");
			xout.writeAttribute("title","Chromosomes");
			xout.writeStartElement("h2");
			xout.writeCharacters("chromosomes");
			xout.writeEndElement();
			for(String k:stats_per_chr.keySet())
				{
				xout.writeStartElement("div");
				

				
				xout.writeAttribute("j:name",k);
				xout.writeAttribute("j:chromosome",k);
				if(header.getSequenceDictionary()!=null)
					{
					SAMSequenceRecord ssr=header.getSequenceDictionary().getSequence(k);
					xout.writeAttribute("length",String.valueOf(ssr.getSequenceLength()));
					}
				xout.writeEmptyElement("a");
				xout.writeAttribute("name","chrom_"+k);

				xout.writeStartElement("h3");
				xout.writeCharacters(k);
				xout.writeEndElement();				

				stats_per_chr.get(k).xml(xout);
				xout.writeEndElement();//div
				}
			xout.writeEndElement();//div
			
			
			xout.writeEndElement();//body
			xout.writeEndElement();//html
			xout.flush();
			ostream.flush();
			ostream.close();
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(xout);
			CloserUtil.close(ostream);
			CloserUtil.close(iter);
			CloserUtil.close(vcfInputStream);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfStats().instanceMainWithExit(args);
		}
	}
