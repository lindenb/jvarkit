package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree.Term;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
/**

## Example

```bash
$ java -jar dist/vcfcompare.jar \
    "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" \
    "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf"
```
<div>
<dl>
<dt>Sample(s) for https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf.</dt>
<dd>
<ol/>
</dd>
<dt>Sample(s) for https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf.</dt>
<dd>
<ol/>
</dd>
<dt>Common Sample(s).</dt>
<dd>
<ol>
<li>M10475</li>
<li>M10478</li>
<li>M10500</li>
<li>M128215</li>
</ol>
</dd>
</dl>
</div>
<div>
<table>
<thead>
<caption>Differences</caption>
<tr>
<th>Key</th>
<th>Count</th>
</tr>
</thead>
<tbody>
<tr>
<td>variation.unique.to.second.file</td>
<td>6</td>
</tr>
<tr>
<td>variation.unique.to.first.file</td>
<td>6</td>
</tr>
<tr>
<td>variation.in.both.files</td>
<td>3</td>
</tr>
<tr>
<td>count.variations</td>
<td>15</td>
</tr>
</tbody>
</table>
</div>
<div>
<table>
<thead>
<caption>M10475</caption>
<tr>
<th>Key</th>
<th>Count</th>
</tr>
</thead>
<tbody>
<tr>
<td>same.genotype</td>
<td>3</td>
</tr>
<tr>
<td>same.genotype.HOM_REF</td>
<td>2</td>
</tr>
<tr>
<td>same.genotype.HOM_VAR</td>
<td>1</td>
</tr>
</tbody>
</table>
</div>
<div>
<table>
<thead>
<caption>M10478</caption>
<tr>
<th>Key</th>
<th>Count</th>
</tr>
</thead>
<tbody>
<tr>
<td>same.genotype</td>
<td>3</td>
</tr>
<tr>
<td>same.genotype.HOM_REF</td>
<td>1</td>
</tr>
<tr>
<td>same.genotype.HOM_VAR</td>
<td>2</td>
</tr>
</tbody>
</table>
</div>
<div>
<table>
<thead>
<caption>M10500</caption>
<tr>
<th>Key</th>
<th>Count</th>
</tr>
</thead>
<tbody>
<tr>
<td>same.genotype</td>
<td>3</td>
</tr>
<tr>
<td>same.genotype.HOM_REF</td>
<td>1</td>
</tr>
<tr>
<td>same.genotype.HOM_VAR</td>
<td>2</td>
</tr>
</tbody>
</table>
</div>
<div>
<table>
<thead>
<caption>M128215</caption>
<tr>
<th>Key</th>
<th>Count</th>
</tr>
</thead>
<tbody>
<tr>
<td>same.genotype.HET</td>
<td>1</td>
</tr>
<tr>
<td>same.genotype</td>
<td>3</td>
</tr>
<tr>
<td>same.genotype.HOM_REF</td>
<td>1</td>
</tr>
<tr>
<td>same.genotype.HOM_VAR</td>
<td>1</td>
</tr>
</tbody>
</table>
</div>





*/
@Program(name="vcfcompare",description="Compares two VCF files",keywords={"vcf","compare"})
public class VCFCompare extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFCompare.class).make();


	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;
	@ParametersDelegate
	WritingSortingCollection writingSortingCollection=new WritingSortingCollection();

	private Input inputs[]=new Input[]{null,null};
	
	private abstract class Venn0
		{
		final String title;
		Venn0(final String title)
			{
			this.title=title;
			}
		VariantContext filter(final VariantContext ctx,int file_id)
			{
			return ctx;
			}
		abstract void  visit(final VariantContext ctx[]);
		
		
		void write(XMLStreamWriter w) throws XMLStreamException
			{
			w.writeStartElement("div");
			w.writeStartElement("table");
			w.writeStartElement("thead");
			
			w.writeStartElement("caption");
			w.writeCharacters(String.valueOf(title));
			w.writeEndElement();//caption	
			
			w.writeStartElement("tr");
			for(String k:new String[]{"Key","Count"})
				{
				w.writeStartElement("th");
				w.writeCharacters(k);
				w.writeEndElement();
				}
			w.writeEndElement();//tr	
			w.writeEndElement();//thead
			w.writeStartElement("tbody");
			/*
			for(String k:count.keySet())
				{
				w.writeStartElement("tr");
				
				w.writeStartElement("td");
				w.writeCharacters(k);
				w.writeEndElement();
					
				w.writeStartElement("td");
				w.writeCharacters(String.valueOf(count.count(k)));
				w.writeEndElement();
				
				w.writeEndElement();//tr
				}
			w.writeEndElement();//tbody 
			*/
			
			w.writeEndElement();//table
			w.writeEndElement();//div
			}
		}
	
	private class Venn1 extends Venn0
		{
		int uniq[]=new int[]{0,0};
		int comm=0;
		
		Venn1(String title)
			{
			super(title);
			}
		
		
		void  visit(VariantContext ctx[])
			{
			VariantContext ctx0=filter(ctx[0], 0);
			VariantContext ctx1=filter(ctx[1], 1);
			if(ctx0!=null)
				{
				if(ctx1!=null)
					{
					comm++;
					}
				else
					{
					uniq[0]++;
					}
				}
			else if(ctx1!=null)
				{
				uniq[1]++;
				}
			}
		String getSampleName()
			{
			return "*";
			}
		void write(XMLStreamWriter w) throws XMLStreamException
			{
			w.writeStartElement("tr");
			
			w.writeStartElement("th");
			w.writeCharacters(String.valueOf(title));
			w.writeEndElement();
			
			w.writeStartElement("th");
			w.writeCharacters(getSampleName());
			w.writeEndElement();
			
			
			
			w.writeEndElement();
			}
		}
	
	private  class VennType extends Venn1
		{
		private VariantContext.Type type;
		public VennType(VariantContext.Type type)
			{
			super(type.name());
			this.type=type;
			}
		@Override
		VariantContext filter(VariantContext ctx, int file_id)
			{
			return ctx==null || !this.type.equals(ctx.getType())?null:ctx;
			}
		}
	
	private abstract class VennPred extends Venn1
		{
		private SequenceOntologyTree.Term term;
		public VennPred(
				String prefix,
				SequenceOntologyTree.Term term)
			{
			super(prefix+" "+term.getAcn());
			this.term=term;
			}
		abstract Set<SequenceOntologyTree.Term> terms(VariantContext ctx,int file_id);
		@Override
		VariantContext filter(VariantContext ctx, int file_id)
			{
			if(ctx==null) return null;
			return terms(ctx,file_id).contains(this.term)?ctx:null;
			}
		}
	
	
	private static class Transition<X>
		{
		X from;
		X to;
		Transition(X from,X to)
			{
			this.from=from;
			this.to=to;
			}
		@Override
		public int hashCode()
			{
			final int prime = 31;
			int result = 1;
			result = prime * result + ((from == null) ? 0 : from.hashCode());
			result = prime * result + ((to == null) ? 0 : to.hashCode());
			return result;
			}
		@Override
		public boolean equals(Object obj)
			{
			if (this == obj) return true;
			if (obj == null) return false;
			if (!(obj instanceof Transition)) return false;
			@SuppressWarnings("unchecked")
			Transition<X> other = (Transition<X>) obj;
			return from.equals(other.from) &&
				   to.equals(other.to);
			}
		@Override
		public String toString()
			{
			return from.toString()+" "+to.toString();
			}
		
		}
	
	private abstract class Venn2 extends Venn0
		{
		private String sample;
		Venn2(String title,String sample)
			{
			super(title);
			this.sample=sample;
			}
		
		abstract void visit(Genotype[] g);
		@Override
		void visit(VariantContext[] ctx)
			{
			VariantContext v0=filter(ctx[0],0);
			if(v0==null) return;
			Genotype g0=v0.getGenotype(this.sample);
			if(g0==null) return;
			
			VariantContext v1=filter(ctx[1],1);
			if(v1==null) return;
			Genotype g1=v1.getGenotype(this.sample);
			if(g1==null) return;
			visit(new Genotype[]{g0,g1});
			}
		}
	private class VennGType extends Venn2
		{
		Counter<Transition<GenotypeType>> count=new Counter<Transition<GenotypeType>>(); 
		VennGType(String sample)
			{
			super("genotypes",sample);
			}
		@Override
		void visit(Genotype[] g)
			{
			count.incr(new Transition<GenotypeType>(g[0].getType(), g[1].getType()));
			}
		}
	
	private class Input
		{
		String filename;
		VCFHeader header;
		AbstractVCFCodec codec;
		SnpEffPredictionParser snpEffPredictionParser=null;
		VepPredictionParser vepPredictionParser=null;
		AnnPredictionParser annPredictionParser=null;
		}
	
	private class LineAndFile
		{
		int fileIdx;
		String line;
		
		private VariantContext _ctx=null;
		VariantContext getContext()
			{
			if(this._ctx==null)
				{
				this._ctx=inputs[this.fileIdx].codec.decode(this.line);
				}
			return this._ctx;
			}
		ContigPosRef getContigPosRef()
			{
			return new ContigPosRef(getContext());
			}
		}
	
	private class LineAndFileCodec extends AbstractDataCodec<LineAndFile>
		{
		@Override
		public LineAndFile decode(DataInputStream dis) throws IOException
			{
			LineAndFile v=new LineAndFile();
			try {
				v.fileIdx=dis.readInt();
			} catch (Exception e) {
				return null;
				}
			v.line=readString(dis);
			return v;
			}
		@Override
		public void encode(DataOutputStream dos, LineAndFile v)
				throws IOException
			{
			dos.writeInt(v.fileIdx);
			writeString(dos,v.line);
			}
		@Override
		public AbstractDataCodec<LineAndFile> clone() {
			return new LineAndFileCodec();
			}
		}

	private class LineAndFileComparator implements Comparator<LineAndFile>
		{
		@Override
		public int compare(final LineAndFile v1, final LineAndFile v2)
			{
			final ContigPosRef ctx1=v1.getContigPosRef();
			final ContigPosRef ctx2=v2.getContigPosRef();
			return ctx1.compareTo(ctx2);
			}
		}
	
	
	
	private VCFCompare()
		{
		}
	
	@Override
	public int doWork(final List<String> args) {
	
		if(args.isEmpty())
			{
			LOG.error("VCFs missing.");
			return -1;
			}
		
		if(args.size()!=2)
			{
			System.err.println("Illegal number or arguments. Expected two VCFs");
			return -1;
			}
		
		
		

		
		PrintWriter pw =null;
		XMLStreamWriter w=null;
		InputStream in=null;
		SortingCollection<LineAndFile> variants=null;
		try
			{
			LineAndFileComparator varcmp=new LineAndFileComparator();
			
			
			variants=SortingCollection.newInstance(LineAndFile.class, new LineAndFileCodec(),
					varcmp,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			variants.setDestructiveIteration(true);

			
			for(int i=0;i< 2;++i)
				{
				this.inputs[i]= new Input();
				this.inputs[i].codec=VCFUtils.createDefaultVCFCodec();
				this.inputs[i].filename= args.get(i);
				LOG.info("Opening "+this.inputs[i].filename);
				in=IOUtils.openURIForReading(this.inputs[i].filename);
				final LineReader lr= new SynchronousLineReader(in);
				final LineIterator li=new LineIteratorImpl(lr);
				this.inputs[i].header=(VCFHeader)this.inputs[i].codec.readActualHeader(li);
				
				this.inputs[i].vepPredictionParser=new VepPredictionParserFactory(this.inputs[i].header).get();
				this.inputs[i].snpEffPredictionParser=new SnpEffPredictionParserFactory(this.inputs[i].header).get();
				this.inputs[i].annPredictionParser=new AnnPredictionParserFactory(this.inputs[i].header).get();
				
				while(li.hasNext())
					{
					LineAndFile laf=new LineAndFile();
					laf.fileIdx=i;
					laf.line=li.next();
					variants.add(laf);
					}
				LOG.info("Done Reading "+this.inputs[i].filename);
				CloserUtil.close(li);
				CloserUtil.close(lr);
				CloserUtil.close(in);
				}
			variants.doneAdding();
			LOG.info("Done Adding");
			
			Set<String> commonSamples=new TreeSet<String>(this.inputs[0].header.getSampleNamesInOrder());
			commonSamples.retainAll(this.inputs[1].header.getSampleNamesInOrder());
			
			
			
			
			
			List<Venn0> venn1List=new ArrayList<VCFCompare.Venn0>();
			venn1List.add(new Venn1("ALL"));
			venn1List.add(new Venn1("having ID")
				{
				@Override
				public VariantContext filter(VariantContext ctx,int fileIndex) {
					return ctx==null || !ctx.hasID()?null:ctx;
					}
				});
			
			venn1List.add(new Venn1("QUAL greater 30")
				{
				@Override
				public VariantContext filter(VariantContext ctx,int fileIndex) {
					return ctx==null || !ctx.hasLog10PError() || ctx.getPhredScaledQual()<30.0?null:ctx;
					}
				});
			
			for(VariantContext.Type t: VariantContext.Type.values())
				{
				venn1List.add(new VennType(t));
				}
			
			
			for(SequenceOntologyTree.Term term:SequenceOntologyTree.getInstance().getTerms())
				{
				venn1List.add(new VennPred("vep",term)
					{	
					@Override
					Set<Term> terms(VariantContext ctx, int file_id)
						{
						Set<Term> tt=new HashSet<SequenceOntologyTree.Term>();
						for(VepPredictionParser.VepPrediction pred:VCFCompare.this.inputs[file_id].vepPredictionParser.getPredictions(ctx))
							{
							tt.addAll(pred.getSOTerms());
							}
						return tt;
						}
					});
				venn1List.add(new VennPred("SnpEff",term)
					{	
					@Override
					Set<Term> terms(VariantContext ctx, int file_id)
						{
						Set<Term> tt=new HashSet<SequenceOntologyTree.Term>();
						for(SnpEffPredictionParser.SnpEffPrediction pred:VCFCompare.this.inputs[file_id].snpEffPredictionParser.getPredictions(ctx))
							{
							tt.addAll(pred.getSOTerms());
							}
						return tt;
						}
					});
				venn1List.add(new VennPred("ANN",term)
					{	
					@Override
					Set<Term> terms(VariantContext ctx, int file_id)
						{
						Set<Term> tt=new HashSet<SequenceOntologyTree.Term>();
						for(AnnPredictionParser.AnnPrediction pred:VCFCompare.this.inputs[file_id].annPredictionParser.getPredictions(ctx))
							{
							tt.addAll(pred.getSOTerms());
							}
						return tt;
						}
					});
				}
			for(String s:commonSamples) 
				{
				venn1List.add(new VennGType(s));
				}


			/* START : digest results ====================== */
			Counter<String> diff=new Counter<String>();
			List<LineAndFile> row=new ArrayList<LineAndFile>();
			CloseableIterator<LineAndFile> iter=variants.iterator();
			for(;;)
				{
				LineAndFile rec=null;
				if(iter.hasNext())
					{
					rec=iter.next();
					}
				if(rec==null || (!row.isEmpty() && varcmp.compare(row.get(0),rec)!=0))
					{
					if(!row.isEmpty())
						{
						diff.incr("count.variations");
						VariantContext contexes_init[]=new VariantContext[]{null,null};
						for(LineAndFile var:row)
							{
							if(contexes_init[var.fileIdx]!=null)
								{
								LOG.error("Duplicate context in "+inputs[var.fileIdx].filename+" : "+var.line);
								continue;
								}
							contexes_init[var.fileIdx]=var.getContext();
							}
						for(Venn0 venn: venn1List)
							{
							venn.visit(contexes_init);
							}
						
						row.clear();
						}
					if(rec==null) break;
					}
				row.add(rec);
				}
			iter.close();
			/* END : digest results ====================== */
			pw = super.openFileOrStdoutAsPrintWriter(outputFile);
			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
			w= xmlfactory.createXMLStreamWriter(pw);
			w.writeStartElement("html");
			w.writeStartElement("body");
			
			/* specific samples */
			w.writeStartElement("div");
			w.writeStartElement("dl");
			for(int i=0;i< 3;++i)
				{
				String title;
				Set<String> samples;
				switch(i)
					{
					case 0:
					case 1:
						title="Sample(s) for "+this.inputs[i].filename+".";
						samples=new TreeSet<String>(this.inputs[i].header.getSampleNamesInOrder());
						samples.removeAll(commonSamples);
						break;
					default:
						title="Common Sample(s).";
						samples=new TreeSet<String>(commonSamples);
						break;
					}	
				w.writeStartElement("dt");
				w.writeCharacters(title);
				w.writeEndElement();
				w.writeStartElement("dd");
				w.writeStartElement("ol");
				for(String s:samples)
					{
					w.writeStartElement("li");
					w.writeCharacters(s);
					w.writeEndElement();
					}
				w.writeEndElement();
				w.writeEndElement();
				}
			w.writeEndElement();//dl
			w.writeEndElement();//div
			
			for(Venn0 v: venn1List)
				{
				v.write(w);
				}
			
			
			
			w.writeEndElement();//body
			w.writeEndElement();//html
			w.writeEndDocument();
			w.close();w=null;
			pw.flush();
			pw.close();
			pw=null;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(pw);
			if(variants!=null) variants.cleanup();
			}
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VCFCompare().instanceMainWithExit(args);
		}

}
