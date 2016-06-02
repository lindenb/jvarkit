package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Level;

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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree.Term;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.MyPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

public class VCFCompare extends AbstractCommandLineProgram
	{
	private Input inputs[]=new Input[]{null,null};
	
	private abstract class Venn0
		{
		String title=null;
		Venn0(String title)
			{
			this.title=title;
			}
		VariantContext filter(VariantContext ctx,int file_id)
			{
			return ctx;
			}
		abstract void  visit(VariantContext ctx[]);
		
		
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
		MyPredictionParser myPredictionParser=null;
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
		public int compare(LineAndFile v1, LineAndFile v2)
			{
			VariantContext ctx1=v1.getContext();
			VariantContext ctx2=v2.getContext();
			int i=ctx1.getContig().compareTo(ctx2.getContig());
			if(i!=0) return i;
			i=ctx1.getStart()-ctx2.getStart();
			if(i!=0) return i;
			i=ctx1.getReference().compareTo(ctx2.getReference());
			if(i!=0) return i;
			return 0;
			}
		}
	
	
	
	private VCFCompare()
		{
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/VcfCompare";
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Compares two VCF files.";
		}
	
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -M (int) Max records in RAM. Optional.");
		out.println(" -T (dir) add temporary directory. Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		SortingCollectionFactory<LineAndFile> factory=new SortingCollectionFactory<LineAndFile>();
		GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, super.getGetOptDefault()+"M:T:"))!=-1)
			{
			switch(c)
				{
				case 'M': factory.setMaxRecordsInRAM(Math.max(1,Integer.parseInt(getopt.getOptArg())));break;
				case 'T': super.addTmpDirectory(new File(getopt.getOptArg()));break;
				default: switch(super.handleOtherOptions(c, getopt, args))
					{
					case EXIT_FAILURE: return -1;
					case EXIT_SUCCESS: return 0;
					default:break;
					}
				}
			}
		if(getopt.getOptInd()==args.length)
			{
			System.err.println("VCFs missing.");
			return -1;
			}
		
		if(getopt.getOptInd()+2!=args.length)
			{
			System.err.println("Illegal number or arguments. Expected two VCFs");
			return -1;
			}
		
		
		

		
		
		XMLStreamWriter w=null;
		InputStream in=null;
		SortingCollection<LineAndFile> variants=null;
		try
			{
			LineAndFileComparator varcmp=new LineAndFileComparator();
			
			factory.setComponentType(LineAndFile.class);
			factory.setComparator(varcmp);
			factory.setTmpDirs(this.getTmpDirectories());
			factory.setCodec(new LineAndFileCodec());
			
			variants=factory.make();
			variants.setDestructiveIteration(true);

			
			for(int i=0;i< 2;++i)
				{
				this.inputs[i]= new Input();
				this.inputs[i].codec=VCFUtils.createDefaultVCFCodec();
				this.inputs[i].filename= args[getopt.getOptInd()+i];
				info("Opening "+this.inputs[i].filename);
				in=IOUtils.openURIForReading(this.inputs[i].filename);
				LineReader lr= new SynchronousLineReader(in);
				LineIterator li=new LineIteratorImpl(lr);
				this.inputs[i].header=(VCFHeader)this.inputs[i].codec.readActualHeader(li);
				
				this.inputs[i].vepPredictionParser=new VepPredictionParser(this.inputs[i].header);
				this.inputs[i].snpEffPredictionParser=new SnpEffPredictionParser(this.inputs[i].header);
				this.inputs[i].myPredictionParser=new MyPredictionParser(this.inputs[i].header);
				
				while(li.hasNext())
					{
					LineAndFile laf=new LineAndFile();
					laf.fileIdx=i;
					laf.line=li.next();
					variants.add(laf);
					}
				info("Done Reading "+this.inputs[i].filename);
				CloserUtil.close(li);
				CloserUtil.close(lr);
				CloserUtil.close(in);
				}
			variants.doneAdding();
			info("Done Adding");
			
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
				venn1List.add(new VennPred("mypred",term)
					{	
					@Override
					Set<Term> terms(VariantContext ctx, int file_id)
						{
						Set<Term> tt=new HashSet<SequenceOntologyTree.Term>();
						for(MyPredictionParser.MyPrediction pred:VCFCompare.this.inputs[file_id].myPredictionParser.getPredictions(ctx))
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
								error("Duplicate context in "+inputs[var.fileIdx].filename+" : "+var.line);
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
			
			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
			w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
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
			}
		catch(Exception err)
			{
			getLogger().log(Level.SEVERE, ""+err.getMessage(),err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
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
