package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Level;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SortingCollection;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.LineReaderUtil;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.MyPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

public class VCFCompare extends AbstractCommandLineProgram
	{
	private Input inputs[]=new Input[]{null,null};
	
	
	
	private static class Venn
		{
		private class Circle
			{
			double x=0;
			int uniq=0;
			
			double weight() { return this.uniq+comm;}
			
			double radius()
				{
				return Math.sqrt(((comm+this.uniq)/weight())/Math.PI);
				}
			double intersectionX(Circle other)
				{
				return 0;
				}
			}
		
		Circle c1=new Circle();
		Circle c2=new Circle();
		int comm=0;
		
		private double intersectionX()
			{
			return 0;
			}
		
		void make()
			{
			double total=comm+c1.uniq+c2.uniq;
			double maxx=(c1.x+c1.radius()+c2.radius());
			double minx=(c1.x+c1.radius()-c2.radius());
			if(comm==0)
				{
				c2.x=maxx;
				}
			else
				{
				for(;;)
					{
					double xmid=intersectionX();
					}
				}
			}
		}
	
	private static interface VariantFilter
		{
		VariantContext filter(VariantContext ctx,int fileIndex);
		public String getSuffix();
		}
	
	private class Input
		{
		String filename;
		VCFHeader header;
		VCFCodec codec;
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
			v.line=dis.readUTF();
			return v;
			}
		@Override
		public void encode(DataOutputStream dos, LineAndFile v)
				throws IOException
			{
			dos.writeInt(v.fileIdx);
			dos.writeUTF(v.line);
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
			int i=ctx1.getChr().compareTo(ctx2.getChr());
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
				default: switch(super.handleOtherOptions(c, getopt))
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
		List<VariantFilter> variantFilters=new ArrayList<VariantFilter>();
		
		variantFilters.add(new VariantFilter()
			{
			@Override
			public String getSuffix()
				{
				return "";//all
				}
			@Override
			public VariantContext filter(VariantContext ctx,int fileIndex) {
				return ctx;
				}
			});
		
		variantFilters.add(new VariantFilter()
			{
			@Override
			public String getSuffix()
				{
				return ".havingID";
				}
			@Override
			public VariantContext filter(VariantContext ctx,int fileIndex) {
				return ctx==null || !ctx.hasID()?null:ctx;
				}
			});
		
		variantFilters.add(new VariantFilter()
			{
			@Override
			public String getSuffix()
				{
				return ".qualgt30";
				}
			@Override
			public VariantContext filter(VariantContext ctx,int fileIndex) {
				return ctx==null || !ctx.hasLog10PError() || ctx.getPhredScaledQual()<30.0?null:ctx;
				}
			});
		
		variantFilters.add(new VariantFilter()
			{
			@Override
			public String getSuffix()
				{
				return ".snp";
				}
			@Override
			public VariantContext filter(VariantContext ctx,int fileIndex) {
				return ctx==null || !ctx.isSNP() ? null:ctx;
				}
			});
		
		variantFilters.add(new VariantFilter()
			{
			@Override
			public String getSuffix()
				{
				return ".indel";
				}
			@Override
			public VariantContext filter(VariantContext ctx,int fileIndex) {
				return ctx==null || !ctx.isIndel() ? null:ctx;
				}
			});
		
		
		XMLStreamWriter w=null;
		InputStream in=null;
		SortingCollection<LineAndFile> variants=null;
		Map<String,Counter<String>> sample2count=new TreeMap<String,Counter<String>>();
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
				this.inputs[i].codec=new VCFCodec();
				this.inputs[i].filename= args[getopt.getOptInd()+i];
				info("Opening "+this.inputs[i].filename);
				in=IOUtils.openURIForReading(this.inputs[i].filename);
				LineReader lr=LineReaderUtil.fromBufferedStream(in);
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
			
			for(String s:commonSamples) sample2count.put(s, new Counter<String>());
			
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
						
						if(contexes_init[0]==null && contexes_init[1]!=null)
							{
							diff.incr("variation.unique.to.second.file");
							}
						for(VariantFilter varfilter:variantFilters)
							{
							VariantContext contexes[]=new VariantContext[]{
									varfilter.filter(contexes_init[0],0),
									varfilter.filter(contexes_init[1],1)
									};
							
							if(contexes[0]==null && contexes[1]==null)
								{
								//both filterered
								continue;
								}
							else if(contexes[0]!=null && contexes[1]==null)
								{
								diff.incr("variation"+varfilter.getSuffix()+".unique.to.first.file");
								}
							else if(contexes[0]==null && contexes[1]==null)
								{
								throw new IllegalStateException();
								}
							else
								{
								diff.incr("variation"+varfilter.getSuffix()+".in.both.files");
								
								for(String sample:commonSamples)
									{
									Counter<String> counter=sample2count.get(sample);
									Genotype g0=contexes[0].getGenotype(sample);
									Genotype g1=contexes[1].getGenotype(sample);
									if(g0.sameGenotype(g1))
										{
										counter.incr("same.genotype"+varfilter.getSuffix());
										counter.incr("same.genotype"+varfilter.getSuffix()+"."+g0.getType().name());
										}
									else
										{
										counter.incr("diff.genotype"+varfilter.getSuffix());
										counter.incr("diff.genotype"+varfilter.getSuffix()+"."+g0.getType().name()+".to."+g1.getType().name());
										}
									
									}
								}
							}
						row.clear();
						}
					if(rec==null) break;
					}
				row.add(rec);
				}
			iter.close();
			
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
			
			/* diff */
			w.writeStartElement("div");
			w.writeStartElement("table");
			w.writeStartElement("thead");
			
			w.writeStartElement("caption");
			w.writeCharacters("Differences");
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
			
			for(String k:diff.keySet())
				{
				w.writeStartElement("tr");
				
				w.writeStartElement("td");
				w.writeCharacters(k);
				w.writeEndElement();
					
				w.writeStartElement("td");
				w.writeCharacters(String.valueOf(diff.count(k)));
				w.writeEndElement();
				
				w.writeEndElement();//tr
				}
			w.writeEndElement();//table
			
			w.writeEndElement();//table
			w.writeEndElement();//div
			
			for(String sample:commonSamples)
				{
				Counter<String> count=sample2count.get(sample);
				w.writeStartElement("div");
				w.writeStartElement("table");
				w.writeStartElement("thead");
				
				w.writeStartElement("caption");
				w.writeCharacters(sample);
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
				w.writeEndElement();//table
				
				w.writeEndElement();//table
				w.writeEndElement();//div
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
