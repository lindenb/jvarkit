package com.github.lindenb.jvarkit.tools.genscan;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReaderUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.genscan.AbstractGeneScan.ChromInfo;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Log;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SortingCollection;

/**
 * GenScan
 *
 */
public class GenScan extends AbstractGeneScan
	{

   // @Option(shortName="IIS",doc="Input is a BAM/SAM file",optional=false)	
 
    
	private List<ChromInfo> chromInfos=new ArrayList<ChromInfo>();
	private Map<String,ChromInfo> chrom2chromInfo=new HashMap<String,ChromInfo>();
	private SortingCollection<DataPoint> dataPoints=null;
	private File faidxFile=null;
	

	
	private class DataPoint implements Comparable<DataPoint>
		{
		int tid;
		int pos;
		double values[]=new double[samples.size()];
		@Override
		public int compareTo(DataPoint other)
			{
			int i=this.tid-other.tid;
			if(i!=0) return i;
			i=this.pos-other.pos;
			if(i!=0) return i;
			return 0;
			}
		}
	
	private class DataPointCodec 
		extends AbstractDataCodec<DataPoint>
		{
		@Override
		public DataPoint decode(DataInputStream dis) throws IOException
			{
			int tid;
			try {  tid=dis.readInt();}
			catch(IOException err) { return null;}
			DataPoint dpt=new DataPoint();
			dpt.tid=tid;
			dpt.pos=dis.readInt();
			for(int i=0;i< dpt.values.length;++i)
				{
				 dpt.values[i]=dis.readDouble();
				}
			return dpt;
			}
		
		@Override
		public void encode(DataOutputStream out, DataPoint o)
				throws IOException
			{
			out.writeInt(o.tid);
			out.writeInt(o.pos);
			for(int i=0;i< o.values.length;++i)
				{
				out.writeDouble(o.values[i]);
				}
			}
		
		@Override
		public DataPointCodec clone()
			{
			return new DataPointCodec();
			}
		}
	
	
	
	/*
	protected int doWork()
		{
		CloseableIterator<DataPoint> iter=null;
		IndexedFastaSequenceFile faidx=null;
		try {
			this.chromInfos=new ArrayList<GenScan.ChromInfo>(dict.size());
			
			
			this.dataPoints= SortingCollection.newInstance(DataPoint.class,
					new DataPointCodec() , 
					new Comparator<DataPoint>()
						{
						@Override
						public int compare(DataPoint p0, DataPoint p1)
							{
							return p0.compareTo(p1);
							}
						},
					super.MAX_RECORDS_IN_RAM); 
			
			if(IN.isEmpty())
				{
				LOG.info("reading from stdin");
				read(System.in);
				}
			else
				{
				for(File f: IN)
					{
					LOG.info("reading "+f);
					InputStream in=IoUtil.openFileForReading(f);
					read(in);
					in.close();
					}
				}
			
			dataPoints.setDestructiveIteration(true);
			dataPoints.doneAdding();
			
			if(DICARD_UNSEEN_CHROM)
				{
				int n=0;
				while(n<this.chromInfos.size())
					{
					if(this.chromInfos.get(n).max_pos < this.chromInfos.get(n).min_pos)
						{
						this.chromInfos.remove(n);
						}
					else
						{
						n++;
						}		
					}
				}
			long genomeSize=0L;
			for(ChromInfo ci:this.chromInfos)
				{
				genomeSize+=ci.getSequenceLength();
				}
			
			double x=this.insets.left;
			for(ChromInfo ci:this.chromInfos)
				{
				ci.x=x;
				ci.width=(ci.getSequenceLength()/(double)genomeSize)*(WIDTH-(insets.left+insets.right));
				x=ci.x+ci.width;
				}
			double log10=Math.pow(10,Math.floor(Math.log10(this.max_value)));
			//max_value=(this.max_value/log10);
			this.max_value=Math.max(this.max_value,Math.ceil(this.max_value/log10)*log10);
			
			log10=Math.pow(10,Math.floor(Math.log10(this.min_value)));
			this.min_value=Math.min(this.min_value,Math.floor(this.min_value/log10)*log10);
			
			
		} catch (Exception e)
			{
			e.printStackTrace();
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(faidx);
			if(dataPoints!=null) dataPoints.cleanup();
			}
		return 0;
		}*/
	
	private Point2D.Double convert(ChromInfo ci,DataPoint dpt,int sample_ids)
		{
		return new Point2D.Double(
			ci.converPosToPixel(dpt.pos),
			(screenSize.getHeight()-insets.bottom)-((dpt.values[sample_ids]-min_value)/(max_value-min_value))*(screenSize.getHeight()-(insets.bottom+insets.top))
			);
		}
	protected void drawPoints(Graphics2D g)
		{
		CloseableIterator<DataPoint> iter=dataPoints.iterator();
		while(iter.hasNext())
			{
			DataPoint dpt=iter.next();
			ChromInfo ci=this.chromInfos.get(dpt.tid);
		
			for(int i=0;i< getSamples().size();++i)
				{
				Point2D.Double  xy=convert(ci,dpt,i);
				
				g.setColor(Color.RED);
				g.draw(new Line2D.Double(xy.x-1, xy.y,xy.x+1, xy.y));
				g.draw(new Line2D.Double(xy.x, xy.y-1,xy.x, xy.y+1));
				}
			}
		iter.close();
		iter=null;
		}
	
	
	boolean firstLineIsHeader=true;
	@SuppressWarnings("resource")
	private void readTsv(InputStream in) throws IOException
		{
		info("reading tsv");
		Pattern tab=Pattern.compile("[\t]");
		
		LineIterator r=new LineIteratorImpl(LineReaderUtil.fromBufferedStream(in));
		boolean foundHeader=false;
		if(firstLineIsHeader)
			{
			if(!r.hasNext()) throw new IOException("First line/header is Missing");
			String line=r.next();
			String tokens[]=tab.split(r.next());
			if(tokens.length<4) throw new IOException("Expected at least 4 columns in "+line);
			for(int c=3;c<tokens.length;++c)
				{
				Sample sample=new Sample();
				sample.sample_id=this.samples.size();
				sample.name=tokens[c];
				this.samples.add(sample);
				}
			foundHeader=true;
			}
		
		
		while(r.hasNext())
			{
			String line=r.next();
			if(line.isEmpty() ) continue;
			String tokens[]=tab.split(line, 4);
			if(tokens.length<4) throw new IOException("Expected at least 4 columns in "+line);
			if(!foundHeader && !firstLineIsHeader)
				{
				for(int c=3;c<tokens.length;++c)
					{
					Sample sample=new Sample();
					sample.sample_id=this.samples.size();
					sample.name="$"+(c-2);
					this.samples.add(sample);
					}
				foundHeader=true;
				continue;
				}
			else
				{
				if(tokens.length!=(3+samples.size()))
					{
					throw new IOException("Expected at least "+(3+samples.size())+" columns in "+line);
					}
				}
			
			DataPoint pt=new DataPoint();
			ChromInfo ci=this.chrom2chromInfo.get(tokens[0]);
			if(ci==null)
				{
				if(this.faidxFile!=null) //dict provided by user
					{
					warning("chromosome "+tokens[0]+" was not defined in dictionary.");
					continue;
					}
				ci=new ChromInfo();
				ci.sequenceName=tokens[0];
				ci.tid=this.chromInfos.size();
				this.chrom2chromInfo.put(tokens[0], ci);
				}
			
			pt.tid=ci.tid;
			
			try {
				pt.pos=Integer.parseInt(tokens[1]);
				if(pt.pos< 0 || pt.pos > ci.getSequenceLength())//TODO only if defined dict
					{
					warning("Position out of range in "+line);
					continue;
					}
				}
			catch (Exception e) {
				warning("bad pos in  "+line);
				continue;
				}
			for(int i=0;i< samples.size();++i)
				{
				try {
					pt.values[i]=Double.parseDouble(tokens[2]);
					if(Double.isNaN(pt.values[i]))
						{
						warning("bad value in "+line);
						continue;
						}
					}
				catch (Exception e) {
					warning("bad value in  "+line);
					continue;//TODO
					}
				samples.get(i).max_value=Math.max(samples.get(i).max_value, pt.values[i]);
				samples.get(i).min_value=Math.min(samples.get(i).min_value, pt.values[i]);
				}
			this.dataPoints.add(pt);
			ci.max_pos=Math.max(ci.max_pos, pt.pos);
			ci.min_pos=Math.min(ci.min_pos, pt.pos);
			}
		CloserUtil.close(r);
		for(Sample sample:this.samples)
			{
			this.max_value=Math.max(sample.max_value, max_value);
			this.min_value=Math.min(sample.min_value, min_value);
			}
		}
	
	@Override
	protected List<ChromInfo> getChromInfos() {
		return this.chromInfos;
		}
	
	/*
	private void readSAM(InputStream in) throws IOException
		{
		info("reading SAM");
		SAMFileReader sfr=new SAMFileReader(in);
		
		SamLocusIterator sli=new SamLocusIterator(sfr);
		sli.setEmitUncoveredLoci(true);
		Iterator<LocusInfo> iter=sli.iterator();
		while(iter.hasNext())
			{
			LocusInfo li=iter.next();
			if(li.getRecordAndPositions().isEmpty()) continue;
			DataPoint pt=new DataPoint();
			pt.tid=li.getSequenceIndex();
			ChromInfo ci=this.chromInfos.get(pt.tid);
			pt.pos=li.getPosition();
			pt.value=li.getRecordAndPositions().size();
			
			this.dataPoints.add(pt);
			
			ci.max_pos=Math.max(ci.max_pos, pt.pos);
			ci.min_pos=Math.min(ci.min_pos, pt.pos);
			this.max_value=Math.max(max_value, pt.value);
			this.min_value=Math.min(min_value, pt.value);
			}
		
		}
		*/
	@Override
	public String getProgramDescription() {
		return "Paint a Genome Scan picture from a Tab delimited file (CHROM/POS/VALUE) or a SAM.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (file.jpg) picture filename out. Required");
		out.println(" -R (fasta) reference fasta indexed with samtools. Optional.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File filout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:R:"))!=-1)
			{
			switch(c)
				{
				case 'o': filout=new File(opt.getOptArg());break;
				case 'R': faidxFile=new File(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		

		if(filout==null)
			{
			error("undefined output filename");
			return -1;
			}
		InputStream in=null;
		try
			{
			
			if(opt.getOptInd()==args.length)
				{
				in=null;
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				in=IOUtils.openURIForReading(filename);
				}
			else
				{
				error("Illegal number of arguments");
				return -1;
				}
			if(faidxFile!=null)
				{
				info("Reading "+faidxFile);
				SAMSequenceDictionary dict=new SAMSequenceDictionaryFactory().load(faidxFile);
				for(SAMSequenceRecord rec: dict.getSequences())
					{
					ChromInfo ci=new ChromInfo();
					ci.sequenceLength=rec.getSequenceLength();
					ci.sequenceName=rec.getSequenceName();
					ci.tid=chromInfos.size();
					chromInfos.add(ci);
					}
				}
			readTsv(in);
			BufferedImage img=makeImage();
			ImageIO.write(img, "JPG", filout);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		
		}
	public static void main(String[] args)
		{
		new GenScan().instanceMainWithExit(args);
		}
	
	}
