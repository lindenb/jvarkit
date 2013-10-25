package com.github.lindenb.jvarkit.tools.genscan;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
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
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;

import org.broad.tribble.readers.AsciiLineReader;

import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import net.sf.picard.cmdline.CommandLineProgram;
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
public class GenScan extends CommandLineProgram
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Paint a Genome Scan picture from a Tab delimited file (CHROM/POS/VALUE) or a SAM. ";
    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="indexed referencE. ",optional=false)
	public File REF;
    @Option(shortName= "w", doc="image WIDTH",optional=true)
	public int WIDTH=1000;
    @Option(shortName= "h", doc="image HEIGHT",optional=true)
	public int HEIGHT=300;
    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME,doc="Ouput filename",optional=false)
	public File OUT;
    @Option(shortName="DC",doc="Do not display unseen chromosome",optional=false)	
    public boolean DICARD_UNSEEN_CHROM=false;
    @Option(shortName="DB",doc="ignore chromosomes left/right side if there's no data there.",optional=false)	
    public boolean DISCARD_BOUNDS=false;

    
    @Option(shortName="IIS",doc="Input is a BAM/SAM file",optional=false)	
    public boolean INPUT_IS_SAM=false;
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input Files. default:stdin. ",minElements=0)
	public List<File> IN=new ArrayList<File>();
    @Option(shortName= "alpha", doc="Default opacity",optional=true)
	public double OPACITY=0.2;
 
    
    
    private Hershey hershey=new Hershey();
	private Insets insets=new Insets(20,150, 30, 30);
	private double max_value=-Double.MAX_VALUE;
	private double min_value=Double.MAX_VALUE;
	
	private SAMSequenceDictionary dict=null;
	
	private static Log LOG=Log.getInstance(GenScan.class);
	
	private class ChromInfo
		{
		SAMSequenceRecord ssr;
		double x;
		double width;
		int max_pos=Integer.MIN_VALUE;
		int min_pos=Integer.MAX_VALUE;
		
		public double converPosToPixel(int pos)
			{
			if(DISCARD_BOUNDS) pos-=this.min_pos;
			return this.x+(pos)/((double)ssr.getSequenceLength())*this.width;
			}
		
		Point2D.Double convert(DataPoint dpt)
			{
			return new Point2D.Double(
				converPosToPixel(dpt.pos),
				(HEIGHT-insets.bottom)-((dpt.value-min_value)/(max_value-min_value))*(HEIGHT-(insets.bottom+insets.top))
				);
			}
		
		private Rectangle2D.Double getFrame()
			{
			return new Rectangle2D.Double(
					x,insets.top,
					width,
					HEIGHT-(insets.top+insets.bottom)
					);
			}
		
		void draw(Graphics2D g)
			{
			Rectangle2D.Double r=getFrame();
			g.setColor(Color.GRAY);
			g.draw(r);
			
			double optW=Math.min(this.ssr.getSequenceName().length()*10,this.width);
			Rectangle2D.Double titleRec=new Rectangle2D.Double(
					(x+width/2)-optW/2.0,
					insets.top-22,
					optW,
					20);
			
			hershey.paint(g, this.ssr.getSequenceName(), titleRec);
			}
		
		int getSequenceLength()
			{
			if(DISCARD_BOUNDS)
				{
				return max_pos-min_pos;
				}
			else
				{
				return ssr.getSequenceLength();
				}
			}
		}
	
	private static class DataPoint implements Comparable<DataPoint>
		{
		int tid;
		int pos;
		double value;
		@Override
		public int compareTo(DataPoint other)
			{
			int i=this.tid-other.tid;
			if(i!=0) return i;
			i=this.pos-other.pos;
			if(i!=0) return i;
			return (this.value< other.value?-1:this.value==other.value?0:1);
			}
		}
	
	private static class DataPointCodec 
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
			dpt.value=dis.readDouble();
			return dpt;
			}
		
		@Override
		public void encode(DataOutputStream out, DataPoint o)
				throws IOException
			{
			out.writeInt(o.tid);
			out.writeInt(o.pos);
			out.writeDouble(o.value);
			}
		
		@Override
		public DataPointCodec clone()
			{
			return new DataPointCodec();
			}
		}
	
	private List<ChromInfo> chromInfos;
	private SortingCollection<DataPoint> dataPoints=null;
	
	@Override
	protected int doWork()
		{
		CloseableIterator<DataPoint> iter=null;
		IndexedFastaSequenceFile faidx=null;
		try {
			faidx=new IndexedFastaSequenceFile(REF);
			this.dict=faidx.getSequenceDictionary();
			this.chromInfos=new ArrayList<GenScan.ChromInfo>(dict.size());
			for(SAMSequenceRecord rec:this.dict.getSequences())
				{
				ChromInfo ci=new ChromInfo();
				ci.ssr=rec;
				chromInfos.add(ci);
				
				}
			
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
			
			
			BufferedImage img=new BufferedImage(this.WIDTH, this.HEIGHT, BufferedImage.TYPE_INT_RGB);
			
			Graphics2D g=img.createGraphics();
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, this.WIDTH, this.HEIGHT);
			
			for(int i=0;i< this.chromInfos.size();++i)
				{
				g.setColor(i%2==0?Color.WHITE:Color.LIGHT_GRAY);
				g.fill(this.chromInfos.get(i).getFrame());
				}
			
			
			//draw axis y
			g.setColor(Color.BLACK);
			for(int i=0;i<= 10;++i)
				{
				double v=min_value+ ((max_value-min_value)/10.0)*i;
				double y=(HEIGHT-insets.bottom)-
						((v-min_value)/(max_value-min_value))*(HEIGHT-(insets.bottom+insets.top))
						;
				g.setColor(Color.LIGHT_GRAY);
				g.draw(new Line2D.Double(insets.left-5, y, WIDTH-insets.right, y));
				
				Rectangle2D.Double r=new Rectangle2D.Double(0,y-12,insets.left,24);
				g.setColor(Color.BLACK);
				this.hershey.paint(g, String.valueOf(v), r);
				}

			Composite oldComposite=g.getComposite();
			g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,(float)this.OPACITY));
			iter=dataPoints.iterator();
			while(iter.hasNext())
				{
				DataPoint dpt=iter.next();
				ChromInfo ci=this.chromInfos.get(dpt.tid);
			
				Point2D.Double  xy=ci.convert(dpt);
				
				g.setColor(Color.RED);
				g.draw(new Line2D.Double(xy.x-1, xy.y,xy.x+1, xy.y));
				g.draw(new Line2D.Double(xy.x, xy.y-1,xy.x, xy.y+1));
				
				}
			iter.close();
			iter=null;
			g.setComposite(oldComposite);
			
			for(ChromInfo ci:this.chromInfos)
				{
				ci.draw(g);
				}
			
			
			g.dispose();
			ImageIO.write(img, "PNG", OUT);
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
		}
	
	private void read(InputStream in) throws IOException
		{
		if(INPUT_IS_SAM)
			{
			readSAM(in);
			}
		else
			{
			readTsv(in);
			}
		}
	
	private void readTsv(InputStream in) throws IOException
		{
		LOG.info("reading tsv");
		Pattern tab=Pattern.compile("[\t]");
		AsciiLineReader r=new AsciiLineReader(in);
		String line;
		while((line=r.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			String tokens[]=tab.split(line, 4);
			DataPoint pt=new DataPoint();
			pt.tid=dict.getSequenceIndex(tokens[0]);
			ChromInfo ci=this.chromInfos.get(pt.tid);
			if(pt.tid==-1)
				{
				LOG.warn("Unknown chromosome "+line);
				continue;
				}
			try {
				pt.pos=Integer.parseInt(tokens[1]);
				if(pt.pos< 0 || pt.pos > ci.ssr.getSequenceLength())
					{
					LOG.warn("Position out of range in "+line);
					continue;
					}
				}
			catch (Exception e) {
				System.err.println("bad pos in  "+line);
				continue;
				}
			
				try {
					pt.value=Double.parseDouble(tokens[2]);
					if(Double.isNaN(pt.value))
						{
						LOG.warn("bad value in "+line);
						continue;
						}
					}
				catch (Exception e) {
					System.err.println("bad value in  "+line);
					continue;
					}
				
			this.dataPoints.add(pt);
			ci.max_pos=Math.max(ci.max_pos, pt.pos);
			ci.min_pos=Math.min(ci.min_pos, pt.pos);
			this.max_value=Math.max(max_value, pt.value);
			this.min_value=Math.min(min_value, pt.value);
			}
		}
	private void readSAM(InputStream in) throws IOException
		{
		LOG.info("reading SAM");
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
	public static void main(String[] args)
		{
		new GenScan().instanceMainWithExit(args);
		}
	}
