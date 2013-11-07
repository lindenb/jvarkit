package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import javax.imageio.ImageIO;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.picard.CigarIterator;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloserUtil;

public class Bam2Raster extends AbstractCommandLineProgram
	{
    
    @Override
    public String getProgramName() {
    	return "BAM to raster graphics.";
    	}
	
   private interface Colorizer
    	{
    	public Color getColor(SAMRecord rec);
    	}
   
   private class QualityColorizer implements Colorizer
		{
	   public Color getColor(SAMRecord rec)
			{	
		    int f=rec.getMappingQuality();
		    if(f>255) f=255;
		    return new Color(f,f,f);
			}
		}
   
   private class FlagColorizer implements Colorizer
		{
	   public Color getColor(SAMRecord rec)
			{
		    if(!rec.getReadPairedFlag() || rec.getProperPairFlag()) return Color.BLACK;
		    if(rec.getMateUnmappedFlag()) return Color.BLUE;
		    if(rec.getDuplicateReadFlag()) return Color.GREEN;
		    return Color.ORANGE;
			}
		}
   
	private File bamFile=null;
	private Interval interval=null;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private int minHDistance=2;
	private int minArrowWidth=2;
	private int maxArrowWidth=5;
	private int featureHeight=30;
	private int spaceYbetweenFeatures=4;
	private Hershey hersheyFont=new Hershey();
	private Colorizer fillColorizer=new QualityColorizer();
	private Colorizer strokeColorizer=new FlagColorizer();
	
	protected double convertToX(int genomic)
		{
		return WIDTH*(genomic-interval.getStart())/(double)(interval.getEnd()-interval.getStart()+1);
		}
	
	protected double left(final SAMRecord rec)
		{
		return convertToX(rec.getAlignmentStart());
		}

	protected double right(final SAMRecord rec)
		{
		return convertToX(rec.getAlignmentEnd());
		}

	
	private BufferedImage build(SAMFileReader r)
		{
		List<List<SAMRecord>> rows=new ArrayList<List<SAMRecord>>();
		SAMRecordIterator iter=null;
		if(interval!=null)
			{
			iter=r.queryOverlapping(interval.getSequence(),interval.getStart(), interval.getEnd());
			}
		else
			{
			iter=r.iterator();
			}
		String currChrom=null;
		int min_left=Integer.MAX_VALUE;
		int max_right=0;
		
		int countReads=0;
		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			
			//when reading from stdin and interval is declared check were ar in the right interval
			if(this.bamFile==null && this.interval!=null)
				{
				if(!this.interval.getSequence().equals(rec.getReferenceName())) continue;
				if(rec.getAlignmentEnd() < this.interval.getStart()) continue;
				if(rec.getAlignmentStart() > this.interval.getEnd()) break;
				}
			
			//when interval is not declared, check only one chromosome
			if(currChrom==null)
				{
				currChrom=rec.getReferenceName();
				}
			else if(!currChrom.equals(rec.getReferenceName()))
				{
				warning("breaking after chromosome "+currChrom);
				break;
				}
			countReads++;
			min_left=Math.min(min_left,rec.getAlignmentStart());
			max_right=Math.max(max_right,rec.getAlignmentEnd());
			
			
			
			for(List<SAMRecord> row:rows)
				{
				SAMRecord last=row.get(row.size()-1);
				if(right(last)+ this.minHDistance > left(rec)) continue;
				row.add(rec);
				rec=null;
				break;
				}
			if(rec!=null)
				{
				List<SAMRecord>  row=new ArrayList<SAMRecord>();
				row.add(rec);
				rows.add(row);
				}
			}
		iter.close();
		
		
		info("Reads:"+countReads);
		final int margin_top=50;
		Dimension imageSize=new Dimension(WIDTH,
				margin_top+ rows.size()*(this.spaceYbetweenFeatures+this.featureHeight)+this.spaceYbetweenFeatures
				);
		BufferedImage img=new BufferedImage(
				imageSize.width,
				imageSize.height,
				BufferedImage.TYPE_INT_RGB
				);
		
		if(interval==null)
			{		
			if(currChrom==null || min_left>max_right) return img;
			interval=new Interval(currChrom, min_left, max_right);
			info("setting interval to "+this.interval);
			}
		
		GenomicSequence genomicSequence=null;
		if(this.indexedFastaSequenceFile !=null)
			{
			genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, this.interval.getSequence());
			}
		
		Graphics2D g=img.createGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, imageSize.width, imageSize.height);
		info("image : "+imageSize.width+"x"+imageSize.height);
		
		//draw bases positions
		
		for(int x=min_left;x<=max_right;++x)
			{
			double oneBaseWidth=convertToX(x+1)-convertToX(x);
			
			g.setColor(x%10==0?Color.BLACK:Color.LIGHT_GRAY);
			g.draw(new Line2D.Double(convertToX(x), 0, convertToX(x), imageSize.height));
			if((x-min_left)%10==0)
				{
				g.setColor(Color.BLACK);
				String xStr=String.valueOf(x);
				AffineTransform tr=g.getTransform();
				AffineTransform tr2=new AffineTransform(tr);
				tr2.translate(convertToX(x), 0);
				tr2.rotate(-Math.PI/2.0);
				hersheyFont.paint(g,
						xStr,
						0,
						0,
						xStr.length()*10,
						oneBaseWidth
						);
				g.setTransform(tr);
				}
			if(genomicSequence!=null)
				{
				char c=genomicSequence.charAt(x-1);
				hersheyFont.paint(g,
						String.valueOf(c),
						convertToX(x)+1,
						1,
						oneBaseWidth-2,
						oneBaseWidth-2
						);
				}
			}
		
		
		int y=margin_top+this.spaceYbetweenFeatures;
		for(List<SAMRecord> row:rows)
			{
			for(SAMRecord rec:row)
				{
				
				double x0=left(rec);
				double x1=right(rec);
				double y0=y;
				double y1=y0+this.featureHeight;
				Shape shapeRec=null;
				if(x1-x0 < minArrowWidth)
					{
					shapeRec=new Rectangle2D.Double(x0, y0, x1-x0, y1-y0);
					}
				else
					{
					GeneralPath path=new GeneralPath();
					double arrow=Math.max(this.minArrowWidth,Math.min(this.maxArrowWidth, x1-x0));
					if(!rec.getReadNegativeStrandFlag())
						{
						path.moveTo(x0, y0);
						path.lineTo(x1-arrow,y0);
						path.lineTo(x1,(y0+y1)/2);
						path.lineTo(x1-arrow,y1);
						path.lineTo(x0,y1);
						}
					else
						{
						path.moveTo(x0+arrow, y0);
						path.lineTo(x0,(y0+y1)/2);
						path.lineTo(x0+arrow,y1);
						path.lineTo(x1,y1);
						path.lineTo(x1,y0);
						}
					path.closePath();
					shapeRec=path;
					}
				
				Stroke oldStroke=g.getStroke();
				g.setStroke(new BasicStroke(2f));
				g.setColor(this.fillColorizer.getColor(rec));
				g.fill(shapeRec);
				g.setColor(this.strokeColorizer.getColor(rec));
				g.draw(shapeRec);
				g.setStroke(oldStroke);
				
				Shape oldClip=g.getClip();
				g.setClip(shapeRec);
				
				if(this.printName)
					{
					hersheyFont.paint(g, rec.getReadName(), shapeRec);
					}
				
				
				CigarIterator ci=CigarIterator.create(rec, this.indexedFastaSequenceFile);
				int prevRefPos=-1;
				while(ci.next())
					{
					int refPos=ci.getReferencePosition();
					if(refPos==-1)
						{
						if(prevRefPos!=-1)
							{
							g.setColor(Color.BLUE); 
							g.fill(new Rectangle2D.Double(
										convertToX(prevRefPos),
										y0,
										2,
										y1-y0
										));
							}
						continue;
						}
					prevRefPos=refPos;
					CigarOperator op=ci.getCigarOperator();
					
					double mutW=convertToX(refPos+1)-convertToX(refPos);
					Shape mut= new Rectangle2D.Double(
							convertToX(refPos),
							y0,
							mutW,
							y1-y0
							);
					
					if(ci.getReadBase()!=null)
						{
						g.setColor(Color.YELLOW);
						this.hersheyFont.paint(g,String.valueOf(ci.getReadBase()),mut);
						}
					
					
					if(ci.isBaseMatching() || op==CigarOperator.EQ)
						{
						continue;
						}
					g.setColor(Color.RED);
					
					switch(op)
						{
						case D: mutW=1; g.setColor(Color.GREEN);break;
						default:break;
						}
					
					
					g.fill(mut);
					
					if(ci.getReadBase()!=null)
						{
						g.setColor(Color.YELLOW);
						this.hersheyFont.paint(g,String.valueOf(ci.getReadBase()),mut);
						}
					
					}
				
				g.setClip(oldClip);
				}
			y+=this.featureHeight+this.spaceYbetweenFeatures;
			}
		g.dispose();
		return img;
		}
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of "+Level.class.getName()+" currently:"+getLogger().getLevel());
		out.println(" -b print bases . Optional. currently:"+printBases);
		out.println(" -r (chr:start-end) restrict to that region. Optional.");
		out.println(" -R (path to fasta) indexed fasta reference. Optional.");
		out.println(" -w (int) image width. Optional. " +WIDTH);
		out.println(" -N print Read name.");
		out.println(" -o (filename) output name. Optional. Default: stdout.");
		}
		
	private boolean printBases=false;
	private boolean printName=false;
	private int WIDTH=1000;
	

	
	@Override
	public int doWork(String args[])
		{
		File fileOut=null;
		File referenceFile=null;
		String region=null;
	    GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, "hvL:o:R:r:w:"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();break;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(Level.parse(getopt.getOptArg()));break;
				case 'o': fileOut=new File(getopt.getOptArg());break;
				case 'R': referenceFile=new File(getopt.getOptArg());break;
				case 'r': region=getopt.getOptArg();break;
				case 'w': this.WIDTH=Math.max(100,Integer.parseInt(getopt.getOptArg()));break;
				case ':': System.err.println("Missing argument for option -"+getopt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+getopt.getOptOpt());return -1;
				}
			}
		
		
		if(getopt.getOptInd()==args.length)
			{
			//stdin
			this.bamFile=null;
			}
		else if(getopt.getOptInd()+1==args.length)
			{
			this.bamFile=new File(args[getopt.getOptInd()]);
			}
		else
			{
			System.err.println("illegal number of arguments.");
			return -1;
			}
	    
		SAMFileReader samFileReader=null;
		try
			{
			if(referenceFile!=null)
				{
				info("loading reference");
				this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(referenceFile);
				}
			
			if(this.bamFile==null)
				{
				warning("READING from stdin");
				samFileReader=new SAMFileReader(System.in);
				}
			else
				{
				info("opening:"+this.bamFile);
				samFileReader=new SAMFileReader(this.bamFile);
				}
			
			SAMFileHeader header=samFileReader.getFileHeader();
			if(region!=null)
				{
				this.interval=IntervalUtils.parseOne(
						header.getSequenceDictionary(),
						region);
				if(this.interval==null)
					{
					System.err.println("Cannot parse interval "+region+" or chrom doesn't exists.");
					return -1;
					}
				}
			else
				{
				warning("NO Interval specified: reading all. Beware memory");
				}
	
			samFileReader.setValidationStringency(ValidationStringency.SILENT);
			BufferedImage img=build(samFileReader);
			samFileReader.close();
			samFileReader=null;
			if(fileOut==null)
				{
				ImageIO.write(img, "PNG", System.out);
				}
			else
				{
				info("saving to "+fileOut);
				ImageIO.write(img, "PNG", fileOut);
				}
			}
		catch(IOException err)
			{
			error(err, err);
			return -1;
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			CloserUtil.close(samFileReader);
			}
		return 0;
		}
	
	public static void main(String[] args)
		{
		new Bam2Raster().instanceMainWithExit(args);
		}
	}
