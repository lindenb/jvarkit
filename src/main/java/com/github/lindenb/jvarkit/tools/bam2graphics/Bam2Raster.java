package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import com.github.lindenb.jvarkit.util.picard.CigarIterator;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class Bam2Raster extends CommandLineProgram
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"BAM to raster graphics.";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM files to process.",optional=false)
	public File IN=null;
    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Image name.",optional=false)
	public File OUT=null;

    @Option(shortName= "L", doc="restrict to that region (chr:start-end)",optional=true)
	public String REGION=null;
    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Indexex reference",optional=false)
	public File  REF=null;
    @Option(shortName= "W", doc="image width",optional=true)
    private int WIDTH=1000;
    
	
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
   
	
	private Interval interval=null;
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	private int minHDistance=2;
	private int minArrowWidth=2;
	private int maxArrowWidth=5;
	private int featureHeight=30;
	private int spaceYbetweenFeatures=4;
	
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
		SAMRecordIterator iter=r.queryOverlapping(interval.getSequence(),interval.getStart(), interval.getEnd());

		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
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
		Dimension imageSize=new Dimension(WIDTH,
				rows.size()*(this.spaceYbetweenFeatures+this.featureHeight)+this.spaceYbetweenFeatures
				);
		BufferedImage img=new BufferedImage(
				imageSize.width,
				imageSize.height,
				BufferedImage.TYPE_INT_RGB
				);
		Graphics2D g=img.createGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, imageSize.width, imageSize.height);
		int y=this.spaceYbetweenFeatures;
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
					if(ci.isBaseMatching() || op==CigarOperator.EQ)
						{
						continue;
						}
					g.setColor(Color.RED);
					double mutW=convertToX(refPos+1)-convertToX(refPos);
					switch(op)
						{
						case D: mutW=1; g.setColor(Color.GREEN);break;
						default:break;
						}
					
					Shape mut= new Rectangle2D.Double(
							convertToX(refPos),
							y0,
							mutW,
							y1-y0
							);
					
					g.fill(mut);
					
					}
				
				g.setClip(oldClip);
				}
			y+=this.featureHeight+this.spaceYbetweenFeatures;
			}
		g.dispose();
		return img;
		}
	
	@Override
	protected int doWork()
		{
		SAMFileReader samFileReader=null;
		try
			{
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);

			
			samFileReader=new SAMFileReader(IN);
			
			this.interval=IntervalUtils.parseOne(
					samFileReader.getFileHeader().getSequenceDictionary(),
					REGION);
			if(this.interval==null)
				{
				System.err.println("Cannot parse interval "+REGION+" or chrom doesn't exists.");
				return -1;
				}
	
			samFileReader.setValidationStringency(super.VALIDATION_STRINGENCY);
			BufferedImage img=build(samFileReader);
			samFileReader.close();
			ImageIO.write(img, "PNG", OUT);
			}
		catch(IOException err)
			{
			err.printStackTrace();
			return -1;
			}
		finally
			{
			if(samFileReader!=null) samFileReader.close();
			}
		return 0;
		}
	
	public static void main(String[] args)
		{
		new Bam2Raster().instanceMainWithExit(args);
		}
	}
