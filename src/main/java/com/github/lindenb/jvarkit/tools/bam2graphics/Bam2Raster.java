/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.LinearGradientPaint;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

public class Bam2Raster extends AbstractBam2Raster
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Bam2Raster.class);
    public Bam2Raster()
    	{
    	}
  
		
   private static interface Colorizer
    	{
    	public Color getColor(SAMRecord rec);
    	}
   /*
   private class QualityColorizer implements Colorizer
		{
	   public Color getColor(SAMRecord rec)
			{	
		    int f=rec.getMappingQuality();
		    if(f>255) f=255;
		    return new Color(f,f,f);
			}
		}*/
   
   private static class FlagColorizer implements Colorizer
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
	
	private Color base2color(char c)
		{
		switch(Character.toUpperCase(c))
			{
			case 'N': return Color.BLACK;
			case 'A': return Color.RED;
			case 'T': return Color.GREEN;
			case 'G': return Color.YELLOW;
			case 'C': return Color.BLUE;
			default: return Color.ORANGE;
			}
		}
	
	private BufferedImage build(final SamReader r)
		{
		List<List<SAMRecord>> rows=new ArrayList<List<SAMRecord>>();
		SAMRecordIterator iter=null;
		if(bamFile!=null)//got index
			{
			iter=r.queryOverlapping(interval.getContig(),interval.getStart(), interval.getEnd());
			}
		else //loop until we get the data
			{
			iter=r.iterator();
			}
		
		int countReads=0;
		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			
			//when reading from stdin  we're in the right interval
			if(this.bamFile==null)
				{
				if(!this.interval.getContig().equals(rec.getReferenceName())) continue;
				if(rec.getAlignmentEnd() < this.interval.getStart()) continue;
				if(rec.getAlignmentStart() > this.interval.getEnd()) break;
				}
			
			//when interval is not declared, check only one chromosome
			
			countReads++;			
			
			for(List<SAMRecord> row:rows)
				{
				SAMRecord last=row.get(row.size()-1);
				if(this.interval!=null )
					{
					if(right(last)+ this.minHDistance > left(rec)) continue;
					}
				else
					{
					if(last.getAlignmentEnd()+1> rec.getAlignmentStart()) continue;
					}
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
		
	
		
		LOG.info("Reads:"+countReads);
		final int ruler_height=String.valueOf(this.interval.getEnd()).length()*20;
		final int refw=(int)Math.max(1.0, WIDTH/(double)(1+interval.getEnd()-interval.getStart()));
		LOG.info("refw:"+refw+" "+WIDTH+" "+(1+interval.getEnd()-interval.getStart()));
		final int margin_top=10+(refw*2)+ruler_height;
		Dimension imageSize=new Dimension(WIDTH,
				margin_top+ rows.size()*(this.spaceYbetweenFeatures+this.featureHeight)+this.spaceYbetweenFeatures
				);
		BufferedImage img=new BufferedImage(
				imageSize.width,
				imageSize.height,
				BufferedImage.TYPE_INT_RGB
				);
		
		
		
		CharSequence genomicSequence=null;
		if(this.indexedFastaSequenceFile !=null)
			{
			genomicSequence=new GenomicSequence(
					this.indexedFastaSequenceFile, this.interval.getContig());
			}
		else
			{
			genomicSequence=new AbstractCharSequence()
					{
					@Override
					public int length()
						{
						return interval.getStart()+10;
						}
					
					@Override
					public char charAt(int index)
						{
						return 'N';
						}
				};
			}
		Graphics2D g=img.createGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, imageSize.width, imageSize.height);
		LOG.info("image : "+imageSize.width+"x"+imageSize.height);
		Map<Integer, Counter<Character>> ref2consensus=new HashMap<Integer,  Counter<Character>>();
		//draw bases positions
		
		for(int x=this.interval.getStart();x<=this.interval.getEnd();++x)
			{
			final double oneBaseWidth=convertToX(x+1)-convertToX(x);
			//draw vertical line
			g.setColor(x%10==0?Color.BLACK:Color.LIGHT_GRAY);
			g.draw(new Line2D.Double(convertToX(x), 0, convertToX(x), imageSize.height));
			
			if((x-this.interval.getStart())%10==0)
				{
				g.setColor(Color.BLACK);
				String xStr=String.format("%,d",x);
				AffineTransform tr=g.getTransform();
				AffineTransform tr2=new AffineTransform(tr);
				tr2.translate(convertToX(x), 0);
				tr2.rotate(Math.PI/2.0);
				g.setTransform(tr2);
				hersheyFont.paint(g,
						xStr,
						0,
						0,
						ruler_height,
						oneBaseWidth
						);
				g.setTransform(tr);
				}
			
			//paint genomic sequence
			char c=genomicSequence.charAt(x-1);
			g.setColor(base2color(c));
			hersheyFont.paint(g,
					String.valueOf(c),
					convertToX(x)+1,
					ruler_height,
					oneBaseWidth-2,
					oneBaseWidth-2
					);
				
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
				
				Paint oldpaint=g.getPaint();
				LinearGradientPaint gradient=new LinearGradientPaint(
						0f, (float)shapeRec.getBounds2D().getY(),
						0f, (float)shapeRec.getBounds2D().getMaxY(),
						new float[]{0f,0.5f,1f},
						new Color[]{Color.DARK_GRAY,Color.WHITE,Color.DARK_GRAY}
						);
				g.setPaint(gradient);
				g.fill(shapeRec);
				g.setPaint(oldpaint);
				g.setColor(this.strokeColorizer.getColor(rec));
				g.draw(shapeRec);
				g.setStroke(oldStroke);
				
				Shape oldClip=g.getClip();
				g.setClip(shapeRec);
				
				
				Cigar cigar=rec.getCigar();
				if(cigar!=null)
					{
					byte bases[]=rec.getReadBases();
					int refpos=rec.getAlignmentStart();
					int readpos=0;
					for(CigarElement ce:cigar.getCigarElements())
						{
						switch(ce.getOperator())
							{
							case H: break;
							case S: readpos+=ce.getLength();break;
							case I:
								{
								g.setColor(Color.GREEN); 
								g.fill(new Rectangle2D.Double(
											convertToX(refpos),
											y0,
											2,
											y1-y0
											));
								readpos+=ce.getLength();
									
								
								break;
								}
							case D:
							case N:
							case P:
								{
								g.setColor(Color.ORANGE); 
								g.fill(new Rectangle2D.Double(
											convertToX(refpos),
											y0,
											convertToX(refpos+ce.getLength())-convertToX(refpos),
											y1-y0
											));
								
								refpos+=ce.getLength();
								break;
								}
							case EQ:
							case X:
							case M:
								{
								for(int i=0;i< ce.getLength();++i)
									{
									if(readpos>=bases.length)
										{
										System.err.println(rec.getReadName()+" "+rec.getCigarString()+" "+rec.getReadString());
										}
									
									char c1=(char)bases[readpos];
									
									/* handle consensus */
									Counter<Character> consensus=ref2consensus.get(refpos);
									if(consensus==null)
										{
										consensus=new 	Counter<Character>();
										ref2consensus.put(refpos,consensus);
										}
									consensus.incr(Character.toUpperCase(c1));
									
									
									char c2=genomicSequence.charAt(refpos-1);
									
									double mutW=convertToX(refpos+1)-convertToX(refpos);
									g.setColor(Color.BLACK);
									Shape mut= new Rectangle2D.Double(
											convertToX(refpos),
											y0,
											mutW,
											y1-y0
											);
									if(ce.getOperator()==CigarOperator.X ||
										(c2!='N' && c2!='n' && Character.toUpperCase(c1)!=Character.toUpperCase(c2)))
										{
										g.setColor(Color.RED);
										g.fill(mut);
										g.setColor(Color.WHITE);
										}
									
									//print read name instead of base
									if(isPrintName())
										{
										
										if(readpos<rec.getReadName().length())
											{
											c1=rec.getReadName().charAt(readpos);
											c1=rec.getReadNegativeStrandFlag()?
													Character.toLowerCase(c1):Character.toUpperCase(c1);
											}
										else
											{
											c1=' ';
											}
										}
									else if(!super.disablePrintBases)
										{
										c1=' ';
										}
									this.hersheyFont.paint(g,String.valueOf(c1),mut);
									
									readpos++;
									refpos++;
									}
								break;
								}
							default: LOG.error("cigar element not handled:"+ce.getOperator());break;
							}
						}
					}
				
				
				
				g.setClip(oldClip);
				}
			y+=this.featureHeight+this.spaceYbetweenFeatures;
			}
		
		//print consensus
		for(int x=this.interval.getStart();x<=this.interval.getEnd() ;++x)
			{
			Counter<Character> cons=ref2consensus.get(x);
			if(cons==null || cons.getCountCategories()==0)
				{
				continue;
				}
			final double oneBaseWidth=(convertToX(x+1)-convertToX(x))-1;

			double x0=convertToX(x)+1;
			for(Character c:cons.keySetDecreasing())
				{
				
				double weight=oneBaseWidth*(cons.count(c)/(double)cons.getTotal());
				g.setColor(Color.BLACK);
				
				if(genomicSequence!=null &&
					Character.toUpperCase(genomicSequence.charAt(x-1))!=Character.toUpperCase(c))
					{
					g.setColor(Color.RED);
					}
					
			
				hersheyFont.paint(g,
						String.valueOf(c),
						x0,
						ruler_height+refw,
						weight,
						oneBaseWidth-2
						);
				x0+=weight;
				}
				
			}
		
		
		g.dispose();
		return img;
		}

		
	
	
	@Override
	public Collection<Throwable> call() throws Exception
			{
			final List<String> args = getInputFiles();
			
			if(getRegion()==null)
				{
				return wrapException("Region was not defined.");
				}
			
			if(args.isEmpty())
				{
				//stdin
				this.bamFile=null;
				}
			else if(args.size()==1)
				{
				this.bamFile=new File(args.get(0));
				}
			else
				{
				return wrapException("illegal number of arguments.");
				}
		    if(getWIDTH()<100)
		    	{
		    	LOG.info("adjusting WIDTH to 100");
		    	setWIDTH(100);
		    	}
			
			SamReader samFileReader=null;
			try
				{
				if(getReferenceFile()!=null)
					{
					LOG.info("loading reference");
					this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(getReferenceFile());
					}
				SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
				if(this.bamFile==null)
					{
					LOG.warn("READING from stdin");
					samFileReader=srf.open(SamInputResource.of(stdin()));
					}
				else
					{
					LOG.info("opening:"+this.bamFile);
					samFileReader=srf.open(this.bamFile);
					}
				
				final SAMFileHeader header=samFileReader.getFileHeader();
				this.interval=IntervalUtils.parseOne(
						header.getSequenceDictionary(),
						region);
				if(this.interval==null)
					{
					return wrapException("Cannot parse interval "+region+" or chrom doesn't exists in sam dictionary.");
					}
				
		
				BufferedImage img=build(samFileReader);
				
				
				samFileReader.close();
				samFileReader=null;
				if(img==null)
					{
					return wrapException("No image was generated.");
					}
				if(getOutputFile()==null)
					{
					ImageIO.write(img, "PNG", stdout());
					}
				else
					{
					LOG.info("saving to "+getOutputFile());
					ImageIO.write(img, "PNG", getOutputFile());
					}
				return Collections.emptyList();
				}
			catch(IOException err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(indexedFastaSequenceFile);
				CloserUtil.close(samFileReader);
				indexedFastaSequenceFile=null;			
				}
	
			}
		
	public static void main(String[] args)
		{
		new Bam2Raster().instanceMainWithExit(args);
		}
	}
