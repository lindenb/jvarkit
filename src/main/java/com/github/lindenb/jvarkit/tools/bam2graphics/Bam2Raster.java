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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;

import javax.imageio.ImageIO;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
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
    public String getProgramDescription() {
    	return "BAM to raster graphics.";
    	}
	
   private interface Colorizer
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
	
	private BufferedImage build(SAMFileReader r)
		{
		List<List<SAMRecord>> rows=new ArrayList<List<SAMRecord>>();
		SAMRecordIterator iter=null;
		if(bamFile!=null)//got index
			{
			iter=r.queryOverlapping(interval.getSequence(),interval.getStart(), interval.getEnd());
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
				if(!this.interval.getSequence().equals(rec.getReferenceName())) continue;
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
		
	
		
		info("Reads:"+countReads);
		final int ruler_height=String.valueOf(this.interval.getEnd()).length()*20;
		final int refw=(int)Math.max(1.0, WIDTH/(double)(1+interval.getEnd()-interval.getStart()));
		info("refw:"+refw+" "+WIDTH+" "+(1+interval.getEnd()-interval.getStart()));
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
					this.indexedFastaSequenceFile, this.interval.getSequence());
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
		info("image : "+imageSize.width+"x"+imageSize.height);
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
					int refpos=rec.getUnclippedStart();
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
									if(printName)
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
									this.hersheyFont.paint(g,String.valueOf(c1),mut);
									
									readpos++;
									refpos++;
									}
								break;
								}
							default: error("cigar element not handled:"+ce.getOperator());break;
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
	public void printOptions(PrintStream out) {
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of "+Level.class.getName()+" currently:"+getLogger().getLevel());
		out.println(" -b print bases . Optional. currently:"+printBases);
		out.println(" -r (chr:start-end) restrict to that region. REQUIRED.");
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
		while((c=getopt.getopt(args, "hvL:o:R:r:w:N"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(Level.parse(getopt.getOptArg()));break;
				case 'o': fileOut=new File(getopt.getOptArg());break;
				case 'R': referenceFile=new File(getopt.getOptArg());break;
				case 'r': region=getopt.getOptArg();break;
				case 'N': printName=!printName;break;
				case 'w': this.WIDTH=Math.max(100,Integer.parseInt(getopt.getOptArg()));break;
				case ':': System.err.println("Missing argument for option -"+getopt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+getopt.getOptOpt());return -1;
				}
			}
		
		if(region==null)
			{
			error("Region was not defined.");
			return -1;
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
			this.interval=IntervalUtils.parseOne(
					header.getSequenceDictionary(),
					region);
			if(this.interval==null)
				{
				System.err.println(
						"Cannot parse interval "+region+" or chrom doesn't exists in sam dictionary.");
				return -1;
				}
			
	
			samFileReader.setValidationStringency(ValidationStringency.SILENT);
			BufferedImage img=build(samFileReader);
			
			
			
			samFileReader.close();
			samFileReader=null;
			if(img==null)
				{
				error("No image was generated.");
				return -1;
				}
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
