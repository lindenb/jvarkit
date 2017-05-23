/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

@Program(name="lowresbam2raster",
		description="Low Resolution BAM to raster graphics",
		keywords={"bam","alignment","graphics","visualization","png","knowngene"}
		)
public class LowResBam2Raster extends Launcher {
	private static final Color ALMOST_BLACK = new Color(20,20,20);
	private static final Color ALMOST_WHITE = new Color(240,240,240);

	private static final Logger LOG = Logger.build(LowResBam2Raster.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--groupby"},description="Group Reads by")
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"-r","--region"},description="restrict to that region. REQUIRED",required=true)
	private String regionStr = null;
	@Parameter(names={"-w","--width"},description="Image width")
	private int WIDTH = 1000 ;
	@Parameter(names={"-clip","--clip"},description="Show clipping")
	private boolean useClip=false;
	@Parameter(names={"--limit","--maxrows"},description="Limit number of rows to 'N' lines. negative: no limit.")
	private int maxRows=-1;
	@Parameter(names={"-depth","--depth"},description="Depth size")
	private int depthSize=100;
	@Parameter(names={"-srf","--samRecordFilter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamFilterParser.buildDefault();
	@Parameter(names={"-kg","--knownGene"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String knownGeneUrl =null;
	private final List<KnownGene> knownGenes = new ArrayList<>();

	private int minDistanceBetweenPairs=3;

	private Interval interval = null;
	private int featureHeight=5;
	private int distnceBetweenFeatures=1;
	private int minArrowWidth=3;
	private final Function<SAMRecord, Integer> readLeft = R->useClip?R.getUnclippedStart():R.getAlignmentStart();
	private final Function<SAMRecord, Integer> readRight = R->useClip?R.getUnclippedEnd():R.getAlignmentEnd();
	private final Function<Integer, Double> pos2pixel = POS->((double)(POS - this.interval.getStart())/(double)this.interval.length()) * this.WIDTH;
	private final Function<Integer, Integer> pixel2pos = PIX-> interval.getStart()+(int)(((double)(PIX)/(double)WIDTH) * (double)this.interval.length());
	private final Function<SAMRecord, Double> left2pixel = R->pos2pixel.apply(readLeft.apply(R));
	private final Function<SAMRecord, Double> right2pixel = R->pos2pixel.apply(readRight.apply(R));
	private final Hershey hersheyFont=new Hershey();
	private final Map<String, PartitionImage> key2partition=new TreeMap<>();
	

	private class SamRecordPair
		implements Locatable
		{
		private SAMRecord R1 = null;
		private SAMRecord R2 = null;
		
		public SamRecordPair(final SAMRecord rec)
			{
			this.R1 = rec;
			}
		
		@Override
		public String getContig() {
			return R1.getContig();
			}
		@Override
		public int getStart() {
			int p = readLeft.apply(this.R1);
			if(R2==null)
				{
				if(this.R1.getReadPairedFlag() && 
					!this.R1.getMateUnmappedFlag() &&
					this.R1.getReferenceName().equals(this.R1.getMateReferenceName()))
					{
					p=Math.min(p, this.R1.getMateAlignmentStart());
					}
				return p;
				}
			else
				{
				return Math.min(p,readLeft.apply(this.R2));
				}
			}
		@Override
		public int getEnd() {
			int p = readRight.apply(this.R1);
			if(R2==null)
				{
				if(this.R1.getReadPairedFlag() && 
					!this.R1.getMateUnmappedFlag() &&
					this.R1.getReferenceName().equals(this.R1.getMateReferenceName()))
					{
					p=Math.max(p, this.R1.getMateAlignmentStart());
					}
				return p;
				}
			else
				{
				return Math.max(p,readRight.apply(this.R2));
				}
			}
		public float getAlpha()
			{
			return 1.f;
			}
		public Color getColor()
			{
			if(samRecordFilter.filterOut(R1)) return Color.PINK;
			if(R1.getReadPairedFlag())
				{
				if(R2==null)
					{
					return Color.RED;
					}
				if(samRecordFilter.filterOut(R2)) return Color.PINK;
				if(!R1.getProperPairFlag())
					{
					return Color.ORANGE;
					}
				}
			return Color.GRAY;
			}
		@Override
		public String toString() {
			final StringBuilder sb=new StringBuilder(R1.getReferenceName()+":"+R1.getStart()+"-"+R1.getEnd());
			if(R2!=null)
				{
				sb.append("/").append(R2.getReferenceName()+":"+R2.getStart()+"-"+R2.getEnd());
				}
			
			return sb.toString();
			}
		}

	
		private class PartitionImage
			{
			private final Map<String,List<SamRecordPair>> readName2pairs = new HashMap<>();
			private final String partitionName;
			private BufferedImage image=null;
			PartitionImage(final String partitionName)
				{
				this.partitionName = partitionName;
				}
			
			void visit(final SAMRecord rec)
				{
				if(rec==null || rec.getReadUnmappedFlag() || !rec.getContig().equals(interval.getContig())) return; 
				String readName = rec.getReadName();
				if(readName.endsWith("/1") || readName.endsWith("/2")) readName=readName.substring(0, readName.length()-1);
				
				
				List<SamRecordPair> rpairList = this.readName2pairs.get(readName);
				
				if( rpairList == null)
					{
					rpairList= new ArrayList<>(1);
					rpairList.add(new SamRecordPair(rec));
					this.readName2pairs.put(readName,rpairList);
					}
				else
					{
					int x=0;
					for(x=0;x < rpairList.size() && rec.getReadPairedFlag() ;++x) {
						final SamRecordPair rp = rpairList.get(x);
						if(rp.R2!=null) continue;
						if(!rp.R1.getReadPairedFlag()) continue;
						if(rp.R1.getFirstOfPairFlag() && rec.getFirstOfPairFlag()) continue;//same R1
						if(rp.R1.getSecondOfPairFlag() && rec.getSecondOfPairFlag()) continue;//same R2
						if(rp.R1.getMateAlignmentStart()!=rec.getAlignmentStart()) continue;
						if(!rp.R1.getMateReferenceName().equals(rec.getReferenceName())) continue;
						
						if(readLeft.apply(rp.R1) < readLeft.apply(rec)) {
							rp.R2= rec;
							}
						else
							{
							rp.R2= rp.R1;
							rp.R1 = rec;
							}
						break;
						}
					if( x == rpairList.size()) {
						rpairList.add(new SamRecordPair(rec));
						}
					}
				}
			void make() {
				
				final List<List<SamRecordPair>> rows = new ArrayList<>();
				for(final SamRecordPair rp : this.readName2pairs.values().
							stream().
							flatMap(L->L.stream()).
							collect(Collectors.toList()))
					{
					if(LowResBam2Raster.this.samRecordFilter.filterOut(rp.R1))
						{
						if(rp.R2==null || LowResBam2Raster.this.samRecordFilter.filterOut(rp.R2)) {
							continue;
							}
						}
					
					int y=0;
					for(y=0;y< rows.size();++y)
						{
						final List<SamRecordPair> row = rows.get(y);
						final SamRecordPair last = row.get(row.size()-1);
						if( last.getEnd() + minDistanceBetweenPairs >= rp.getStart()) continue;
						break;
						}
					if(y==rows.size() && (LowResBam2Raster.this.maxRows<0 || rows.size()<LowResBam2Raster.this.maxRows))
						{
						final List<SamRecordPair> row = new ArrayList<>();
						row.add(rp);
						rows.add(row);
						}
					}
				this.readName2pairs.clear();
				
				//get coverage
				final int coverage[] = new int[LowResBam2Raster.this.interval.length()];
				Arrays.fill(coverage, 0);
				rows.stream().flatMap(R->R.stream()).forEach(RP->{
					for(int R=0;R<2;++R) {
						final SAMRecord rec= (R==0?RP.R1:RP.R2);
						if(rec==null) continue;
						final Cigar cigar = rec.getCigar();
						if(cigar==null || samRecordFilter.filterOut(rec)) continue;							
						int ref = rec.getAlignmentStart();
						for(final CigarElement ce:cigar.getCigarElements())
							{
							final CigarOperator op = ce.getOperator();
							if(op.consumesReferenceBases())
								{
								if(op.consumesReadBases())
									{
									for(int x=0;x<ce.getLength();++x)
										{	
										final int pos=ref+x;
										if(pos<interval.getStart() || pos>interval.getEnd()) continue;
										coverage[ pos - interval.getStart() ]++;
										}
									}
								ref+=ce.getLength();
								}
							}
						}
					});
				
				
				int y=0;
				final String positionFormat="%,d";
				final int ruler_height = String.format(positionFormat,LowResBam2Raster.this.interval.getEnd()).length()*featureHeight;


				final Dimension imageDimension = new Dimension(WIDTH,0);
				
				imageDimension.height += featureHeight*2+distnceBetweenFeatures;//title
				imageDimension.height += featureHeight*2+distnceBetweenFeatures;//interval
				imageDimension.height += ruler_height+distnceBetweenFeatures;//ruler
				
				imageDimension.height += (knownGenes.size()*(featureHeight+distnceBetweenFeatures));
				imageDimension.height += (rows.size()*(featureHeight+distnceBetweenFeatures));
				
				final BufferedImage img=new BufferedImage(
						imageDimension.width,
						imageDimension.height,
						BufferedImage.TYPE_INT_RGB
						);
				final Graphics2D g= img.createGraphics();
				g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
				g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
				g.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);

				g.setColor(ALMOST_WHITE);
				g.fillRect(0, 0, imageDimension.width, imageDimension.height);
					
				//title
				g.setColor(ALMOST_BLACK);
				hersheyFont.paint(g, this.partitionName, new Rectangle2D.Double(1,y,this.partitionName.length()*featureHeight*2,featureHeight*2));
				y+=featureHeight*2+distnceBetweenFeatures;
				//interval
				g.setColor(ALMOST_BLACK);
				hersheyFont.paint(g, interval.getName(), new Rectangle2D.Double(1,y, interval.getName().length()*featureHeight*2,featureHeight*2));
				y+=featureHeight*2+distnceBetweenFeatures;
				//ruler
				final int ruler_shift=interval.length()/10;
				int ruler_start= interval.getStart() - (ruler_shift==0?0:interval.getStart()%ruler_shift);
				while(ruler_start<interval.getEnd())
					{
					
					ruler_start+=ruler_shift;
					g.setColor(Color.GRAY);
					g.draw(new Line2D.Double(
							pos2pixel.apply( ruler_start ),
							y,
							pos2pixel.apply( ruler_start ),
							imageDimension.height*2
							));
					
					final String xStr=String.format(positionFormat,ruler_start);
					final AffineTransform tr=g.getTransform();
					final AffineTransform tr2=new AffineTransform(tr);
					tr2.translate(pos2pixel.apply( ruler_start + 1 ), y);
					tr2.rotate(Math.PI/2.0);
					g.setTransform(tr2);
					g.setColor(ALMOST_BLACK);
					hersheyFont.paint(g,
							xStr,
							0,
							0,
							ruler_height,
							featureHeight*2
							);
					g.setTransform(tr);
					if(ruler_shift==0) break;
					}
				y+=ruler_height+distnceBetweenFeatures;

				
				for(final KnownGene gene:LowResBam2Raster.this.knownGenes) {
					final double cdsHeigh= LowResBam2Raster.this.featureHeight*0.9;
					final double y0 = y;
					final double y1 = y0+cdsHeigh;
					final double midY= (y0+y1)/2.0;
			
					g.setColor(ALMOST_BLACK);
					g.draw(new Line2D.Double(
						pos2pixel.apply(gene.getTxStart()+1),
						midY,
						pos2pixel.apply(gene.getTxEnd()),
						midY)
						);
					
					/* strand symbols */
					for(int pixX=0;
						pixX< LowResBam2Raster.this.WIDTH;
						pixX+=30)
						{
						int pos1= LowResBam2Raster.this.pixel2pos.apply(pixX);
						if(pos1< gene.getTxStart()+1) continue;
						if(pos1> gene.getTxEnd()) break;
						final GeneralPath ticks = new GeneralPath();
						if(gene.isPositiveStrand())
							{
							
							}
						else
							{
							
							}
						g.draw(ticks);
						}
				
					/* exons */
					for(KnownGene.Exon exon:gene.getExons())
						{
						final Rectangle2D exonRect= new Rectangle2D.Double(
							pos2pixel.apply(exon.getStart()+1),
							y0,
							pos2pixel.apply(exon.getEnd())-pos2pixel.apply(exon.getStart()+1),
							cdsHeigh
							);
						g.draw(exonRect);
						g.fill(exonRect);
						}
				
					/* coding line */
					if(!gene.isNonCoding())
						{
						final Rectangle2D cdsRect= new Rectangle2D.Double(
							pos2pixel.apply(gene.getCdsStart()+1),
							y0,
							pos2pixel.apply(gene.getCdsEnd())-pos2pixel.apply(gene.getCdsStart()+1),
							cdsHeigh
							);
						g.draw(cdsRect);
						g.fill(cdsRect);
						}
					y+=featureHeight;
					y+=distnceBetweenFeatures;
					}
				
				
				// print depth
				if(LowResBam2Raster.this.depthSize>0)
					{
					final double depth_array[]=new double[LowResBam2Raster.this.WIDTH];
					Arrays.fill(depth_array, 0);
					for(int x=0;x< depth_array.length;++x)
						{
						final int chromStart = pixel2pos.apply(x);
						final int chromEnd = pixel2pos.apply(x+1);
						double sum=0;
						double count=0;
						for(int chromPos=chromStart;chromPos<chromEnd;++chromPos)
							{
							final int idx = chromPos-interval.getStart();
							if(idx<0 || idx>= coverage.length) continue;
							sum += coverage[idx];
							++count;
							}
						depth_array[x]=sum/count;
						}
					
					final double minDepth = Arrays.stream(depth_array).min().orElse(0);
					final double maxDepth = Arrays.stream(depth_array).max().orElse(minDepth+1);
					
					
					final GeneralPath generalPath = new GeneralPath();
					generalPath.moveTo(0, y+LowResBam2Raster.this.depthSize);
					for(int x=0;x< depth_array.length;++x)
						{
						final double d= depth_array[x];
						final double h= ((d-minDepth)/(maxDepth-minDepth))*LowResBam2Raster.this.depthSize;
						generalPath.lineTo(x, y + LowResBam2Raster.this.depthSize - h);
						}
					generalPath.lineTo(WIDTH, y+LowResBam2Raster.this.depthSize);
					generalPath.closePath();;
					
					g.setColor(Color.BLUE);
					g.fill(generalPath);
					g.setColor(ALMOST_BLACK);
					g.draw(generalPath);
					
					final String label="Depth ["+(int)minDepth+" - "+(int)maxDepth+"]";
					for(int x=0;x<2;++x) {	
						g.setColor(x==0?ALMOST_WHITE:ALMOST_BLACK);
						hersheyFont.paint(g,
								label,
								new Rectangle2D.Double(
									1+x,
									y +x + LowResBam2Raster.this.depthSize-10,
									label.length()*10,
									10
									)
								);
							}
					y+= LowResBam2Raster.this.depthSize;
					y+= LowResBam2Raster.this.distnceBetweenFeatures;
					}

				
				for(final List<SamRecordPair> row: rows)
					{
					for(final SamRecordPair rp: row) {
						final double y0=y;
						final double y1=y0+featureHeight;
						final double ymid = y0+featureHeight/2.0;
						final Composite oldComposite = g.getComposite();
						g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,rp.getAlpha()));
						
						g.setColor(rp.getColor());
						g.draw(new Line2D.Double(
								pos2pixel.apply(rp.getStart()),
								ymid,
								pos2pixel.apply(rp.getEnd()),
								ymid
								));
						// debug : print position
						//hersheyFont.paint(g,rp.toString(),(pos2pixel.apply((int)rp.getEnd())).intValue(),y0,rp.toString().length()*12,featureHeight);
						final Set<Integer> refposOfInsertions = new HashSet<>();

						for(int R=0;R<2;++R)
							{
							final SAMRecord rec = (R==0?rp.R1:rp.R2);
							if(rec==null) continue;
							
							final Cigar cigar = rec.getCigar();
							if(cigar!=null)
								{
								int ref = rec.getAlignmentStart();
								for(final CigarElement ce:cigar.getCigarElements())
									{
									final CigarOperator op = ce.getOperator();
									if(op.equals(CigarOperator.INSERTION))
										{
										refposOfInsertions.add(ref);
										}
									if(op.consumesReferenceBases())
										{
										ref+=ce.getLength();
										}
									}
								}
							double x0 = left2pixel.apply(rec);
							double x1 = right2pixel.apply(rec);
							final Shape shapeRec;
							if(x1-x0 < LowResBam2Raster.this.minArrowWidth)
								{
								shapeRec=new Rectangle2D.Double(x0, y0, x1-x0, y1-y0);
								}
							else
								{
								final GeneralPath path=new GeneralPath();
								double arrow=Math.max(LowResBam2Raster.this.minArrowWidth,Math.min(LowResBam2Raster.this.minArrowWidth, x1-x0));
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
							final Shape oldClip=g.getClip();
							g.setColor(rp.getColor());
							g.setClip(shapeRec);
							g.draw(shapeRec);
							g.fill(shapeRec);
							
							if(cigar!=null && LowResBam2Raster.this.useClip)
								{
								int ref=rec.getAlignmentStart();
								for(final CigarElement ce:cigar.getCigarElements()) {
									final CigarOperator op = ce.getOperator();
									if(op.isClipping()) {
										final Rectangle2D clipRect= new Rectangle2D.Double(
												pos2pixel.apply(ref),
												y0,
												pos2pixel.apply(ref+ce.getLength())-pos2pixel.apply(ref),
												featureHeight
												);
										g.setColor(Color.YELLOW);
										g.fill(clipRect);
										}
									if(op.consumesReferenceBases()) ref+=ce.getLength();
									}
								}
							
							g.setClip(oldClip);
							}
						g.setComposite(oldComposite);
						}
					y+= LowResBam2Raster.this.featureHeight;
					y+= LowResBam2Raster.this.distnceBetweenFeatures;
					}
				
				

				g.dispose();
				this.image = img;
				}
			}
		
		private void scan(final SamReader r) {
			final SAMRecordIterator iter=r.query(
						interval.getContig(),
						interval.getStart(),
						interval.getEnd(),
						false
						);
			while(iter.hasNext())
				{
				final SAMRecord rec=iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				
				if(this.samRecordFilter.filterOut(rec)) 
					{
					//don't dicard now, we need to build pairs of reads
					if(!rec.getReadPairedFlag()) continue;
					if(rec.getMateUnmappedFlag()) continue;
					if(!this.interval.getContig().equals(rec.getMateReferenceName())) continue;
					}
			
				if(!this.interval.getContig().equals(rec.getReferenceName())) continue;
								
				if(this.readLeft.apply(rec) < this.interval.getStart()) 
					{
					continue;
					}
				if(this.readRight.apply(rec) > this.interval.getEnd()) {
					break;
					}
				
				final String group=this.groupBy.apply(rec.getReadGroup());
				if(group==null) continue;
				PartitionImage partition =  this.key2partition.get(group);
				if( partition == null)
					{
					partition=new PartitionImage(group);
					this.key2partition.put(group,partition);
					}
				partition.visit(rec);
				}
			iter.close();
			}

		
		@Override
		public int doWork(final List<String> args) {
				if(this.regionStr==null)
					{
					LOG.error("Region was not defined.");
					return -1;
					}
			
			    if(this.WIDTH<100)
			    	{
			    	LOG.info("adjusting WIDTH to 100");
			    	this.WIDTH=100;
			    	}
				
			    
			    
				SamReader samFileReader=null;
				try
					{
					
					
					
					final SamReaderFactory srf = super.createSamReaderFactory();
					
					

					this.interval = new IntervalParser().parse(this.regionStr);
					if(this.interval==null)
						{
						LOG.error("Cannot parse interval "+regionStr+" or chrom doesn't exists in sam dictionary.");
						return -1;
						}
					LOG.info("Interval is "+this.interval );
					
					if(this.knownGeneUrl!=null)
						{
						IntervalTreeMap<List<KnownGene>> map=KnownGene.loadUriAsIntervalTreeMap(this.knownGeneUrl,
								(KG)->(KG.getContig().equals(this.interval.getContig()) && !(KG.getEnd()<this.interval.getStart() || KG.getStart()+1>this.interval.getEnd()
								)));
						
						this.knownGenes.addAll(map.values().stream().flatMap(L->L.stream()).collect(Collectors.toList()));
						}
					

					for(final String bamFile: IOUtils.unrollFiles(args))
						{
						samFileReader = srf.open(SamInputResource.of(bamFile));
						final SAMFileHeader header=samFileReader.getFileHeader();
						final SAMSequenceDictionary dict=header.getSequenceDictionary();

						
						if(dict==null) {
							LOG.error("no dict in "+bamFile);
							return -1;
							}
						if(dict.getSequence(this.interval.getContig())==null){
							LOG.error("no such chromosome in "+bamFile+" "+this.interval);
							return -1;
							}
						scan(samFileReader);
						samFileReader.close();
						samFileReader=null;
						}
					
					if(this.key2partition.isEmpty())
						{
						LOG.error("No data was found.(no Read-Group specified ?");
						return -1;
						}
			
					this.key2partition.values().stream().forEach(P->P.make());
					
					//assemble everything
					int image_width= this.key2partition.values().stream().mapToInt(P->P.image.getWidth()).max().getAsInt();
					int image_height= this.key2partition.values().stream().mapToInt(P->P.image.getHeight()).sum();
					
					final BufferedImage img= new BufferedImage(image_width, image_height, BufferedImage.TYPE_INT_RGB);
					final Graphics2D g=img.createGraphics();
					g.setRenderingHint(
							RenderingHints.KEY_RENDERING,
							RenderingHints.VALUE_RENDER_QUALITY
							);

					int y=0;
					for(final String key:this.key2partition.keySet()) {
						BufferedImage subImg = this.key2partition.get(key).image;
						g.drawImage(subImg,0,y,null);
						y+=subImg.getHeight();
						}
					g.dispose();


				
					if(this.outputFile==null)
						{
						ImageIO.write(img, "PNG", stdout());
						}
					else
						{
						LOG.info("saving to "+this.outputFile);
						final String format=(this.outputFile.getName().toLowerCase().endsWith(".png")?"PNG":"JPG");
						ImageIO.write(img,format, this.outputFile);
						}
					return RETURN_OK;
					}
				catch(Exception err)
					{
					LOG.error(err);
					return -1;
					}
				finally
					{
					CloserUtil.close(samFileReader);
					}
		
				}
			
		public static void main(String[] args)
			{
			new LowResBam2Raster().instanceMainWithExit(args);
			}
		
}
