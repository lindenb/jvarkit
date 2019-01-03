/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.BufferedList;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
/**

BEGIN_DOC




### Screenshot

![img](https://pbs.twimg.com/media/B3in9wrIAAElLz8.jpg)


```
$ java -jar distBamCmpCoverage.jar  -o out.png file1.bam file2.bam fileN.bam
```




END_DOC
*/
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;


/**

BEGIN_DOC


### Screenshot

![img](https://pbs.twimg.com/media/B3in9wrIAAElLz8.jpg)

### Example

```
$ java -jar distBamCmpCoverage.jar  -o out.png file1.bam file2.bam fileN.bam
```


END_DOC
*/

@Program(
	name="bamcmpcoverage",
	description="Creates the figure of a comparative view of the depths sample vs sample. Memory consideration: the tool alloc an array of bits which size is: (MIN(maxdepth-mindepth,pixel_width_for_one_sample) * count_samples)^2",
	keywords={"sam","bam","visualization","coverage"}
	)
public class BamCmpCoverage extends Launcher
	{
	private static final Logger LOG = Logger.build(BamCmpCoverage.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-w","--width"},description="image width")
	private int imgageSize = 1000 ;
	@Parameter(names={"-m","--minDepth"},description="min depth")
	private int minDepth = 0 ;
	@Parameter(names={"-M","--maxDepth"},description="max depth")
	private int maxDepth = 1000 ;
	@Parameter(names={"-r","--region"},description="restrict to region")
	private String regionStr = null;
	@Parameter(names={"-b","--bed"},description="restrict to region")
	private File bedFile = null;
	@Parameter(names={"--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition samRecordPartition = SAMRecordPartition.sample;

	
	
	private double sampleWidth=0;
	private double marginWidth=0;
	private BufferedImage image=null;
	private Map<String,Integer> sample2column=new HashMap<>();
	/** restrict to BED */
	private IntervalTreeMap<Boolean> intervals=null;

		
	
	/** delegate  bitmap for a matrix comparing two samples */
	private class BitDepthMatrix
		{
		/** owner bitmap */
		private BitSampleMatrix owner;
		/** index sample x in owner matrix */
		private int sample_x;
		/** index sample y in owner matrix */
		private int sample_y;
		BitDepthMatrix(BitSampleMatrix owner)
			{
			this.owner=owner;
			}
		public boolean getXY(int x,int y)
			{
			int left= sample_x*owner.bitSize +x;
			int top = sample_y*owner.bitSize +y ;
			
			int bit_index= left + top *owner.getWidth();
			if(bit_index> owner.getWidth()*owner.getWidth()) throw new RuntimeException();
			return this.owner.bitSet.get(bit_index);
			}
		private int convertDepthToPix(int depth)
			{	
			return (int)((((double)depth)/(double)(maxDepth-minDepth))*owner.bitSize);
			}

		public void setDepth(int depthx,int depthy)
			{
			if(depthx>=(maxDepth-minDepth)) return;
			if(depthy>=(maxDepth-minDepth)) return;
			int scaledx=convertDepthToPix(depthx);
			int scaledy=convertDepthToPix(depthy);
			int tx = this.sample_x*owner.bitSize + scaledx ;
			int ty =(this.sample_y*owner.bitSize + scaledy) *owner.getWidth(); 
			
			
			int bit_index= (int)(tx+ty);
			if(bit_index> owner.getWidth()*owner.getWidth())
				{
				throw new RuntimeException(
						" sample "+this.sample_x+"/"+this.sample_y+" "+
						"bit index:"+bit_index+" for "+depthx+"/"+depthy+" "+owner.bitSet.size());
				}
			
			this.owner.bitSet.set(bit_index);
			}
		}
	
	/** stores a bitmap for a matrix of 'N' sample */
	private class BitSampleMatrix
		{
		private BitSet bitSet;
		/** number of pixels required for a square in the bitmap */
		private int bitSize;
		/** number of samples */
		private int n_samples;
		/** bitmap */
		private BitDepthMatrix depthMatrix=null;
		
		
		BitSampleMatrix(int n_samples)
			{
			/* max-min coverage */
			int diffDepth=  BamCmpCoverage.this.maxDepth-BamCmpCoverage.this.minDepth;
			LOG.info("diffDepth:"+diffDepth);
			
			/* square size according to final image size */
			int squarePixel = (int)Math.ceil(BamCmpCoverage.this.sampleWidth);
			LOG.info("squarePixel:"+squarePixel+" ("+BamCmpCoverage.this.sampleWidth+")");
			
			/* we use the minimal size */
			this.bitSize= Math.min(diffDepth,squarePixel);
			if(this.bitSize <1) this.bitSize=1;
			LOG.info("bitSize:"+bitSize);
			this.n_samples=n_samples;
			int matrix_size =(n_samples*n_samples)*(this.bitSize*this.bitSize); 
			LOG.info("Alloc memory for biset size:"+matrix_size);
			this.bitSet = new BitSet(matrix_size);
			this.depthMatrix=new BitDepthMatrix(this);
			}
		
		public int getWidth()
			{
			return this.n_samples*this.bitSize;
			}
		
		
		public BitDepthMatrix get(int sample_x,int sample_y)
			{
			this.depthMatrix.sample_x=sample_x;
			this.depthMatrix.sample_y=sample_y;
			return this.depthMatrix;
			}
		}
	
	private class Depth
		{
		int tid=0;
		int pos=0;
		int depths[]=new int[sample2column.size()];
		public String toString()
			{
			return "("+(tid+1)+"):"+pos;
			}
		}

	private void paint(Graphics2D g,BitDepthMatrix matrix)
		{
		final Line2D.Double segment=new Line2D.Double();
		for(int i=0;i< matrix.owner.bitSize;++i)
			{
			for(int j=0;j<  matrix.owner.bitSize;++j)
				{	
				if(!matrix.getXY(i, j)) continue;
				
				double x=this.marginWidth+matrix.sample_x*this.sampleWidth;
				x+= this.sampleWidth*(i)/((double)matrix.owner.bitSize);

				
				double y=this.marginWidth+matrix.sample_y*this.sampleWidth;
				y+= this.sampleWidth*(j)/((double)matrix.owner.bitSize);
				
				segment.x1= x;
				segment.y1= y;
				segment.x2= x+0.01;
				segment.y2= y+0.01;
				g.draw(segment);
				}
			}
		}
	
	private void paint(final BitSampleMatrix g,final Depth depth)
		{
		final int sampleUnit= this.maxDepth-this.minDepth;
		
		for(int i=0;i< depth.depths.length;++i)
			{		
			int di = depth.depths[i]-this.minDepth;
			if(di<0 || di>=sampleUnit) continue;
			
			
			
			for(int j=0;j< depth.depths.length;++j)
				{	
				int dj = depth.depths[j]-this.minDepth;
				if(dj<0 || dj>=sampleUnit) continue;
				
				g.get(i, j).setDepth(di, dj);
				g.get(j, i).setDepth(dj, di);
				
				}
			}
		}
	
	
	private void readBedFile(final File bedFile)
		{
		if(this.intervals==null)
			{
			intervals=new IntervalTreeMap<Boolean>();
			}
		try
			{
			LOG.info("Reading "+bedFile);
			final BedLineCodec bedLineCodec = new BedLineCodec();
			BufferedReader r=IOUtils.openFileForBufferedReading(bedFile);
			String line;
			while((line=r.readLine())!=null)
				{
				final BedLine bedLine = bedLineCodec.decode(line);
				if(bedLine==null) continue;
				this.intervals.put(bedLine.toInterval(),true);
				}
			CloserUtil.close(r);
			}
		catch(final IOException err)
			{
			LOG.error(err);
			throw new RuntimeException(err);	
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		if(outputFile==null)
			{
			LOG.error("output image file not defined");
			return -1;
			}
		
		if(this.imgageSize<1)
			{
			LOG.error("Bad image size:" +this.imgageSize);
			return -1;
			}
		
		if(this.minDepth<0)
			{
			LOG.error("Bad min depth : "+this.minDepth);
			return -1;
			}
		if(this.minDepth>=this.maxDepth)
			{
			LOG.error("Bad min<max depth : "+this.minDepth+"<"+this.maxDepth);
			return 1;
			}
		
		if(this.bedFile!=null)
			{
			readBedFile(this.bedFile);
			}
		
		if(regionStr!=null && this.intervals!=null)
			{
			LOG.error("bed and interval both defined.");
			return -1;
			}
		
		try
			{
			final ConcatSam.Factory concatSamFactory = new ConcatSam.Factory();
			final SamReaderFactory srf= concatSamFactory.getSamReaderFactory();
			srf.disable(SamReaderFactory.Option.EAGERLY_DECODE);
			srf.disable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS);
			srf.disable(SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS);
			
			if(this.regionStr!=null)
				{
				concatSamFactory.addInterval(this.regionStr);
				}
			
			
			
			ConcatSam.ConcatSamIterator concatIter = concatSamFactory.open(args);
			
			
			
			final SAMSequenceDictionary dict=concatIter.getFileHeader().getSequenceDictionary();
			
			final Set<String> samples= 
					concatIter.getFileHeader().
					getReadGroups().
					stream().
					map(RG->this.samRecordPartition.apply(RG,"N/A")).
					collect(Collectors.toSet());
						
			LOG.info("Samples:"+samples.size());
			for(String sample:samples)
				{
				this.sample2column.put(sample, this.sample2column.size());
				}
			
			//create image
			LOG.info("Creating image "+this.imgageSize+"x"+this.imgageSize);
			this.image=new BufferedImage(this.imgageSize, this.imgageSize, BufferedImage.TYPE_INT_RGB);
			Graphics2D g=this.image.createGraphics();
			this.marginWidth=this.imgageSize*0.05;
			double drawingWidth=(this.imgageSize-1)-marginWidth;
			this.sampleWidth=drawingWidth/samples.size();
			//g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, this.imgageSize, this.imgageSize);
			g.setColor(Color.BLACK);
			Hershey hershey =new Hershey();
			for(final String sample_x:samples)
				{
				double labelHeight=marginWidth;
				if(labelHeight>50) labelHeight=50;
				
				g.setColor(Color.BLACK);
				hershey.paint(g,
						sample_x,
						marginWidth + sample2column.get(sample_x)*sampleWidth,
						marginWidth - labelHeight,
						sampleWidth*0.9,
						labelHeight*0.9
						);
				
        		AffineTransform old=g.getTransform();
        		AffineTransform tr= AffineTransform.getTranslateInstance(
        				marginWidth ,
        				marginWidth + sample2column.get(sample_x)*sampleWidth
        				);
        		tr.rotate(Math.PI/2);
        		g.setTransform(tr);
        		hershey.paint(g,
						sample_x,
						0.0,
						0.0,
						sampleWidth*0.9,
						labelHeight*0.9
						);        		//g.drawString(this.tabixFile.getFile().getName(),0,0);
        		g.setTransform(old);
				
				for(String sample_y:samples)
					{
					
					Rectangle2D rect=new Rectangle2D.Double(
							marginWidth + sample2column.get(sample_x)*sampleWidth,
							marginWidth + sample2column.get(sample_y)*sampleWidth,
							sampleWidth,
							sampleWidth
							);
					g.setColor(Color.BLUE);
					g.draw(new Line2D.Double(
							rect.getMinX(),rect.getMinY(),
							rect.getMaxX(),rect.getMaxY())
							);
					g.setColor(Color.BLACK);
					g.draw(rect);
					}
				}
			
			
			//ceate bit-array
			BitSampleMatrix bitMatrix=new BitSampleMatrix( samples.size());

						
			//preivous chrom
			//int prev_tid=-1;
			BufferedList<Depth> depthList=new BufferedList<Depth>();
			g.setColor(Color.BLACK);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict).logger(LOG);
			LOG.info("Scanning bams...");
			while(concatIter.hasNext())
				{
				final SAMRecord rec=progress.watch(concatIter.next());
				if(this.samRecordFilter.filterOut(rec)) continue;
				final String sample= this.samRecordPartition.getPartion(rec, "N/A");
				final int sample_id= this.sample2column.get(sample);
				
				final Cigar cigar= rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				int refPos=rec.getAlignmentStart();
				
				
				/* cleanup front pos */
				while(!depthList.isEmpty())
					{
					final Depth front=depthList.getFirst();
					
					if(front.tid!=rec.getReferenceIndex().intValue() ||
						front.pos < refPos )
						{
						paint(bitMatrix,front);
						depthList.removeFirst();
						continue;
						}
					else
						{
						break;
						}		
					}
				
				for(final CigarElement ce:cigar.getCigarElements())
					{
					final CigarOperator op=ce.getOperator();
					if(!op.consumesReferenceBases()) continue;
					if(op.consumesReadBases())
						{
						for(int i=0;i< ce.getLength();++i)
							{
							Depth depth=null;
							int pos=refPos+i;
							
							//ignore non-overlapping BED
							if(this.intervals!=null &&
								!this.intervals.containsOverlapping(new Interval(rec.getReferenceName(),pos,pos)))
								{
								continue;
								}
							else if(depthList.isEmpty())
								{
								depth=new Depth();
								depth.pos=pos;
								depth.tid=rec.getReferenceIndex();
								depthList.add(depth);
								}
							else if(depthList.getLast().pos< pos)
								{
								Depth prev=depthList.getLast();

								while(prev.pos< pos)
									{
									depth=new Depth();
									depth.pos=prev.pos+1;
									depth.tid=rec.getReferenceIndex();
									depthList.add(depth);
									prev=depth;
									}
								depth=prev;
								}
							
							else
								{
								int lastPos=depthList.get(depthList.size()-1).pos;
								int distance= lastPos-pos;
								int indexInList=(depthList.size()-1)-(distance);
								if(indexInList<0)
									{
									//can appen when BED declared and partially overlap the read
									continue;
									}
								
								depth = depthList.get((depthList.size()-1)-(distance));
								if(depth.pos!=pos)
									{
									LOG.error(" "+pos+" vs "+depth.pos+" "+lastPos);
									return -1;
									}
								}
							depth.depths[sample_id]++;
							}
						}
					refPos+=ce.getLength();
					}
				}
			while(!depthList.isEmpty())
				{
				//paint(g,depthList.remove(0));
				paint(bitMatrix,depthList.remove(0));
				}
			progress.finish();
			concatIter.close();concatIter=null;
			
			//paint bitset

		
			for(int x=0;x< bitMatrix.n_samples;++x)
				{
				for(int y=0;y< bitMatrix.n_samples;++y)
					{
					LOG.info("Painting...("+x+"/"+y+")");
					paint(g, bitMatrix.get(x, y));
					}
				}
			
			g.dispose();
			
			//save file
			LOG.info("saving " + this.outputFile);
			if(this.outputFile.getName().toLowerCase().endsWith(".png"))
				{
				ImageIO.write(this.image, "PNG",this.outputFile);
				}
			else
				{
				ImageIO.write(this.image, "JPG", this.outputFile);
				}
			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		}
			 	
	public static void main(String[] args) {
		new BamCmpCoverage().instanceMainWithExit(args);
		}
	}
