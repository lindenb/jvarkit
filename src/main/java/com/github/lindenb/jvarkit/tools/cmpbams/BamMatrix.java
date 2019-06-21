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

*/
package com.github.lindenb.jvarkit.tools.cmpbams;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;

/**
BEGIN_DOC

## Example

```
java -jar dist/bammatrix.jar -o out.png -r "chr1:2345-6789"  -bx NOVASEQ/Sample/outs/phased_possorted_bam.bam
```

![https://twitter.com/yokofakun/status/1038060108373286912](https://pbs.twimg.com/media/Dmft0cSXoAAp78l.jpg)

END_DOC
*/

@Program(name="bammatrix",description="Bam matrix, inspired from 10x/loupe ",
keywords={"sam","bam","compare","matrix"},
creationDate="20190620",
modificationDate="20190621",
generate_doc=false
)
public class BamMatrix  extends Launcher
	{
	private static final Logger LOG = Logger.build(BamMatrix.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx = null;
	@Parameter(names={"-s","--size"},description="matrix size")
	private int matrix_size = 1000;
	@Parameter(names={"-r","-r1","--region"},description="first region." + IntervalParser.OPT_DESC,required=true)
	private String region1Str=null;
	@Parameter(names={"-r2","--region2"},description="2nd region. Default: use first region. " + IntervalParser.OPT_DESC)
	private String region2Str=null;
	@Parameter(names={"--name"},description="use 'BX:Z:' attribute from 10x genomics  as the read name. \"Chromium barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences.\". See https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam")
	private NameExtractor nameExtractor = NameExtractor.READ_NAME;
	@Parameter(names={"-sa","--sa"},description="Use other canonical alignemts from the 'SA:Z:*' attribute")
	private boolean use_sa_align = false;
	@Parameter(names={"--color_-scale"},description="Color scale")
	private ColorScale color_scale = ColorScale.LOG;
	@Parameter(names={"--mapq"},description="min mapping quality")
	private int min_mapq = 30;
	
	/* actual SamReader */
	private SamReader samReader = null;
	/* actual sam dict */
	private SAMSequenceDictionary dict;
	/* user interval X axis */
	private Interval userIntervalX = null;
	/* user interval Y axis */
	private Interval userIntervalY = null;
	
	private enum NameExtractor {
		READ_NAME,BX;
		String getName(final SAMRecord rec) {
			switch(this) {
				case READ_NAME :  return rec.getReadName();
				case BX: return rec.hasAttribute("BX")?String.class.cast(rec.getAttribute("BX")):null;
				default: throw new IllegalStateException(this.name());
				}
			}
		};
	
	private enum ColorScale {
		LINEAR,LOG
	}
	
	
	private abstract class ReadCounter
		{
		protected ReadCounter() {}
		
		boolean accept(final SAMRecord rec) {
			if(rec.getReadUnmappedFlag()) return false;
			if(rec.getReadFailsVendorQualityCheckFlag())  return false;
			if(rec.isSecondaryOrSupplementary())  return false;
			if(rec.getMappingQuality()  < min_mapq) return false;
			return true;
		}
		/** return the names of the Read names in the interval */
		protected abstract Set<String> getNamesMatching(final Interval r) throws IOException;
		
		}
	
	private class MemoryReadCounter extends ReadCounter
		{
		private final IntervalTreeMap<List<Interval>> treeMap  = new IntervalTreeMap<>();
		MemoryReadCounter()  throws IOException{
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samReader.getFileHeader());
			final QueryInterval qIntervalX = new QueryInterval(dict.getSequenceIndex(userIntervalX.getContig()), userIntervalX.getStart(), userIntervalX.getEnd());
			final QueryInterval qIntervalY = new QueryInterval(dict.getSequenceIndex(userIntervalY.getContig()), userIntervalY.getStart(), userIntervalY.getEnd());
			
			final QueryInterval[] qArray = QueryInterval.optimizeIntervals(new QueryInterval[] {qIntervalX,qIntervalY});
			try(final SAMRecordIterator iter= BamMatrix.this.samReader.query(qArray, false))
				{
				while(iter.hasNext()) {
					final SAMRecord rec = iter.next();
					if(!accept(rec)) continue;
					final String name = nameExtractor.getName(rec);
					if(StringUtils.isBlank(name)) continue;
					
					for(final AlignmentBlock ab:rec.getAlignmentBlocks())
						{
						final Interval r = new Interval(
								rec.getReferenceName(),
								ab.getReferenceStart(),
								ab.getReferenceStart()+ab.getLength(),
								rec.getReadNegativeStrandFlag(),
								name
								);
						List<Interval> list = this.treeMap.get(r);
						if(list==null) {
							list=new ArrayList<>();
							this.treeMap.put(r,list);
							}
						list.add(r);
						}
					if( use_sa_align) {
						for(final SAMRecord rec2:SAMUtils.getOtherCanonicalAlignments(rec))
							{
							if(!(userIntervalX.overlaps(rec2) || userIntervalY.overlaps(rec2)))
								{
								continue;
								}
														
							final Interval r = new Interval(
									rec2.getReferenceName(),
									rec2.getStart(),
									rec2.getStart(),
									rec2.getReadNegativeStrandFlag(),
									name
									);
							List<Interval> list = this.treeMap.get(r);
							if(list==null) {
								list=new ArrayList<>();
								this.treeMap.put(r,list);
								}
							list.add(r);
							}
						}
					
					}
				}
			LOG.debug("treeMap.size="+treeMap.size());
			}

		@Override
		protected Set<String> getNamesMatching(final Interval r) throws IOException
			{
			return this.treeMap.getOverlapping(r).
					stream().
					flatMap(L->L.stream()).
					map(R->R.getName()).
					collect(Collectors.toCollection(()->new HashSet<>(10_000)));
			}

		}
	
	@SuppressWarnings("unused")
	private class DiskBackedReadCounter extends ReadCounter {
		@Override
		protected HashSet<String> getNamesMatching(Interval r) throws IOException {
			final HashSet<String> set = new HashSet<>(10_000);
			
			try(final SAMRecordIterator iter=samReader.query(r.getContig(),r.getStart(),r.getEnd(),false))
				{
				while(iter.hasNext()) {
					final SAMRecord rec = iter.next();
					if(!accept(rec)) continue;
					final String name = nameExtractor.getName(rec);
					if(StringUtils.isBlank(name)) continue;
					set.add(name);
					}
				}
			return set;
			}
	}
	
	@Override
	public int doWork(List<String> args) {
		if(StringUtils.isBlank(region2Str)) {
			this.region2Str = region1Str;
		}
		
		try {
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT)
					;
			if(this.faidx!=null) srf.referenceSequence(this.faidx);
			
			this.samReader = srf.open(SamInputResource.of(oneAndOnlyOneFile(args)));
			if(!this.samReader.hasIndex()) {
				LOG.error("Input is not indexed");
				return -1;
				}
			this.dict = SequenceDictionaryUtils.extractRequired(this.samReader.getFileHeader());
			final IntervalParser intervalParser = new IntervalParser(dict);
			this.userIntervalX = intervalParser.parse(this.region1Str);
			if(this.userIntervalX==null) {
				LOG.error("Cannot parse interval "+this.region1Str);
				return -1;
				}
			this.userIntervalY = intervalParser.parse(this.region2Str);
			if(this.userIntervalY==null) {
				LOG.error("Cannot parse interval "+this.region2Str);
				return -1;
				}
			final ReadCounter counter = new MemoryReadCounter();
			
			final int distance= Math.max(
					this.userIntervalX.getLengthOnReference(),
					this.userIntervalY.getLengthOnReference()
					);
			final double pixel2base = distance/(double)matrix_size;
			short max_count=1;
			short counts[]=new short[this.matrix_size*this.matrix_size];

			final ProgressFactory.Watcher<Interval> progress  = ProgressFactory.newInstance().logger(LOG).dictionary(this.dict).build();
			/* loop over each pixel 1st axis */
			for(int pixY=0;pixY< this.matrix_size;pixY++)
				{				
				final int start1 = (int)(this.userIntervalY.getStart() + pixY * pixel2base);
				final int end1 = start1 + (int)pixel2base;
				final Interval qy = new Interval(this.userIntervalY.getContig(), start1, end1);
				progress.apply(qy);
				if(!qy.overlaps(this.userIntervalY)) continue;
				final Set<String> set1 = counter.getNamesMatching(qy);
				if(set1.isEmpty()) continue;
				
				/* loop over each pixel 2nd axis */
				for(int pixX=0;pixX< this.matrix_size;pixX++)
					{
					final int start2 = (int)(this.userIntervalX.getStart() + pixX* pixel2base);
					final int end2 = start2 + (int)pixel2base;
					final Interval qx = new Interval(this.userIntervalX.getContig(), start2, end2);
					if(!qx.overlaps(this.userIntervalX)) continue;
					
					final int count_common;
					if(qx.compareTo(qy)==0) {
						count_common = set1.size();
						}
					else
						{
						final HashSet<String> common = new HashSet<>(set1);
						common.retainAll(counter.getNamesMatching(qx));
						count_common = common.size();
						}
					short count =  count_common>Short.MAX_VALUE?Short.MAX_VALUE:(short)count_common;
					max_count = (short)Math.max(count, max_count);
					counts[pixY*this.matrix_size+pixX] = count;
					}
				}
			progress.close();
			
			final int cov_height = 50;
			final Insets margins = new Insets(10, 100+cov_height, 100+cov_height, 10);
			
			final Dimension drawingAreaDim = new Dimension(
					this.matrix_size+margins.left+margins.right,
					this.matrix_size+margins.top+margins.bottom
					);
			
			final BufferedImage img = new BufferedImage(
					drawingAreaDim.width,
					drawingAreaDim.height,
					BufferedImage.TYPE_INT_RGB
					);
			final Graphics2D g = img.createGraphics();
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, drawingAreaDim.width,drawingAreaDim.height);
			
			g.translate(margins.left, margins.top);
			
			final double logMaxV =Math.log(max_count); 

			for(int pix1=0;pix1< this.matrix_size;pix1++)
				{
				for(int pix2=0;pix2< this.matrix_size;pix2++)
					{
					short count = counts[pix1*this.matrix_size+pix2];
					if(count==0) continue;
					final int gray;
					switch(color_scale) {
						case LINEAR:
							gray = 255-(int)(255*(count/(double)max_count));
							break;
						case LOG:
							gray = 255-(int)(255*((Math.log(count))/logMaxV));
							break;
						default: throw new IllegalStateException(color_scale.name());
						}
					g.setColor(new Color(gray,gray,gray));
					g.fillRect(pix1, pix2, 1, 1);
					}
				}
			// draw frame
			g.setColor(Color.GRAY);
			g.drawRect(0, 0, this.matrix_size, this.matrix_size);

			g.translate(-margins.left,-margins.top);
			
			
			for(int side=0;side< 2;++side) {
				System.err.println(side);
				final Interval r = (side==0?this.userIntervalX:this.userIntervalY);
				final AffineTransform oldtr = g.getTransform();
				AffineTransform tr;
				if(side==1) {
					tr = AffineTransform.getTranslateInstance(0, margins.top);
					tr.concatenate(AffineTransform.getRotateInstance(Math.PI/2.0));
				} else {
					tr = AffineTransform.getTranslateInstance(margins.left, margins.top + matrix_size);
				}
				g.setTransform(tr);
				// draw depth
				final int cov[]=new int[matrix_size];
				int max_cov=1;
				final IntervalList intervalList = new IntervalList(dict);
				intervalList.add(r);
				final SamLocusIterator sli = new SamLocusIterator(this.samReader,intervalList,true);
				while(sli.hasNext()) {
					final LocusInfo locusInfo = sli.next();
					final int pos = locusInfo.getPosition();
					if(pos < r.getStart() || pos > r.getEnd()) continue;
					
					final int depth = locusInfo.getRecordAndOffsets().size();
					final int array_index = (int)(((pos-r.getStart())/(double)r.getLengthOnReference())*matrix_size);
										
					cov[array_index] = Math.max(cov[array_index],depth);
					max_cov = Math.max(max_cov, depth);
					}
				sli.close();
				
				// draw ruler
				
				final GeneralPath gp =new GeneralPath();
				gp.moveTo(0, cov_height);
				for(int x=0;x < cov.length;++x) {
				gp.lineTo(x, cov_height-(cov[x]/(double)max_cov)*cov_height);	
				}
				gp.lineTo(cov.length, cov_height);
				gp.closePath();
				g.setColor(Color.GREEN);
				g.fill(gp);
				
				
				//draw label
				final Hershey herschey = new Hershey();
				g.setColor(Color.DARK_GRAY);
				final int font_size=10;
				String label = StringUtils.niceInt(r.getStart());
				herschey.paint(
						g,
						label,
						new Rectangle2D.Double(0, 0, label.length()*font_size, font_size)
						);
				label = StringUtils.niceInt(r.getEnd());
				herschey.paint(
						g,
						label,
						new Rectangle2D.Double(matrix_size-(label.length()*font_size),0, label.length()*font_size, font_size)
						);
				
				label = r.getContig();
				herschey.paint(
						g,
						label,
						new Rectangle2D.Double(matrix_size/2.0-(label.length()*font_size)/2.0,0, label.length()*font_size, font_size)
						);
				
				g.setTransform(oldtr);
				}
			
			
			g.dispose();
			try {
				if(this.outputFile==null)
					{
					ImageIO.write(img,"PNG",stdout());
					}
				else
					{
					ImageIO.write(img,this.outputFile.getName().endsWith(".png")?"PNG":"JPG", this.outputFile);
					}
				} catch(final IOException err)
				{
					throw new RuntimeIOException(err);
				}
			
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.samReader);
			}
		}
	public static void main(String[] args) {
		new BamMatrix().instanceMainWithExit(args);
		}
}
