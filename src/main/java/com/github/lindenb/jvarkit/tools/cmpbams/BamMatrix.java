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
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
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
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;

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
modificationDate="20190620",
generate_doc=false
)
public class BamMatrix  extends Launcher
	{
	private static final Logger LOG = Logger.build(CompareBams.class).make();
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
	@Parameter(names={"-bx"},description="use 'BX:Z:' attribute from 10x genomics  as the read name. \"Chromium barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences.\". See https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam")
	private boolean use_bx=false;
	@Parameter(names={"-sa","--sa"},description="Use other canonical alignemts from the 'SA:Z:*' attribute")
	private boolean use_sa_align = false;
	@Parameter(names={"--linear"},description="Use linear colors (default is log)")
	private boolean linear_colors = false;


	private abstract class ReadCounter
		{
		final SamReader sr;
		final QueryInterval qInterval1;
		final QueryInterval qInterval2;
		ReadCounter(final SamReader sr,QueryInterval qInterval1,final QueryInterval qInterval2) {
			this.sr= sr;
			this.qInterval1 = qInterval1;
			this.qInterval2 = qInterval2;
			}
		String getReadName(final SAMRecord rec) {
			if(use_bx) {
				if(rec.hasAttribute("BX")) {
					return  (String)rec.getAttribute("BX");
				}
				else {
					return null;
				}
			} else
			{
				return rec.getReadName();
			}
		}
		
		private int[] getDepth(final Interval r) throws IOException {
			final int cov[] =  new int[r.getLengthOnReference()];
			
			try(final SAMRecordIterator iter=this.sr.query(
					r.getContig(),r.getStart(),r.getEnd(),false))
				{
				while(iter.hasNext()) {
					final SAMRecord rec = iter.next();
					if(!accept(rec)) continue;
					for(final AlignmentBlock ab:rec.getAlignmentBlocks()) {
						for(int len=0;len< ab.getLength();len++)
							{
							final int pos= ab.getReferenceStart()+len;
							if(pos < r.getStart() || pos>r.getEnd()) continue;
							final int array_index1 = (pos-r.getStart());
							cov[array_index1]++;
							}
						}
					}
				}
			return cov;
			}
		
		boolean accept(final SAMRecord rec) {
			if(rec.getReadUnmappedFlag()) return false;
			if(rec.getReadFailsVendorQualityCheckFlag())  return false;
			if(rec.isSecondaryOrSupplementary())  return false;
			if(rec.getMappingQuality()  < 30) return false;
			return true;
		}
		
		protected abstract Set<String> getNamesMatching(final Interval r) throws IOException;
		
		}
	
	private class MemoryReadCounter extends ReadCounter
		{
		private final IntervalTreeMap<List<Interval>> treeMap  = new IntervalTreeMap<>();
		MemoryReadCounter(final SamReader sr,QueryInterval qInterval1,final QueryInterval qInterval2)  throws IOException{
			super(sr,qInterval1,qInterval2);
			final QueryInterval[] qArray = QueryInterval.optimizeIntervals(new QueryInterval[] {qInterval1,qInterval2});
			try(final SAMRecordIterator iter=this.sr.query(qArray, false))
				{
				while(iter.hasNext()) {
					final SAMRecord rec = iter.next();
					if(!accept(rec)) continue;
					final String name = getReadName(rec);
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
							final int tid2 = rec2.getReferenceIndex();
							if(tid2== this.qInterval1.referenceIndex && CoordMath.overlaps(this.qInterval1.start, this.qInterval1.end, rec2.getStart(),rec2.getEnd()))
								{
								//ok
								}
							else if(tid2== this.qInterval2.referenceIndex && CoordMath.overlaps(this.qInterval2.start, this.qInterval2.end, rec2.getStart(),rec2.getEnd()))
								{
								//ok
								}
							else
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

		protected Set<String> getNamesMatching(final Interval r) throws IOException
			{
			return this.treeMap.getOverlapping(r).
					stream().
					flatMap(L->L.stream()).
					map(R->R.getName()).
					collect(Collectors.toSet());
			}

		}
	
	@SuppressWarnings("unused")
	private class DiskBackedReadCounter extends ReadCounter {
		DiskBackedReadCounter(final SamReader sr,QueryInterval qInterval1,final QueryInterval qInterval2)  throws IOException{
			super(sr,qInterval1,qInterval2);
			}
		@Override
		protected HashSet<String> getNamesMatching(Interval r) throws IOException {
			final HashSet<String> set = new HashSet<>(10_000);
			
			try(final SAMRecordIterator iter=this.sr.query(r.getContig(),r.getStart(),r.getEnd(),false))
				{
				while(iter.hasNext()) {
					final SAMRecord rec = iter.next();
					if(!accept(rec)) continue;
					final String name = getReadName(rec);
					if(StringUtils.isBlank(name)) continue;
					set.add(name);
					}
				}
			return set;
			}
	}
	
	@Override
	public int doWork(List<String> args) {
		SamReader sr = null;
		if(StringUtils.isBlank(region2Str)) {
			this.region2Str = region1Str;
		}
		
		try {
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT)
					;
			if(this.faidx!=null) srf.referenceSequence(this.faidx);
			
			sr = srf.open(SamInputResource.of(oneAndOnlyOneFile(args)));
			if(!sr.hasIndex()) {
				LOG.error("Input is not indexed");
				return -1;
				}
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(sr.getFileHeader());
			final IntervalParser intervalParser = new IntervalParser(dict);
			final Interval r1 = intervalParser.parse(this.region1Str);
			if(r1==null) {
				LOG.error("Cannot parse interval "+this.region1Str);
				return -1;
				}
			final Interval r2 = intervalParser.parse(this.region2Str);
			if(r2==null) {
				LOG.error("Cannot parse interval "+this.region2Str);
				return -1;
				}
			final ReadCounter counter = new MemoryReadCounter(
					sr,
					new QueryInterval(dict.getSequenceIndex(r1.getContig()), r1.getStart(), r1.getEnd()), 
					new QueryInterval(dict.getSequenceIndex(r2.getContig()), r2.getStart(), r2.getEnd())
					);
			
			final int distance= Math.max(r1.getLengthOnReference(),r2.getLengthOnReference());
			final double pixel2base = distance/(double)matrix_size;
			short max_count=1;
			short counts[]=new short[this.matrix_size*this.matrix_size];

			/* loop over each pixel 1st axis */
			for(int pix1=0;pix1< this.matrix_size;pix1++)
				{
				final int start1 = (int)(r1.getStart() + pix1 * pixel2base);
				final int end1 = start1 + (int)pixel2base;
				final Interval q1 = new Interval(r1.getContig(), start1, end1);
				if(!q1.overlaps(r1)) continue;
				final Set<String> set1 = counter.getNamesMatching(r1);
				if(set1.isEmpty()) continue;
				
				/* loop over each pixel 2nd axis */
				for(int pix2=0;pix2< this.matrix_size;pix2++)
					{
					final int start2 = (int)(r2.getStart() + pix2 * pixel2base);
					final int end2 = start2 + (int)pixel2base;
					final Interval q2 = new Interval(r2.getContig(), start2, end2);
					if(!q2.overlaps(r2)) continue;
					
					final int count_common;
					if(q1.compareTo(q2)==0) {
						count_common = set1.size();
						}
					else
						{
						final HashSet<String> common = new HashSet<>(set1);
						common.retainAll(counter.getNamesMatching(r2));
						count_common = common.size();
						}
					short count =  count_common>Short.MAX_VALUE?Short.MAX_VALUE:(short)count_common;
					max_count = (short)Math.max(count, max_count);
					counts[(this.matrix_size - pix1)*this.matrix_size+pix2] = count;
					}
				}
			final Insets margins = new Insets(10, 100, 100, 10);
			
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
					if( linear_colors) {
						gray = 255-(int)(255*(count/(double)max_count));
					} else {
						gray = 255-(int)(255*((Math.log(count))/logMaxV));
						}
					g.setColor(new Color(gray,gray,gray));
					g.fillRect(pix1, pix2, 1, 1);
					}
				}
			// draw frame
			g.setColor(Color.GRAY);
			g.drawRect(0, 0, this.matrix_size, this.matrix_size);

			/*
			final Hershey herschey = new Hershey();
			herschey.paint(
					g,
					r1.getContig()+":"+StringUtils.niceInt(r1.getStart())+"-"+StringUtils.niceInt(r1.getEnd()),
					new Rectangle2D.Double(0, this.matrix_size, 10, 10)
					);
			*/
			
			g.translate(-margins.left,-margins.top);
			
			
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
			CloserUtil.close(sr);
			}
		}
	public static void main(String[] args) {
		new BamMatrix().instanceMainWithExit(args);
		}
}
