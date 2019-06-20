package com.github.lindenb.jvarkit.tools.cmpbams;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
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
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;

@Program(name="bammatrix",description="Bam matrix ",
keywords={"sam","bam","compare","matix"},
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
	@Parameter(names={"-r","-r1","--region"},description="first region",required=true)
	private String region1Str=null;
	@Parameter(names={"-r2","--region2"},description="2nd region. Default: use first region.")
	private String region2Str=null;
	@Parameter(names={"-bx"},description="use 'BX:Z:' attribute for read name.")
	private boolean use_bx=false;


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
		
		boolean accept(final SAMRecord rec) {
			if(rec.getReadUnmappedFlag()) return false;
			if(rec.getReadFailsVendorQualityCheckFlag())  return false;
			if(rec.isSecondaryOrSupplementary())  return false;
			if(rec.getMappingQuality()  < 30) return false;
			return true;
		}
		
		/** count number of read overlapping both interval */
		public abstract short count(Interval q1,final Interval q2) throws IOException;
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
					for(final SAMRecord rec2:SAMUtils.getOtherCanonicalAlignments(rec))
						{
						//TODO check in qInterval1 & 2
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
			LOG.debug("treeMap.size="+treeMap.size());
			}

		private HashSet<String> getNamesMatching(final Interval r)
			{
			return this.treeMap.getOverlapping(r).
					stream().
					flatMap(L->L.stream()).
					map(R->R.getName()).
					collect(Collectors.toCollection(HashSet::new));
			}

		public short count(final Interval q1,final Interval q2)
			{
			final HashSet<String> set1 = getNamesMatching(q1);
			final HashSet<String> set2 = getNamesMatching(q2);
			set1.retainAll(set2);
			int n= set1.size();
			//LOG.debug("count"+q1+" "+q2+" ="+set1.size()+" "+set2.size()+" "+ n);
			return n>Short.MAX_VALUE?Short.MAX_VALUE:(short)n;
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

			for(int pix1=0;pix1< this.matrix_size;pix1++)
				{
				final int start1 = (int)(r1.getStart() + pix1 * pixel2base);
				final int end1 = start1 + (int)pixel2base;
				final Interval q1 = new Interval(r1.getContig(), start1, end1);
				if(!q1.overlaps(r1)) continue;
				for(int pix2=0;pix2< this.matrix_size;pix2++)
					{
					final int start2 = (int)(r2.getStart() + pix2 * pixel2base);
					final int end2 = start2 + (int)pixel2base;
					final Interval q2 = new Interval(r2.getContig(), start2, end2);
					if(!q2.overlaps(r2)) continue;
					
					short count = counter.count(q1, q2);
					max_count = (short)Math.max(count, max_count);
					counts[pix1*this.matrix_size+pix2] = count;
					}
				}
			final BufferedImage img = new BufferedImage(this.matrix_size, this.matrix_size, BufferedImage.TYPE_INT_RGB);
			final Graphics2D g = img.createGraphics();
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, this.matrix_size, this.matrix_size);
			g.setColor(Color.GRAY);

			final double logMaxV =Math.log(max_count); 

			for(int pix1=0;pix1< this.matrix_size;pix1++)
				{
				for(int pix2=0;pix2< this.matrix_size;pix2++)
					{
					short count = counts[pix1*this.matrix_size+pix2];
					if(count==0) continue;
					final int gray = 255-(int)(255*((Math.log(count))/logMaxV));
					g.setColor(new Color(gray,gray,gray));
					g.fillRect(pix1, pix2, 1, 1);
					}
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
			CloserUtil.close(sr);
			}
		}
	public static void main(String[] args) {
		new BamMatrix().instanceMainWithExit(args);
		}
}
