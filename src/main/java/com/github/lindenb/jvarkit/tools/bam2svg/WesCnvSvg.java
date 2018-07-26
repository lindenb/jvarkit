/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bam2svg;

import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter.OnNotFound;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Input

input is a set of bam file or a file with the '*.list' suffix containing the path to the bam files.

## Example

```
$ find dir -name "*.bam"  > bam.list
$ java -jar dist/wescnvsvg.jar -R ref.fasta -B cnv.bed bam.list > cnv.svg 
```

END_DOC
 */
@Program(name="wescnvsvg",
description="SVG visualization of bam DEPTH for multiple regions",
keywords={"bam","alignment","graphics","visualization","svg","wes","bed","capture","exome"}
)
public class WesCnvSvg  extends Launcher {
	private static final Logger LOG = Logger.build(WesCnvSvg.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-B","--bed","-b","--capture"},description="BED Capture. Regions to be observed.",required=true)
	private File bedFile = null;
	@Parameter(names={"-R","--ref"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidxFile = null;
	@Parameter(names={"-w","--width"},description="Page width")
	private int drawinAreaWidth = 1000 ;
	@Parameter(names={"-smooth","--smooth"},description="Smoothing DEPTH window size. Negative=don't smooth")
	private int smoothSize = 100 ;


	
	
	private class BamInput implements Closeable
		{
		int index;
		File bamFile;
		SamReader samReader=null;
		String sample;
		double minDepth=0;
		double maxDepth=0;
		@Override
		public void close() throws IOException {
			CloserUtil.close(samReader);
			}
		double getPixelHeight() {
			return 100.0;
		}
		}
	private static class SampleInfo
		{
		double coverage[];
		}
	
	private class CaptureInterval
		implements Locatable
		{
		final QueryInterval queryInterval;
		final List<SampleInfo> sampleInfos = new ArrayList<>(bamInputs.size());
		double pixelx=0.0;
		CaptureInterval(final QueryInterval queryInterval)
			{
			this.queryInterval = queryInterval;
			}
		@Override
		public int getStart() {
			return queryInterval.start;
			}
		@Override
		public int getEnd() {
			return queryInterval.end;
			}
		
		public int getBaseLength() {
			return 1+(getEnd()-getStart());
		}
		
		@Override
		public String getContig() {
			return refDict.getSequence(this.queryInterval.referenceIndex).getSequenceName();
			}
		public double getPixelX1() {
			return pixelx;
			}
		
		public double getPixelWidth(){
			return (this.getBaseLength()/(double)WesCnvSvg.this.countBasesToBeDisplayed)*WesCnvSvg.this.drawinAreaWidth;
		}
		
		String getName() {
			return this.getContig()+":"+format(this.getStart())+"-"+format(this.getEnd());
			}
		}
	
	private final List<CaptureInterval> intervals = new ArrayList<>();
	private final List<BamInput> bamInputs = new ArrayList<>();
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	private SAMSequenceDictionary refDict; 
	private DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private double globalMaxDepth = 0.0;
	private int countBasesToBeDisplayed = 0;

	/** convert double to string */
	private String format(double v)
		{
		return this.decimalFormater.format(v);
		}
	
	@Override
	public int doWork(final List<String> args) {
		XMLStreamWriter w = null;
		BufferedReader r = null;
		FileOutputStream fout=null;
		try
			{
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidxFile);
			this.refDict = this.indexedFastaSequenceFile.getSequenceDictionary();
			if(this.refDict==null || this.refDict.isEmpty()) throw new JvarkitException.FastaDictionaryMissing(this.faidxFile);
			final ContigNameConverter contigNameConverter=ContigNameConverter.fromOneDictionary(this.refDict);
			contigNameConverter.setOnNotFound(OnNotFound.RAISE_EXCEPTION);
			
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT).
					referenceSequence(this.faidxFile);
			for(final File bamFile:IOUtils.unrollFiles2018(args)) {
				final BamInput bi = new BamInput();
				bi.index = this.bamInputs.size();
				bi.bamFile = bamFile;
				bi.samReader = srf.open(bamFile);
				JvarkitException.BamHasIndex.verify(bi.samReader);
				final SAMSequenceDictionary dict2= JvarkitException.BamDictionaryMissing.mustHaveDictionary(bi.samReader);
				SequenceUtil.assertSequenceDictionariesEqual(this.refDict, dict2);
				
				bi.sample = bi.samReader.getFileHeader().getReadGroups().stream().
					map(V->V.getSample()).
					filter(S->!StringUtil.isBlank(S)).findFirst().orElse(bamFile.getName());
				this.bamInputs.add(bi);
				}
			if(this.bamInputs.isEmpty()) {
				LOG.error("no bam input");
				return -1;
			}
			
			
			r = IOUtils.openFileForBufferedReading(this.bedFile);
			final BedLineCodec bedCodec = new BedLineCodec();
			String line;
			List<QueryInterval> listQueryIntervals = new ArrayList<>();
			while((line=r.readLine())!=null)
				{
				if(BedLine.isBedHeader(line)) continue;
				final BedLine bed = bedCodec.decode(line);
				if(bed==null || bed.getStart()>=bed.getEnd()) {
					LOG.warn("Ignoring "+line);
					continue;
				}
				Interval interval = bed.toInterval();
				final int tid = this.refDict.getSequenceIndex(contigNameConverter.apply(interval.getContig()));
				if(tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(interval.getContig(), refDict);
				final QueryInterval qInterval = new QueryInterval(tid, interval.getStart(), interval.getEnd());
				listQueryIntervals.add(qInterval);
				}
			this.intervals.addAll( Arrays.stream(QueryInterval.optimizeIntervals(listQueryIntervals.toArray(new QueryInterval[listQueryIntervals.size()]))).
					map(Q->new CaptureInterval(Q)).
					collect(Collectors.toList()));
			this.countBasesToBeDisplayed = intervals.stream().
					mapToInt(R->1+(R.getEnd()-R.getStart())).
					sum();
			if(this.countBasesToBeDisplayed<1) {
				LOG.error("Nothing to display. BED count==0");
				return -1;
				}
			else
				{
				double x1=0;
				for(int i=0;i< this.intervals.size();++i) {
					final CaptureInterval ci = this.intervals.get(i);
					ci.pixelx = x1;
					x1+=ci.getPixelWidth();
					}
				}
			
			
			for(final CaptureInterval ci:this.intervals) {
				final QueryInterval singletonInterval[]=new QueryInterval[]{ci.queryInterval};
				for(final BamInput bi:this.bamInputs)
					{
					final SampleInfo si = new SampleInfo();
					LOG.info("get cov "+ci.getName()+" for "+bi.bamFile);
					ci.sampleInfos.add(si);
					si.coverage= new double[ci.getBaseLength()];
					Arrays.fill(si.coverage, 0f);
					final SAMRecordIterator iter=bi.samReader.queryOverlapping(singletonInterval);
					while(iter.hasNext())
						{
						final SAMRecord rec = iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						final Cigar cigar=rec.getCigar();
						if(cigar==null) continue;
						int ref1=rec.getAlignmentStart();
						for(final CigarElement ce:cigar) {
							final CigarOperator op = ce.getOperator();
							if(op.consumesReferenceBases())
								{
								if(op.consumesReadBases()){
									for(int x=0;x< ce.getLength();++x){
										final int pos=ref1+x;
										if(pos< ci.getStart()) continue;
										if(pos> ci.getEnd()) break;
										si.coverage[pos-ci.getStart()]++;
										}
									}
								ref1+=ce.getLength();
								}
							}
						}
					iter.close();
					
					if(this.smoothSize>0)
						{
						final double newcov[]=new double[si.coverage.length];
						for(int x=0;x<si.coverage.length;++x) {
							double sum=0;
							int count=0;
							for(int y=-this.smoothSize;y<=this.smoothSize;++y)
								{
								int array_index = x+y;
								if(array_index<0 || array_index>= si.coverage.length) continue;
								sum+=si.coverage[array_index];
								count++;
								}
							newcov[x]=(sum/count);
							}
						System.arraycopy(newcov, 0, si.coverage, 0, newcov.length);
						}
					
					}
				}
			// compute min/max depth for each sample
			for(final BamInput bi:this.bamInputs)
				{
				bi.minDepth = this.intervals.stream().flatMapToDouble(CI->
						DoubleStream.of(CI.sampleInfos.get(bi.index).coverage)
						).min().orElse(0);
				bi.maxDepth = this.intervals.stream().flatMapToDouble(CI->
					DoubleStream.of(CI.sampleInfos.get(bi.index).coverage)
					).max().orElse(1);
				}
			
			this.globalMaxDepth = Math.max(1.0,this.bamInputs.stream().mapToDouble(B->B.maxDepth).max().orElse(0));
			LOG.debug("global max depth "+this.globalMaxDepth);
			
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			if(this.outputFile==null)
				{
				w=xof.createXMLStreamWriter(stdout(), "UTF-8");
				}
			else
				{
				fout = new FileOutputStream(this.outputFile);
				w=xof.createXMLStreamWriter(fout, "UTF-8");
				}
			
			w.writeStartDocument("UTF-8", "1.0");
			
			
			final Dimension dim=new Dimension(this.drawinAreaWidth,0);
			final int bed_header_height=20;
			dim.height += bed_header_height;
			dim.height+=(int)this.bamInputs.stream().
					mapToDouble(B->B.getPixelHeight()).
					sum();
			
			LOG.debug("drawing area: "+dim);
			
			w.writeStartElement("svg");
			w.writeAttribute("width", String.valueOf(dim.width));
			w.writeAttribute("height", String.valueOf(dim.height));
			w.writeDefaultNamespace(SVG.NS);
			w.writeNamespace("xlink", XLINK.NS);
			
			
			w.writeStartElement("style");
			w.writeCharacters(
					"g.maing {stroke:black;stroke-width:0.5px;fill:none;}\n"+
					"text.sampleLabel {stroke:none;stroke-width:0.5px;fill:blue;}" +
					"text.captureLabel {stroke:none;stroke-width:0.5px;fill:salmon;text-anchor:middle;}" +
					"polygon.area {stroke:darkgray;stroke-width:0.5px;fill:lightgray;}" +
					"line.linedp {stroke:green;stroke-width:0.5px;}"
					);
			w.writeEndElement();//style
			
			w.writeStartElement("title");
			w.writeCharacters(getProgramName());
			w.writeEndElement();

			w.writeStartElement("g");
			w.writeAttribute("class", "maing");
			
			int y=0;
			w.writeStartElement("g");
			w.writeComment("interval background");
			for(final CaptureInterval ci:this.intervals)
				{
				w.writeStartElement("text");
				w.writeAttribute("class", "captureLabel");
				w.writeAttribute("x",String.valueOf(ci.getPixelX1()+ci.getPixelWidth()/2.0));
				w.writeAttribute("y",String.valueOf(bed_header_height-2));
				w.writeCharacters(ci.getName());
				w.writeEndElement();//text

				
				w.writeStartElement("rect");
					w.writeAttribute("style","stroke:black;fill:none;");
					w.writeAttribute("x",String.valueOf(ci.getPixelX1()));
					w.writeAttribute("y", "0");
					w.writeAttribute("width",String.valueOf(ci.getPixelWidth()));
					w.writeAttribute("height",String.valueOf(dim.height));
					w.writeStartElement("title");
						w.writeCharacters(ci.getName());
					w.writeEndElement();
				w.writeEndElement();
				}
			w.writeEndElement();//interval background
			
			
			y+= bed_header_height;
			
			for(final BamInput bi:this.bamInputs)
				{	
				w.writeComment(bi.bamFile.getPath());
				w.writeStartElement("g");
				w.writeAttribute("transform","translate(0,"+y+")");
				
				w.writeStartElement("text");
				w.writeAttribute("class", "sampleLabel");
				w.writeAttribute("x","5");
				w.writeAttribute("y","12");
				w.writeCharacters(bi.sample);
				w.writeEndElement();//text
				
				for(final CaptureInterval ci:this.intervals)
					{
					final SampleInfo si = ci.sampleInfos.get(bi.index);
					double x=ci.getPixelX1();
					w.writeStartElement("g");
					w.writeAttribute("transform","translate("+x+",0)");
					
					double sum_interval[]=new double[1+(int)ci.getPixelWidth()];
					LOG.info("interval.lenghth : "+sum_interval.length);
					Arrays.fill(sum_interval, 0);
					int count_interval[]=new int[sum_interval.length];
					
					for(int pos=ci.getStart();pos<=ci.getEnd();++pos)
						{
						int array_index=(int)(((pos-ci.getStart())/(double)ci.getBaseLength())*(double)ci.getPixelWidth());
						//LOG.info("pos="+pos+"->"+array_index+"/"+sum_interval.length+" "+ci.getName()+" "+ci.getPixelWidth());
						if(array_index>=sum_interval.length )  throw new IllegalStateException("boum "+array_index);
						
						
						sum_interval[array_index]+=si.coverage[pos-ci.getStart()];
						count_interval[array_index]++;
						}
					for(int z=0;z< count_interval.length;++z)
						{
						if(count_interval[z]==0) 
							{
							sum_interval[z]=0;
							LOG.info("index error  sum==0 for"+z+"/"+sum_interval.length+" "+ci.getName()+" "+ci.getPixelWidth());
							}
						else
							{
							sum_interval[z]/=count_interval[z];
							}
						}
					
					final List<Point2D.Double> points = new ArrayList<>(sum_interval.length);
					points.add(new Point2D.Double(0,bi.getPixelHeight()));

					for(int z=0;z< sum_interval.length;++z)
						{
						double new_y = bi.getPixelHeight()-(sum_interval[z]/this.globalMaxDepth)*bi.getPixelHeight();
						points.add(new Point2D.Double(z,new_y));
						}
					points.add(new Point2D.Double(ci.getPixelWidth(),bi.getPixelHeight()));
					
					
					for(int z=1;z+1<points.size();++z) {
						if(points.get(z).getY()==points.get(z+1).getY())
							{
							points.get(z).x = points.get(z+1).x;
							points.remove(z+1);
							}
						}
				
					points.add(new Point2D.Double(x,bi.getPixelHeight()));
					w.writeEmptyElement("polygon");
					w.writeAttribute("class","area");
					w.writeAttribute("points",
							points.
								stream().
								map(S->format(S.getX())+","+format(S.getY())).
								collect(Collectors.joining(" "))
							);
					w.writeEndElement();//g
					
					int depthshift=10;
					if(this.globalMaxDepth<=10) depthshift=1;
					else if(this.globalMaxDepth>=100) depthshift=100;
					
					int depth=depthshift;
					while(depth< bi.maxDepth)
						{
						double new_y = bi.getPixelHeight()-(depth/this.globalMaxDepth)*bi.getPixelHeight();
						w.writeStartElement("line");
							w.writeAttribute("class","linedp");
							w.writeAttribute("stroke-dasharray","4");
							w.writeAttribute("x1", "0");
							w.writeAttribute("x2", String.valueOf(dim.width));
							w.writeAttribute("y1", String.valueOf(new_y));
							w.writeAttribute("y2", String.valueOf(new_y));
							w.writeStartElement("title");
								w.writeCharacters(String.valueOf(depth));
							w.writeEndElement();
						w.writeEndElement();//line
						depth+=depthshift;
						}
					
					
					}
				// name for this sample
				w.writeStartElement("rect");
					w.writeAttribute("x","5");
					w.writeAttribute("y", "12");
					w.writeCharacters(bi.sample);
				w.writeEndElement();//rect
				
				//frame for this sample
				w.writeStartElement("rect");
					w.writeAttribute("x","0");
					w.writeAttribute("y", "0");
					w.writeAttribute("width",String.valueOf(dim.width));
					w.writeAttribute("height",String.valueOf(bi.getPixelHeight()));
				w.writeEndElement();//rect

				
				w.writeEndElement();//g
				y+=bi.getPixelHeight();
				}
			
			w.writeStartElement("g");
			w.writeComment("interval lines");
			for(int n=0;n<=this.intervals.size();n++)
				{
				w.writeEmptyElement("line");
				String x1=
						n<this.intervals.size()?
						String.valueOf(this.intervals.get(n).getPixelX1()):
						String.valueOf(dim.width)
						;
				w.writeAttribute("x1", x1);
				w.writeAttribute("y1", "0");
				w.writeAttribute("x2", x1);
				w.writeAttribute("y2", String.valueOf(dim.height));
				}
			w.writeEndElement();//interval lines

			
			w.writeEndElement();//g
			w.writeEndElement();//svg
			w.writeEndDocument();
			w.flush();
			w.close();
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(fout);
			CloserUtil.close(r);
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(this.bamInputs);
			}
		}

	
public static void main(final String[] args)
	{
	new WesCnvSvg().instanceMainWithExit(args);
	}
}
