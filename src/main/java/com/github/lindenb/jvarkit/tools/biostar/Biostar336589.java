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
package com.github.lindenb.jvarkit.tools.biostar;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
/**
BEGIN_DOC

## input

input is a BED file (file or stdin). https://genome.ucsc.edu/FAQ/FAQformat.html#format1

  * column 1: chrom
  * column 2: start (you'll get faster results if the input is sorted on chrom/start )
  * column 3: end
  * column 4 is the name of the feature
  * column 5 is the score [0-1000] or '.'
  * column 6 strand +/-
  * column 7 ignored
  * column 8 ignored
  * column 9 is '.' or R,G,B (as in the bed specification) or it's treated as a full svg:style (e.g: `fill:red;stroke:blue;` ) 


multiple bed files are splitted into 'tracks'.

## Example


https://gist.github.com/lindenb/b6debad569dcb5112e76da893d68dd81

```
$ wget -O - -q  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" |\
	gunzip -c | awk '{printf("%s\t%s\t%s\t%s\t%d\t+\t.\t.\t%s\n",$2,$3,$4,$8,rand()*1000,NR%20==0?"255,0,250":".");}' |\
	java -jar dist/biostar336589.jar -R src/test/resources/human_b37.dict   --url 'http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__' --title gaps -mr 300 -fh 20 > ~/jeter.svg 
```

```
$ wget -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" |\
	gunzip -c | cut -f 3,5,6 |\
	sort -t $'\t' -k1,1V -k2,2n |\
	bedtools merge |\
	java -jar dist/biostar336589.jar -md 10000 \
	 	-R src/test/resources/human_b37.dict > out.svg
```

https://gist.github.com/lindenb/5250750014441cc36586dd1f47ed0e37

## Example

```
cat src/test/resources/rotavirus_rf.fa.fai  | awk '{if(NR==1) srand(); printf("%s\t0\t%s\t%s\t%d\t%s\n",$1,$2,$1,rand()*1000,rand()<0.5?rand()<0.5?".":"+":"-");}'   > tmp1.bed
cat src/test/resources/rotavirus_rf.fa.fai  | awk '{if(NR==1) srand(); printf("%s\t0\t%s\t%s\t%d\t%s\n",$1,$2,$1,rand()*1000,rand()<0.5?rand()<0.5?".":"+":"-");}'   > tmp2.bed
cat src/test/resources/rotavirus_rf.fa.fai  | awk '{if(NR==1) srand(); printf("%s\t0\t%s\t%s\t%d\t%s\n",$1,$2,$1,rand()*1000,rand()<0.5?rand()<0.5?".":"+":"-");}'   > tmp3.bed


java -jar dist/biostar336589.jar -R src/test/resources/rotavirus_rf.fa -a 60   tmp1.bed tmp2.bed tmp3.bed > out.svg
```

https://gist.github.com/lindenb/7fd1fa1d3dedcfe38e009387d9f8579c


## Screenshot


https://twitter.com/yokofakun/status/1038060108373286912

![https://twitter.com/yokofakun/status/1038060108373286912](https://pbs.twimg.com/media/Dmft0cSXoAAp78l.jpg)

https://twitter.com/yokofakun/status/1039054467889590272

![https://pbs.twimg.com/media/Dmt2GyvWsAAHfvY.jpg](https://pbs.twimg.com/media/Dmt2GyvWsAAHfvY.jpg)


END_DOC
 */
@Program(name="biostar336589",
description="displays circular map as SVG from BED and REF file",
keywords= {"genome","browser","circular","bed","svg"},
biostars= {336589,367522}
)

public class Biostar336589 extends Launcher{
	private static final Logger LOG = Logger.build(Biostar336589.class).make();

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names="-R",description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx=null;
	@Parameter(names={"-css","--css"},description="custom svg css file")
	private File customCssFile=null;
	@Parameter(names="-md",description="min distance in bp between two features on the same arc.")
	private int min_distance_bp =100;
	@Parameter(names="-mr",description="min internal radius")
	private double min_internal_radius =100;
	@Parameter(names="-fh",description="arc height")
	private double feature_height =10;
	@Parameter(names="-ch",description="contig height")
	private double contig_height =10;
	@Parameter(names="-as",description="arrow size")
	private double arrow_size =10;
	@Parameter(names="-da",description="distance between arcs ")
	private double distance_between_arc =10;
	@Parameter(names="-ms",description="skip chromosome reference length lower than this value. ignore if <=0")
	private int skip_chromosome_size = -1;
	@Parameter(names="-a",description="rotate for 'x' seconds. ignore if <=0")
	private int animation_seconds = -1;
	@Parameter(names={"-hist","--histogram"},description="histogram mode: score of each bed item must be defined. Items must not overlap")
	private boolean histogram_mode = false;
	@Parameter(names={"--score-start-color"},description="When coloring with bed/score this is the color start. "+ColorUtils.Converter.OPT_DESC)
	private String colorScoreStartStr = "white";
	@Parameter(names={"--score-end-color"},description="When coloring with bed/score this is the color end. "+ColorUtils.Converter.OPT_DESC)
	private String colorScoreEndStr = "black";
	@Parameter(names={"--width"},description="Linear SVG size. When defined, this option switches to a linear(!=circular) view")
	private int linear_width = -1;
	@Parameter(names={"--bedpe","-bpe","--interact"},description="Input is UCSC interact format https://genome.ucsc.edu/goldenPath/help/interact.html. I will draw arcs.",hidden=true)
	private boolean input_is_bedpe = false;

	@Parameter(names={"-u","--url","--hyperlink"},description=
			"creates a hyperlink when 'click' in an area. "
			+ "The URL must contains __CHROM__, __START__ and __END__ that will be replaced by their values. "
			+ "IGV : \"http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__\" , "
			+ "UCSC: \"http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__\" "
			)
	private String hyperlinkType = "none";
	@Parameter(names={"--title"},description="document title")
	private String domTitle = Biostar336589.class.getSimpleName();

	private class Arc implements Locatable {
		int tid = -1;
		int start;
		int end;
		String name;
		byte strand=(byte)0;
		int score  =0;
		String cssStyle = null;
		@Override
		public String getContig() {
			return dict.getSequence(this.tid).getSequenceName();
			}
		@Override
		public int getStart() {
			return start;
			}
		@Override
		public int getEnd() {
			return end;
			}
		}
	private class ArcInterval extends Arc {
		int targetTid = -1;
		int targetStart;
		int targetEnd;
		byte targetStrand;
		
		public String getTargetContig() {
			return dict.getSequence(this.targetTid).getSequenceName();
			}
		public int getTargetStart() {
			return targetStart;
			}
		public int getTargetEnd() {
			return targetEnd;
			}
		}
	
	private static class TrackContig {
		final List<List<Arc>> rows = new ArrayList<>();
		TrackContig(int tid) {
			
			}
		}
	
	private class Track {
		private final int id;
		private String name = "";
		private List<TrackContig> contigs ;
		Track(int id) {
			this.id=id;
			this.contigs = new ArrayList<>(Biostar336589.this.dict.size());
			for(int tid=0;tid < Biostar336589.this.dict.size();++tid) {
				this.contigs.add(new TrackContig(tid));
				}
			}
		TrackContig get(int tid) {
			return this.contigs.get(tid);
			}
		int maxRows() {
			return this.contigs.stream().mapToInt(C->C.rows.size()).max().orElse(0);
			}
		}
	
	private SAMSequenceDictionary dict;
	private long reference_length;
	private long tid_to_start[];
	private final ColorUtils colorUtils = new ColorUtils();
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private final DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private Color scoreColorStart = Color.WHITE;
	private Color scoreColorEnd = Color.BLACK;

	
	private String format(double v) {
		return this.decimalFormater.format(v);
	}
	
	private static final Pattern RGB_PATTERN = Pattern.compile("\\d+,\\d+,\\d+");

	private byte toStrand(final String strandStr) {
		if(strandStr.equals("+")) {
			return (byte)1;
			}
		else if(strandStr.equals("-")) {
			return  (byte)-1;
		}else  {
			return (byte)0;
			} 
	}
	
	private String toCss(final String s) {
		if(StringUtil.isBlank(s)) return null;
		if(s.contains(":") || s.contains(";"))
			{
			return s.trim();
			}
		if(!RGB_PATTERN.matcher(s).matches()) return null;
		final String tokens[] = CharSplitter.COMMA.split(s);
		if(tokens.length!=3) return null;
		
		try 
			{
			final int r = Integer.parseInt(tokens[0]);
				final int g = Integer.parseInt(tokens[1]);
				final int b = Integer.parseInt(tokens[2]);
				if(r>=0 && r<256 && g>=0 && g<256 && b>=0 && b<256)
					{
					return "fill:rgb("+r+","+g+","+b+")";
					}
			}
		catch(final NumberFormatException err)
			{
			return null;
			}
		return null;
		}
	
	private Color between(final Color i,final Color j,float f)
		{
		if(f<0) f=0f;
		if(f>1) f=1f;
		final int r= i.getRed()+(int)((j.getRed()-i.getRed())*f);
		final int g= i.getGreen()+(int)((j.getGreen()-i.getGreen())*f);
		final int b= i.getBlue()+(int)((j.getBlue()-i.getBlue())*f);
		return new Color(r,g,b);
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.histogram_mode && this.input_is_bedpe) {
			LOG.error("cannot use both histogram and interact format");
			return -1;
		}
		
		if(this.faidx==null) {
			LOG.error("undefined REF");
			return -1;
		}
		if(this.min_internal_radius < 10) {
			this.min_internal_radius = 10;
		}
		if(this.feature_height < 2) {
			this.feature_height = 2;
		}
		if(this.arrow_size<=0) {
			this.arrow_size=0.01;
		}
		if(this.contig_height<=0) 
			{
			this.contig_height = 2.0;
			}
		
		
		try {
			this.scoreColorStart = this.colorUtils.parse(this.colorScoreStartStr);
			if(this.scoreColorStart==null) this.scoreColorStart = Color.WHITE;
			} 
		catch(final Exception err) 
			{
			this.scoreColorStart = Color.WHITE;
			}
		try {
			this.scoreColorEnd = this.colorUtils.parse(this.colorScoreEndStr);
			if(this.scoreColorEnd==null) this.scoreColorEnd = Color.BLACK;
			} 
		catch(final Exception err) 
			{
			this.scoreColorEnd = Color.BLACK;
			}
		
		int maxScore=0;
		BufferedReader br = null;
		PrintWriter out = null;
		try
			{
			this.dict = SAMSequenceDictionaryExtractor.extractDictionary(this.faidx.toPath());
			if(this.dict==null ) {
				throw new JvarkitException.DictionaryMissing(String.valueOf(this.faidx.toString()));
				}
			if(this.skip_chromosome_size > 0 )
				{
				this.dict = new SAMSequenceDictionary(
					this.dict.getSequences().stream().
					filter(S->S.getSequenceLength()> this.skip_chromosome_size).
					collect(Collectors.toList()));
				}
			if(this.dict.isEmpty()) {
				throw new JvarkitException.DictionaryMissing(String.valueOf(this.faidx.toString()));
				}
			
			this.reference_length = this.dict.getReferenceLength();
			this.tid_to_start = new long[this.dict.size()];
			Arrays.fill(this.tid_to_start,0L);
			
			
			long n = 0;
			for(int i=0;i< dict.size();i++) {
				this.tid_to_start[i] = n;
				n += dict.getSequence(i).getSequenceLength();
				}

			final List<Track> tracks = new ArrayList<>(1+args.size());
			final Set<String> skipped_contigs = new HashSet<>();
			final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(this.dict);
			final BedLineCodec codec = new BedLineCodec();
			for(final String filename:args.isEmpty()?Collections.singletonList((String)null):args)
				{
				final Track currentTrack = new Track(tracks.size());
				if(!StringUtil.isBlank(filename)) currentTrack.name =  filename;
				tracks.add(currentTrack);
				
				br = super.openBufferedReader(filename);
				String line;
				while((line=br.readLine())!=null)
					{
					if(StringUtil.isBlank(line) || BedLine.isBedHeader(line)) continue;
					final BedLine bedLine = codec.decode(line);
					final String newCtg = converter.apply(bedLine.getContig());
					if(StringUtil.isBlank(newCtg)) {
						if(skipped_contigs.add(bedLine.getContig())) {
							LOG.warn("unknown contig "+bedLine.getContig()+". Skipping.");
							}
						continue;
						}
					final SAMSequenceRecord ssr = this.dict.getSequence(newCtg);
					if(ssr==null) continue;
					if(bedLine.getStart() > ssr.getSequenceLength()) continue;
					if(bedLine.getEnd() < 1) continue;
					final Arc arc;
					
					if(!this.input_is_bedpe) {
						arc = new Arc();
						arc.tid = ssr.getSequenceIndex();
						arc.start = Math.max(bedLine.getStart(),0);
						arc.end = Math.min(bedLine.getEnd(),ssr.getSequenceLength());
						arc.strand = toStrand(bedLine.getOrDefault(5,"."));
						}
					else
						{
						final ArcInterval arci = new ArcInterval();
						arc = arci;
						for(int side=0;side<2;++side) {
							final int shift = (side==0?8:13);
							final String ctg2 = bedLine.getOrDefault(shift+0, null);
							if(StringUtil.isBlank(ctg2)) {
								if(skipped_contigs.add(ctg2)) {
									LOG.warn("unknown contig "+(side+1)+"/2 "+bedLine+". Skipping.");
									}
								continue;
								}
							final SAMSequenceRecord ssr2 = this.dict.getSequence(ctg2);
							if(ssr2==null) continue;
							int tid2 = ssr2.getSequenceIndex();
							int start2 = -1;
							int end2 = -1;
							byte strand2 = -1;
							try
								{
								start2 =Math.max(0, Integer.parseInt(bedLine.getOrDefault(shift+1, "")));
								end2 =Math.min(ssr2.getSequenceLength(), Integer.parseInt(bedLine.getOrDefault(shift+2, "")));
								strand2 = toStrand(bedLine.getOrDefault(shift+4, "."));
								}
							catch(final NumberFormatException err2)
								{
								tid2=-1;
								start2=-1;
								end2=-1;
								}
							if(side==0)
								{
								arci.tid = tid2;
								arci.start = start2;
								arci.end = end2;
								arci.strand = strand2;
								}
							else
								{
								arci.targetTid = tid2;
								arci.targetStart = start2;
								arci.targetEnd = end2;
								arci.targetStrand = strand2;
								}
							}
						if(	arci.tid<0 ||
							arci.targetTid<0 || 
							arci.end < arci.start ||
							arci.targetEnd< arci.targetStart
							) {
							LOG.warn("bad interval bed record  "+bedLine+". Skipping.");
							continue;
							}
						}
					
					
					arc.name = bedLine.getOrDefault(3,"");
					
					final String scoreStr = bedLine.getOrDefault(4,"0");
					if(StringUtil.isBlank(scoreStr)|| scoreStr.equals("."))
						{
						if(!this.input_is_bedpe && this.histogram_mode )
							{
							LOG.warn("no score defined for "+line+" in histogram mode. skipping.");
							continue;
							}
						}
					else
						{
						try {
							arc.score= Math.min(1000,Math.max(0,Integer.parseInt(scoreStr)));
							maxScore = Math.max(maxScore, arc.score);
							}
						catch(final NumberFormatException err)
							{
							LOG.warn("bad score for "+line);
							if(this.histogram_mode )
								{
								LOG.warn("skipping.");
								continue;
								}
							arc.score=0;
							}
						}
					
					//color
					arc.cssStyle = toCss(bedLine.getOrDefault(this.input_is_bedpe?7:8,""));
					
					
					final TrackContig currTrackContig = currentTrack.get(arc.tid);
					if(this.histogram_mode) // only one row for histograms
						{
						if(currTrackContig.rows.isEmpty())
							{
							currTrackContig.rows.add(new LinkedList<>());
							}
						currTrackContig.rows.get(0).add(arc);
						}
					else
						{
						int y=0;
						for(y=0;y< currTrackContig.rows.size();++y)
							{
							final List<Arc> row = currTrackContig.rows.get(y);
							if(row.stream().noneMatch(A->A.withinDistanceOf(arc,
									this.histogram_mode?0:this.min_distance_bp
									)))
								{
								row.add(0,arc);//add in front, should be faster if data are sorted
								break;
								}
							
							}
						if(y==currTrackContig.rows.size())
							{
							final List<Arc> row = new LinkedList<>();
							currTrackContig.rows.add(row);
							row.add(arc);	
							}
						}
					}
				
				br.close();
				}
			LOG.info("number of arcs : "+ tracks.stream().mapToInt(T->T.maxRows()).sum());
			

			
			out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			final XMLOutputFactory xof = XMLOutputFactory.newInstance();
			final XMLStreamWriter w = xof.createXMLStreamWriter(out);
			w.writeStartDocument("UTF-8", "1.0");
			w.writeStartElement("svg");
			if(this.linear_width<=0)
				{
				writeCircular(w,tracks,maxScore);
				}
			else
				{
				writeLinear(w,tracks,maxScore);
				}
			w.writeEndElement();//svg
			w.writeEndDocument();
			
			
			w.flush();
			w.close();
			CloserUtil.close(w);
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(br);
			}
		}
	
	private void writeStyle(final XMLStreamWriter w) throws XMLStreamException,IOException {
		w.writeStartElement("style");
		if(this.customCssFile!=null) {
			w.writeCharacters(IOUtil.slurp(this.customCssFile));
			}
		else
			{
			w.writeCharacters(
					"g.maing {stroke:black;stroke-width:0.5px;fill:whitesmoke;font-size:10pt;}\n"+
					".feature {stroke:lightcoral;stroke-width:0.3px;fill:rosybrown;opacity:0.8;pointer-events:all;cursor:crosshair;}"+
					".contig0 {stroke:dimgray;stroke-width:0.8px;fill:gainsboro;}"+
					".contig1 {stroke:dimgray;stroke-width:0.8px;fill:lightgrey;}"+
					".contig {stroke:none;fill:steelblue;}"+
					"circle.track {stroke:lightgray;fill:none;stroke-width:0.8px;}"
					);
			}
		w.writeEndElement();//style		
		}
	
	private void writeScript(final XMLStreamWriter w) throws IOException,XMLStreamException
		{
		w.writeStartElement("script");
		final StringBuilder openBrowserFunction = new StringBuilder(
				"function openGenomeBrowser(contig,chromStart,chromEnd) {\n"
				);
		if( hyperlinkType.contains("__CHROM__") &&
			hyperlinkType.contains("__START__")	&&
			hyperlinkType.contains("__END__") &&
			!hyperlinkType.contains("\"")
			)
			{
			openBrowserFunction.append("var url=\""+this.hyperlinkType+"\".replace(/__CHROM__/g,contig).replace(/__START__/g,chromStart).replace(/__END__/g,chromEnd);\n");
			openBrowserFunction.append("window.open(url,'_blank');\n");

			}
		else
			{
			//nothing
			}
		openBrowserFunction.append("}\n");
		
		w.writeCData(
			openBrowserFunction.toString() +
			"function clicked(evt,contig,chromStart,chromEnd){\n" +
		    "   openGenomeBrowser(contig,chromStart,chromEnd);\n" +
		    "}\n");                
		w.writeEndElement();//script
		}
	private void writeTitle(final XMLStreamWriter w) throws XMLStreamException {
		w.writeStartElement("title");
		w.writeCharacters(this.domTitle);
		w.writeEndElement();
		}
	
	private void writeArcFeatures(final XMLStreamWriter w,final Arc arc,int maxScore) throws XMLStreamException
		{
		final String clickedAttribute = "clicked(evt,\""+arc.getContig()+"\","+arc.getStart()+","+arc.getEnd()+")";

		w.writeAttribute("onclick", clickedAttribute);
		
		w.writeAttribute("class", "feature");
		if(!StringUtil.isBlank(arc.cssStyle))
			{
			w.writeAttribute("style",arc.cssStyle);
			}
		else if(maxScore>0)
			{
			final float f = ((arc.score/(float)maxScore));
			final Color c = between(scoreColorStart,scoreColorEnd,f);
			
			w.writeAttribute("style", "fill:rgb(" +
					c.getRed()+ ","+ c.getGreen()+","+ c.getBlue()+
					")");
			}		
		
		w.writeStartElement("title");
		if(!StringUtil.isBlank(arc.name))
			{
			w.writeCharacters(arc.name);
			}
		else
			{
			w.writeCharacters(arc.getContig()+":"+
					niceIntFormat.format(arc.getStart())+"-"+
					niceIntFormat.format(arc.getEnd())
					);
			}
		w.writeEndElement();
		}

	private void writeLinear(
			final XMLStreamWriter w,
			final List<Track> tracks,
			int maxScore
			)throws XMLStreamException,IOException
		{
		
		final int margin_left= 100;
		final int longest = this.dict.getSequences().stream().mapToInt(SSR->SSR.getSequenceLength()).max().orElse(0);
		if(longest<=0) return;
		
		final Function<Integer, Double> tid2y = (TID)->this.distance_between_arc+(TID*(this.feature_height+this.distance_between_arc));
		final Function<Integer, Double> pos2x = (P)-> (P/((double)longest))*this.linear_width;
		
		final int img_width = this.linear_width  +  margin_left*2 +1;
		final int img_height = (int)((this.dict.size()+1)*(this.feature_height+this.distance_between_arc));
		w.writeAttribute("width", String.valueOf(img_width));
		w.writeAttribute("height", String.valueOf(img_height));
		w.writeDefaultNamespace(SVG.NS);
		w.writeNamespace("xlink", XLINK.NS);
		writeStyle(w);
		writeScript(w);
		writeTitle(w);
		w.writeStartElement("g");
		w.writeAttribute("class", "maing");
		w.writeAttribute("transform", "translate("+margin_left+",0)");
		
		w.writeStartElement("g");
		for(int tid=0;tid< this.dict.size();++tid)
			{
			final SAMSequenceRecord ssr=this.dict.getSequence(tid);
			w.writeStartElement("g");
			
			w.writeStartElement("rect");
			w.writeAttribute("class", "contig"+(tid%2));
			w.writeAttribute("x","0");
			w.writeAttribute("y",format(tid2y.apply(tid)));
			w.writeAttribute("width",format(pos2x.apply(ssr.getSequenceLength())));
			w.writeAttribute("height",format(this.feature_height));
			w.writeEndElement();
			
			w.writeStartElement("text");
			w.writeAttribute("class", "contig");
			w.writeAttribute("style", "text-anchor:end;");
			w.writeAttribute("x",String.valueOf(0-1));
			w.writeAttribute("y",format(tid2y.apply(tid)+this.feature_height));
			w.writeCharacters(ssr.getSequenceName());
			w.writeEndElement();
			
			w.writeEndElement();//g
			}
		w.writeEndElement();//end g contigs
		if(tracks.size()>1) {
			LOG.warn("in linear mode, all tracks are merged");
			}
		for(final Track currentTrack: tracks) {
			w.writeStartElement("g");
			w.writeComment("track "+String.valueOf(currentTrack.name));
			for(final TrackContig trackcontig : currentTrack.contigs)
				{
				for(final List<Arc> arcs: trackcontig.rows)
					{
					if(arcs.size()>1) {
						LOG.warn("in linear mode, all tracks are merged");
						}
					for(final Arc arc: arcs)
						{
						double x1= pos2x.apply(arc.start);
						double x2= pos2x.apply(arc.end);
						final double y1 = tid2y.apply(arc.tid);
						final double y2 = y1 + this.feature_height;
						final double midy = this.feature_height/2.0 + y1 ;
						if(x2<=x1) x2+=0.5;
						if(this.input_is_bedpe)
							{
							final ArcInterval arci = ArcInterval.class.cast(arc);
							
							double tx1= pos2x.apply(arci.targetStart);
							double tx2= pos2x.apply(arci.targetEnd);
							
							
							w.writeStartElement("path");
							final StringBuilder sb =  new StringBuilder();
							if(arci.tid==arci.targetTid)
								{
								//TODO
								}
							else
								{
								final double ty1 = tid2y.apply(arci.targetTid);
								final double ty2 = ty1 + this.feature_height;
								final double tmidy = (ty1+ty2)/2.0 ;
								
								sb.append("M ").append(format(x1)).append(" ").append(format(midy));
								sb.append(" L ").append(format(x2)).append(" ").append(format(midy));
								sb.append(" L ").append(format(tx2)).append(" ").append(format(tmidy));
								sb.append(" L ").append(format(tx1)).append(" ").append(format(tmidy));
								sb.append("Z");
								}
							w.writeAttribute("d",sb.toString());
							}
						else if(arc.strand==0 || (x2-x1)<this.arrow_size*2)
							{
							w.writeStartElement("rect");
							w.writeAttribute("x",format(x1));
							w.writeAttribute("y",format(y1));
							w.writeAttribute("width",format(x2-x1));
							w.writeAttribute("height",format(this.feature_height));
							}
						else if(arc.strand==-1) {
							w.writeStartElement("path");
							final StringBuilder sb =  new StringBuilder();
							sb.append("M ").append(format(x1)).append(" ").append(format(midy));
							sb.append(" L ").append(format(x1+this.arrow_size)).append(" ").append(format(y1));
							sb.append(" L ").append(format(x2)).append(" ").append(format(y1));
							sb.append(" L ").append(format(x2)).append(" ").append(format(y2));
							sb.append(" L ").append(format(x1+this.arrow_size)).append(" ").append(format(y2));
							sb.append("Z");
							w.writeAttribute("d",sb.toString());
							}
						else
							{
							w.writeStartElement("path");
							final StringBuilder sb =  new StringBuilder();
							sb.append("M ").append(format(x1)).append(" ").append(format(y1));
							sb.append(" L ").append(format(x2-this.arrow_size)).append(" ").append(format(y1));
							sb.append(" L ").append(format(x2)).append(" ").append(format(midy));
							sb.append(" L ").append(format(x2-this.arrow_size)).append(" ").append(format(y2));
							sb.append(" L ").append(format(x1)).append(" ").append(format(y2));
							sb.append("Z");
							w.writeAttribute("d",sb.toString());
							}
						writeArcFeatures(w,arc,maxScore);
						
						w.writeEndElement();
						}
					}
				}
			w.writeEndElement();//g of track
			}
		
		w.writeEndElement();//g@maing
		}
	
	
	
	private void writeCircular(
			final XMLStreamWriter w,
			final List<Track> tracks,
			int maxScore
			)throws XMLStreamException,IOException
		{
		final double img_radius = 
				this.min_internal_radius +
				(tracks.stream().mapToInt(T->T.maxRows()).sum()+3)*(this.feature_height+this.distance_between_arc)+
				(this.contig_height+this.distance_between_arc)
				;
		double radius = this.min_internal_radius;

		LOG.info("image radius : "+img_radius);

	
		w.writeAttribute("width", String.valueOf(Math.ceil(img_radius*2)));
		w.writeAttribute("height", String.valueOf(Math.ceil(img_radius*2)));
		w.writeDefaultNamespace(SVG.NS);
		w.writeNamespace("xlink", XLINK.NS);
		
		writeStyle(w);
		writeScript(w);
		writeTitle(w);
		


		
		
		w.writeStartElement("g");
		w.writeAttribute("class", "maing");
		
		if(this.animation_seconds >0) {
			w.writeStartElement("animateTransform");
			w.writeAttribute("attributeName", "transform");
			w.writeAttribute("attributeType", "XML");
			w.writeAttribute("type", "rotate");
			w.writeAttribute("from","0 "+format(img_radius)+" "+format(img_radius));
			w.writeAttribute("to","360 "+format(img_radius)+" "+format(img_radius));
			w.writeAttribute("dur",String.valueOf(this.animation_seconds)+"s");
			w.writeAttribute("repeatCount","indefinite");
			w.writeEndElement();
			}
		
		w.writeStartElement("g");
		w.writeAttribute("transform", "translate("+format(img_radius)+","+format(img_radius)+")");
		w.writeStartElement("g");
		for(final Track currentTrack: tracks) {
			w.writeStartElement("g");
			w.writeComment("track "+String.valueOf(currentTrack.name));
			if(currentTrack.id>0)
				{
				w.writeEmptyElement("circle");
				w.writeAttribute("class","track");
				w.writeAttribute("cx","0");
				w.writeAttribute("cy","0");
				w.writeAttribute("r",format(radius));
				}
			
			for(final TrackContig contigArc : currentTrack.contigs)
				{
				double track_contig_radius = radius;
				for(final List<Arc> row:contigArc.rows)
						{
						w.writeStartElement("g");
						for(final Arc arc: row)
							{
							
							w.writeStartElement("path");
							if(this.input_is_bedpe)
								{
								//TODO
								}
							else if(this.histogram_mode)
								{
								final double height = (arc.score/((double)(maxScore<=0?1:maxScore)))*this.feature_height;
								w.writeAttribute("d", arc(arc.tid,track_contig_radius,
										track_contig_radius+height,
										arc.start,
										arc.end,(byte)0));
								}
							else
								{
								w.writeAttribute("d", arc(arc.tid,track_contig_radius,track_contig_radius+this.feature_height,arc.start,arc.end,arc.strand));
								}
							writeArcFeatures(w,arc,maxScore);
							w.writeEndElement();//path
							
							}
						w.writeEndElement();//g
						track_contig_radius += this.distance_between_arc+this.feature_height;
						}
					}
			w.writeEndElement();//g
			radius += currentTrack.maxRows()*(this.distance_between_arc+this.feature_height);
			}
		
		w.writeEndElement();//g
		w.writeStartElement("g");
		radius+=this.feature_height;
		for(int tid=0;tid< dict.size();++tid)
			{
			final SAMSequenceRecord ssr = this.dict.getSequence(tid);
			w.writeStartElement("path");
			w.writeAttribute("class", "contig"+(tid%2));
			w.writeAttribute("d", arc(tid,radius,radius+this.contig_height,0,ssr.getSequenceLength(),(byte)0));

			w.writeStartElement("title");
			w.writeCharacters(ssr.getSequenceName());
			w.writeEndElement();
			
			w.writeEndElement();//path
			
			w.writeStartElement("text");
			w.writeAttribute("class", "contig");
			w.writeAttribute("x", format(radius+this.distance_between_arc+this.contig_height));
			w.writeAttribute("y", "0");
			w.writeAttribute("transform","rotate("+format(Math.toDegrees(
					refpos2ang(this.tid_to_start[tid]+(long)(ssr.getSequenceLength()/2.0))
					))+")");
			w.writeCharacters(ssr.getSequenceName());
			w.writeEndElement();
			}
		w.writeEndElement();//g
		
		w.writeEndElement();//g
		
		w.writeEndElement();//g		
		}
	
	
	
	private String pointToStr(final Point2D.Double p)
		{
		return format(p.x)+" "+format(p.y);
		}
	
	private double refpos2ang(long pos) {
		return (pos/((double)this.reference_length))*(2.0*Math.PI);
		}
	
	private String arc(
		final int tid,
		final double radius1,
		final double radius2,
		final int start,final int end,
		byte strand
		) {
		final long index_at_start = this.tid_to_start[tid];
		final double mid_radius= (radius1+radius2)/2.0;
		final double r_start = refpos2ang(index_at_start+ start);
		final double r_end = refpos2ang(index_at_start+ end);
		final double arc_length = (r_end-r_start)*Math.PI;
		final double angle_arrow= this.arrow_size / mid_radius;

		final Point2D.Double p1 = polarToCartestian(radius1, r_start);
		final Point2D.Double p2 = polarToCartestian(radius1, r_end);
		final Point2D.Double p3 = polarToCartestian(radius2, r_end);
		final Point2D.Double p4 = polarToCartestian(radius2, r_start);
		final StringBuilder sb = new StringBuilder();
		
		if(this.histogram_mode || strand==0|| 2.0*this.arrow_size<=arc_length)  {
			sb.append("M ").append(pointToStr(p1));
			
			sb.append(" A ").
				append(format(radius1)).
				append(" ").
				append(format(radius1)).
				append(" 0"). //X axis rotation
				append(" 0").// large arc
				append(" 1 ").// sweep flag (positive angle direction)
				append(pointToStr(p2));
			
			sb.append(" L").append(pointToStr(p3));
			
			sb.append(" A ").
				append(format(radius2)).
				append(" ").
				append(format(radius2)).
				append(" 0"). //X axis rotation
				append(" 0").// large arc
				append(" 0 ").// sweep flag (positive angle direction)
				append(pointToStr(p4));
			
			sb.append(" L ").append(pointToStr(p1));
			}
		else if(strand==-1) {
			sb.append("M ").append(pointToStr(polarToCartestian(mid_radius, r_start)));
			sb.append("L ").append(pointToStr(polarToCartestian(radius1, r_start+angle_arrow)));
			
			sb.append(" A ").
				append(format(radius1)).
				append(" ").
				append(format(radius1)).
				append(" 0"). //X axis rotation
				append(" 0").// large arc
				append(" 1 ").// sweep flag (positive angle direction)
				append(pointToStr(p2));

			sb.append("L ").append(pointToStr(p3));	
			
			sb.append(" A ").
				append(format(radius2)).
				append(" ").
				append(format(radius2)).
				append(" 0"). //X axis rotation
				append(" 0").// large arc
				append(" 0 ").// sweep flag (positive angle direction)
				append(pointToStr(polarToCartestian(radius2, r_start+angle_arrow)));
			
			}
		else if(strand==1)
			{			
			sb.append("M ").append(pointToStr(p1));
			
			sb.append(" A ").
				append(format(radius1)).
				append(" ").
				append(format(radius1)).
				append(" 0"). //X axis rotation
				append(" 0").// large arc
				append(" 1 ").// sweep flag (positive angle direction)
				append(pointToStr(polarToCartestian(radius1, r_end - angle_arrow)));
			
			sb.append("L ").append(pointToStr(polarToCartestian(mid_radius, r_end)));
			sb.append("L ").append(pointToStr(polarToCartestian(radius2, r_end-angle_arrow)));	
			
			sb.append(" A ").
				append(format(radius2)).
				append(" ").
				append(format(radius2)).
				append(" 0"). //X axis rotation
				append(" 0").// large arc
				append(" 0 ").// sweep flag (positive angle direction)
				append(pointToStr(p4));
			}
		sb.append(" Z");
		
		return sb.toString();
		}
	
	private Point2D.Double polarToCartestian(double radius,double angle)
		{
		return new Point2D.Double(
				Math.cos(angle)*radius,
				Math.sin(angle)*radius
				);
		}
	
	public static void main(final String[] args) {
		new Biostar336589().instanceMainWithExit(args);
	}
}
