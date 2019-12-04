package com.github.lindenb.jvarkit.tools.sashimi;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.DoubleFunction;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Result;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Exon;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Intron;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

@Program(name="plotsashimi",
description="Print Sashimi plots from Bam",
keywords={"bam","visualization","svg","rna","exon"},
modificationDate="20191104",
creationDate="20191104"
)
public class PlotSashimi extends Launcher {
private static final Logger LOG = Logger.build(PlotSashimi.class).make();

@Parameter(names={"-g","--gtf"},description=GtfReader.OPT_DESC)
private Path gtfPath = null;

@Parameter(names={"-r","--region","--interval"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,required=true)
private IntervalListProvider intervalListProvider= IntervalListProvider.empty();

@Parameter(names={"-R","--reference"},description="For Reading CRAM. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
private Path faidx = null;

@Parameter(names={"-o","--out"},description=ArchiveFactory.OPT_DESC,required=true)
private Path outputFile=null;

@SuppressWarnings("serial")
@DynamicParameter(names = "-D", description = "Other parameters.")
public Map<String, String> dynamicParams = new HashMap<String,String>() {{{
	put("width","1000");
	put("coverage.height","300");
	}}};

private final IntervalTreeMap<Transcript> transcriptMap = new IntervalTreeMap<>();
private Document document;
private final DecimalFormat decimalFormater = new DecimalFormat("##.##");

private Element element(final String tag) {
	return this.document.createElementNS(SVG.NS, tag);
	}
private Text text(final Object o) {
	return this.document.createTextNode(o==null?"":String.valueOf(o));
	}
private Element element(final String tag,final Object content) {
	final Element E = element(tag);
	E.appendChild(text(content));
	return E;
	}
/** convert double to string */
private String format(double v)
	{
	return this.decimalFormater.format(v);
	}

private void plotSashimi(
	final ArchiveFactory archive,
	final SamReader samReader,
	final Locatable interval
	) {
	final int drawing_width = Math.max(100,Integer.parseInt(this.dynamicParams.getOrDefault("width","1000")));
	final int coverageHeight =  Math.max(100,Integer.parseInt(this.dynamicParams.getOrDefault("coverage.height","300")));
	final double pixelperbase =  drawing_width/(double)interval.getLengthOnReference();
	
	final Function<Integer, Double> pos2pixel = POS-> (POS - interval.getStart())/(double)interval.getLengthOnReference() * drawing_width;
	
	final Counter<SimpleInterval> gaps = new  Counter<>();
	final int coverage[] = new int[interval.getLengthOnReference()];
	try(SAMRecordIterator iter=samReader.queryOverlapping(interval.getContig(), interval.getStart(), interval.getEnd())) {
		while(iter.hasNext()) {
			final SAMRecord rec = iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(rec.isSecondaryOrSupplementary()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			final Cigar cigar =rec.getCigar();
			if(cigar==null || cigar.isEmpty()) continue;
			int ref= rec.getAlignmentStart();
			for(final CigarElement ce:cigar) {
				if(ref> interval.getEnd()) break;
				final CigarOperator op = ce.getOperator();
				
				if(op.equals(CigarOperator.N))  {
					gaps.incr(new SimpleInterval(rec.getContig(),ref,ref+ce.getLength()-1));
					}

				
				if(op.consumesReferenceBases()) {
					if(op.consumesReadBases()) {
						for(int x=0;x<ce.getLength();++x) {
							final int pos1 = ref+x;
							if(pos1< interval.getStart()) continue;
							if(pos1> interval.getEnd()) break;
							coverage[pos1-interval.getStart()]++;
							}
						}
					ref+=ce.getLength();
					}
				}
			}
		}
	
		final int max_coverage = Math.max(1,Arrays.stream(coverage).max().orElse(0));


	
		while(this.document.hasChildNodes()) {
			this.document.removeChild(this.document.getFirstChild());
			}
	
		final Element svgRoot = element("svg");
		this.document.appendChild(svgRoot);
		
		final Element title = element("title");
		svgRoot.appendChild(title);
		title.appendChild(text(interval.toString()));
		
		final Element style = element("style");
		svgRoot.appendChild(style);
		style.appendChild(text(
				".coverage { background-color:red;}" +
				".frame { fill:none; color: darkgray;}"
				));

		
		final Element descr = element("desc");
		svgRoot.appendChild(descr);
		descr.appendChild(text("Author: Pierre Lindenbaum"));
	
		final Element maing = element("g");
		svgRoot.appendChild(maing);

		
		
		final Element covPath = element("path");
		covPath.setAttribute("class", "coverage");
		maing.appendChild(covPath);
		
		int y=0;
		final StringBuilder sb = new StringBuilder();
		sb.append( "M 0 "+format(y+coverageHeight));
		for(int k=0;k< coverage.length;k++)
			{
			if(k+1< coverage.length && coverage[k]==coverage[k+1]) continue; 
			final double dpy= y+ coverageHeight - coverageHeight*(coverage[k]/(double)max_coverage);
			sb.append(" L "+format(pixelperbase*k)+" "+format(dpy));
			}
		sb.append(" L "+format(drawing_width)+" "+format(y+ coverageHeight));
		sb.append(" Z");
		covPath.setAttribute("d", sb.toString());
		covPath.appendChild(element("title","Coverage. Max:"+StringUtils.niceInt(max_coverage)));
		y+=coverageHeight;
		y+=2;
		
		 boolean drawAbove  = true;
		/* plot arc */
			if(!gaps.isEmpty()) {
				int max_occurence = (int)gaps.count(gaps.getMostFrequent());
				for(final SimpleInterval intron:gaps.keySet()) {
					final int junctionStart = intron.getStart()-1;
	                final int junctionEnd = intron.getEnd()+1;
					if(!CoordMath.encloses(interval.getStart(), interval.getEnd(), junctionStart, junctionEnd)) continue;
	                
					double xstart = pos2pixel.apply(junctionStart);
					double xend = pos2pixel.apply(junctionEnd);
					double ystart=  y+ coverageHeight - coverageHeight*(coverage[junctionStart - interval.getStart()]/(double)max_coverage);
					double yend=  y+ coverageHeight - coverageHeight*(coverage[junctionStart- interval.getStart()]/(double)max_coverage);

					final Element arc= element("path");
					sb.setLength(0);
					
					double x_mid = (xend - xstart)/2.0;
					double x2 = xstart + x_mid;
					double y2= ystart+(yend-ystart)/2.0 + (drawAbove?-1:1)*x_mid;
					
					sb.append("M "+format(xstart)+" "+format(ystart));
					sb.append("Q "+format(x2)+" "+format(y2)+" "+format(xend)+" "+format(yend));
					
					arc.setAttribute("p", sb.toString());
					arc.setAttribute("class","arc");
					arc.setAttribute("style", "stroke-width:5px;");
					arc.appendChild(element("title",new SimpleInterval(interval.getContig(),junctionStart,junctionEnd).toNiceString()));
					covPath.appendChild(arc);
					
	                drawAbove = !drawAbove;
				}
			}

		final Element transcripts_g = element("g");
		maing.appendChild(transcripts_g);
		int transcript_height = Math.max(10, Integer.parseInt(this.dynamicParams.getOrDefault("transcript.height", "12")));
		for(final Transcript transcript: this.transcriptMap.getOverlapping(interval)) {
			final Element transcript_g = element("g");
			transcripts_g.appendChild(transcript_g);
			final Element tr = element("line");
			tr.setAttribute("x1", format(0));
			tr.setAttribute("y1", format(y+transcript_height/2.0));
			tr.setAttribute("x2", format(0));
			tr.setAttribute("y2", format(y+transcript_height/2.0));
			transcript_g.appendChild(tr);
			
			for(final Exon exon:transcript.getExons()) {
				if(!exon.overlaps(interval)) continue;
				final Element exon_rect = element("rect");
				exon_rect.setAttribute("x", format(0));
				exon_rect.setAttribute("y", format(y+transcript_height/2.0));
				exon_rect.appendChild(element("title",exon.getName()));
				transcript_g.appendChild(exon_rect);
				}
			y+=transcript_height+1;
		}
		
		final Element frame_rect = element("frame");
		frame_rect.setAttribute("class", "frame");
		frame_rect.setAttribute("x", "0");
		frame_rect.setAttribute("y", "0");
		frame_rect.setAttribute("width", format(drawing_width));
		frame_rect.setAttribute("height",format(y));
		
		svgRoot.setAttribute("width",format(drawing_width+1));
		svgRoot.setAttribute("height",format(y+1));
		
		
		try {
		
		final Transformer tr = TransformerFactory.newInstance().newTransformer();
		final Result result;
		
		String fname="jeter.svg";
		try(final PrintWriter pw=archive.openWriter(fname)) {
			result = new StreamResult(pw);
			tr.transform(new DOMSource(this.document),result);
			pw.flush();
			}
		
		} catch(Exception err) {
			throw new RuntimeException(err);
		}

	}

@Override
public int doWork(final List<String> args) {
	ArchiveFactory archive=null;
	try
		{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		dbf.setNamespaceAware(true);
		DocumentBuilder db = dbf.newDocumentBuilder();
		this.document = db.newDocument();

		
		final SamReaderFactory srf = super.createSamReaderFactory();
		if(faidx!=null) {
			srf.referenceSequence(this.faidx);
			}
		if(this.gtfPath!=null) {
			try(GtfReader gtfReader=new GtfReader(this.gtfPath)) {
				if(this.faidx!=null) gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(this.faidx)));
				gtfReader.
					getAllGenes().
					stream().
					flatMap(G->G.getTranscripts().stream()).
					forEach(T->this.transcriptMap.put(new Interval(T), T));
				}
			}
		
		archive = ArchiveFactory.open(this.outputFile);
		for(final Path bam: IOUtils.unrollPaths(args)) {
			try(SamReader sr = srf.open(bam)) {
				if(!sr.hasIndex()) {
					LOG.error("Bam is not indexed "+bam);
					return -1;
				}
				final SAMFileHeader header= sr.getFileHeader();
				final ArchiveFactory final_archive = archive;
				this.intervalListProvider.
					dictionary(header.getSequenceDictionary()).
					stream().
					forEach(R->{
						plotSashimi(final_archive,sr,R);
					});
					
				}
			}
		return 0;
		} 
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	finally {
		CloserUtil.close(archive);
		}
	}

public static void main(final String[] args) {
	new PlotSashimi().instanceMainWithExit(args);
	}
}
