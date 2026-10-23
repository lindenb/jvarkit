/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
import java.awt.Insets;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.w3c.dom.Document;
import org.w3c.dom.DocumentFragment;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.hershey.Hershey;
import com.github.lindenb.jvarkit.illumina.IlluminaConstants;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.lang.HasName;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.locatable.SimplePosition;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.XLINK;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.LocatableDelegate;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.svg.SVGDocument;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;

public class BamToDOM extends SVGDocument {
	private static final int HEIGHT_MAIN_TITLE=30;
	private static final int HEIGHT_SAMPLE_NAME=20;
	
	private enum GroupBy {
		
	}
	
	private static final Logger LOG = Logger.of(BamToDOM.class);

	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private File outputFile = null;
	@Parameter(names={"-i","--interval","--region"},description=IntervalParser.OPT_DESC,required=true)
	private String intervalStr = null;
	@Parameter(names={"-w","--width"},description="Page width")
	private int drawinAreaWidth = 1000 ;
	@Parameter(names={"-c","--showclipping","--clip"},description="Show clipping")
	private boolean showClipping = false;
	@Parameter(names={"-S","--vcf"},description="Indexed VCF. the Samples's name must be the same than in the BAM")
	private Path vcf = null;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition samRecordPartition = SAMRecordPartition.sample;
	//@Parameter(names={"--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	//private SamRecordFilter samRecordFilter = SamRecordFilterFactory.getDefault();
	@Parameter(names={"--mapq"},description="min mapping quality")
	private int mappingQuality =  1;
	@Parameter(names={"--prefix"},description="file prefix")
	private String prefix="";
	@Parameter(names={"--bases"},description="print bases in read")
	private boolean print_bases = false;
    @Parameter(names= {"--gff","--gff3"},description="Optional Tabix indexed GFF3 file.")
    private Path gtfFile = null;

	@DynamicParameter(names = "-D", description = "custom parameters. '-Dkey=value'. Undocumented.")
	private Map<String, String> dynaParams = new HashMap<>();

	private final Hershey hershey=new Hershey();
	private Path fastaReference = null;
	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private SAMSequenceDictionary referenceDict=null;
	private GenomicSequence genomicSequence=null;
	private final List<Region> all_regions = new ArrayList<>();
	
	
	private class SamRecordPair implements HasName
		{
		private SAMRecord R1 = null;
		private SAMRecord R2 = null;
		
		public SamRecordPair(final SAMRecord rec)
			{
			this.R1 = rec;
			}
		@Override
		public String getName() {
			return R1.getReadName();
			}
		}
	
	private class Region extends LocatableDelegate<Locatable> {
		double x;
		double width;
		double y=0.0;//may change
		double height=0.0; //may change
		Region(final Locatable loc) {
			super(loc);
			}
		public double getX() {
			return x;
			}
		public double getWidth() {
			return width;
			}
		public double getMaxX() {
			return getX()+getWidth();
			}
		int trimPosition(int genomic1) {
			return Math.min(Math.max(genomic1, this.getStart()),this.getEnd());
			}
		
		double getBaseWidth() {
			return (1.0/(double)getLengthOnReference())*this.getWidth();
			}
		
		double baseToPixel(int genomic1) {
			return this.getX() + ((genomic1-getStart())/(double)getLengthOnReference())*this.getWidth();
			}
		
		double baseToPixelTrim(int genomic1) {
			return Math.min(Math.max(baseToPixel(genomic1), this.getX()),this.getMaxX());
			}
		
		QueryInterval toQueryInterval(int extend) {
			return new QueryInterval(referenceDict.getSequenceIndex(getContig()),
					Math.max(getStart()-extend,1), 
					getEnd()+extend
					);
			}
		
		boolean overlaps(int start2,int end2) {
			return CoordMath.overlaps(getStart(), getEnd(), start2, end2);
			}
		GenomicSequence getGenomicSequence() {
			return BamToDOM.this.getGenomicSequence(getContig());
			}
		
		// BAM stuff
		final Map<String,List<SAMRecord>> group2buffer= new HashMap<>();
		
		void addSAMRecord(SAMFileHeader header, final SAMRecord rec) {
			final String groupName = getGroupForSAMRecord(header, rec);
			List<SAMRecord> L = group2buffer.get(groupName);
			if(L==null) {
				L= new ArrayList<SAMRecord>();
				group2buffer.put(groupName, L);
				}
			L.add(rec);
			}
		
		Node plotBam(SAMFileHeader header) {
			Element g = BamToDOM.this.group();
			
			return g;
			}
		
		
		}
	
	
	
	public String getGroupForSAMRecord(final SAMFileHeader samHeader,final SAMRecord rec) {
		return null;
		}
	
	public void plot(Path path,SamReader samReader) {
		final int extend=150;
		final SAMFileHeader samHeader = samReader.getFileHeader();
		final String sampleName= samHeader.getReadGroups().stream()
				.map(RG->RG.getSample())
				.filter(SN->!StringUtils.isBlank(SN))
				.findFirst()
				.orElse(IOUtils.getFilenameWithoutCommonSuffixes(path));
		final IntervalTreeMap<Region> regionMap=new IntervalTreeMap<BamToDOM.Region>();
		for(Region r: this.all_regions) {
			regionMap.put(new Interval(r), r);
			}
		try(CloseableIterator<SAMRecord> iter=samReader.query(optimizeIntervals(extend), false)) {
			while(iter.hasNext()) {
				final SAMRecord rec=iter.next();
				regionMap.getOverlapping(rec).stream().forEach(R->R.addSAMRecord(samHeader,rec));
				}
			}
		final Set<String> groups = all_regions.stream().flatMap(R->R.group2buffer.keySet().stream()).collect(Collectors.toSet());
		for(String groupName: groups) {
			
			}
		
		for(Region r: this.all_regions) {
			
			}

		}
	
	private QueryInterval[]  optimizeIntervals(int extend) {
		return QueryInterval.optimizeIntervals( this.all_regions
				.stream()
				.map(R->R.toQueryInterval(extend))
				.toArray(n->new QueryInterval[n])
				);
		}
	
	private Collection<SamRecordPair> viewAsPaired(final List<SAMRecord> records) {
		final Map<String,SamRecordPair> hash=new HashMap<String, BamToDOM.SamRecordPair>();
		for(SAMRecord rec:records) {
			final String name=rec.getReadName();
			SamRecordPair p = hash.get(name);
			if(p==null)  {
				hash.put(name, new SamRecordPair(rec));
				}
			else if(p.R2!=null)
				{
				p.R2=rec;
				}
			}
		return hash.values();
		}
	
	public void setShowClipping(boolean showClipping) {
		this.showClipping = showClipping;
		}
	public boolean isShowClipping() {
		return showClipping;
		}
	
	private int left(final SAMRecord rec)
		{
		return (isShowClipping()?rec.getUnclippedStart():rec.getAlignmentStart());
		}

	private int right(final SAMRecord rec)
		{
		return (isShowClipping()?rec.getUnclippedEnd():rec.getAlignmentEnd());
		}
	
	
		
	/** get read name that is displayed in the popup window */
	private String getReadNameForPopup(final SAMRecord rec) {
		final StringBuilder sb = new StringBuilder(rec.getReadName());
		
		if(rec.getReadPairedFlag()) {
			if(rec.getFirstOfPairFlag()) sb.append(" 1");
		if(rec.getSecondOfPairFlag()) sb.append(" 2");
		sb.append("/2");
		sb.append("\n");
		sb.append("pos: ").append(rec.getContig()+":"+StringUtils.niceInt(rec.getStart())+"-"+StringUtils.niceInt(rec.getEnd())).append("\n");
		if(rec.getStart()!=rec.getUnclippedStart() || rec.getEnd()!=rec.getUnclippedEnd()) {
			sb.append("unclipped: ").append(rec.getContig()+":"+StringUtils.niceInt(rec.getUnclippedStart())+"-"+StringUtils.niceInt(rec.getUnclippedEnd())).append("\n");
		}

		sb.append("strand:").append(rec.getReadNegativeStrandFlag()?" -":" +").append("\n"); 
		sb.append("flag:").append(rec.getFlags()).append("\n");
		sb.append("MAPQ:"+rec.getMappingQuality()).append("\n");
		if(rec.getMateUnmappedFlag())  {
			sb.append(" mate:unmapped\n");
			}
		else
			{
			sb.append("mate:"+rec.getMateReferenceName()+":"+rec.getMateAlignmentStart()+" "+(rec.getMateNegativeStrandFlag()?" -":" +")).append("\n");
			sb.append("tlen:"+StringUtils.niceInt(rec.getInferredInsertSize())).append("\n");
			}
		}
		if(rec.getReadFailsVendorQualityCheckFlag())  sb.append("fails-quality\n");
		if(rec.getDuplicateReadFlag()) sb.append("duplicate\n");
		if(rec.isSecondaryAlignment())  sb.append("secondary\n");
		if(rec.getSupplementaryAlignmentFlag())  sb.append("supplementary\n");
		for(SAMTagAndValue stav : rec.getAttributes()) {
			sb.append(stav.tag).append(":").append(stav.value).append("\n");
			}
		sb.append("\n");
		return sb.toString().trim();
		}
	
	
	private void defsForBases(double w,double h) {

		//base
		for(String base:new String[]{"a","A","t","T","g","G","c","C","n","N"})
			{
			final Element path = super.path();
			super.defsElement.appendChild(path);
			path.setAttribute("id","b"+base);
			path.setAttribute("class","b"+base.toUpperCase());
			path.setAttribute("d",this.hershey.svgPath(
					base,
					0,
					0,
					w,
					h
					));
			setTitle(path,base);
			}
		}
	
	private GenomicSequence getGenomicSequence(final String ctg) {
		if(this.genomicSequence==null || !this.genomicSequence.hasName(ctg)) {
			this.genomicSequence= new GenomicSequence(this.indexedFastaSequenceFile,ctg);
			}
		return this.genomicSequence;
		}
	
	private Node writeRef( ) {
		double y=0;
		final Node root= group();
		/* write REFERENCE */
		root.appendChild(comment("BEGIN REFERENCE"));
		for(Region region : this.all_regions) {
			final Node g= group();
			root.appendChild(g);
			for(int pos=region.getStart();
					region.getBaseWidth() > 5 && //ignore if too small
					pos<=region.getEnd();++pos)
				{
				final char c= region.getGenomicSequence().charAt(pos-1);
				final Element use = this.use(
						 region.baseToPixel(pos),
						 y,
						 "#b"+c
						 );
				
				setTitle(use, String.valueOf(c)+" "+StringUtils.niceInt(pos));
				g.appendChild(use);
				}
			}
		root.appendChild(comment("END REFERENCE"));
		return root;
		}
	
	Element ruler() {
		double y=0;
		double height=10;
		final double featureHeight = 10;
		Element root= super.group();
		for(Region interval: this.all_regions) {
			final Element g= super.group();
			root.appendChild(g);
			/* write ruler */
			int prev_ruler_printed=-1;
			double prev_ruler_printed_x=-1;

			for(int pos= interval.getStart(); pos<= interval.getEnd();++pos)
				{
				double x= interval.baseToPixel(pos);
				if(prev_ruler_printed_x!=-1 &&  prev_ruler_printed_x+ Math.max(20.0,featureHeight) >= x)continue;
				prev_ruler_printed_x = x;
				
				Element rect = super.rect(x, y, 0.5, (height-y));
				rect.setAttribute("class", "rulerline");
				setTitle(rect,StringUtils.niceInt(pos));
				g.appendChild(rect);
				
				if(prev_ruler_printed==-1 ||  interval.baseToPixel(prev_ruler_printed)+featureHeight < x)
					{
					x= (interval.baseToPixel(pos)+interval.baseToPixel(pos+1))/2.0;
					Element txt = super.text(0, 0,StringUtils.niceInt(pos));
					txt.setAttribute("class", "rulerlabel");
					txt.setAttribute("transform", "translate("+format(x)+","+y+") rotate(90) ");
					
					g.appendChild(txt);
					
					
					prev_ruler_printed=pos;
					}
				
				}
			}
		
		return root;
		}
	
	
	private DocumentFragment defsMutations(double featureWidth) {
		final DocumentFragment root=fragment();
		//mutations
	
		Element rect = super.rect(0, 0, featureWidth, featureWidth);
		rect.setAttribute("id", "nocall");
		rect.setAttribute("class", "nocall");
		
		root.appendChild(rect);
		
		rect = super.rect(0, 0, featureWidth, featureWidth);
		rect.setAttribute("id", "homref");
		rect.setAttribute("class", "homref");
		root.appendChild(rect);
		
		rect = super.rect(0, 0, featureWidth, featureWidth);
		rect.setAttribute("id", "homvar");
		rect.setAttribute("class", "homvar");
		root.appendChild(rect);
		
		
		Element group = super.group();
		group.setAttribute("id","het");
		group.setAttribute("class","het");
		rect = super.rect(0, 0,featureWidth/2.0,featureWidth/2.0);
		group.appendChild(rect);
		rect = super.rect(featureWidth/2.0,featureWidth/2.0,featureWidth/2.0,featureWidth/2.0);
		group.appendChild(rect);
		rect = super.rect(0, 0, featureWidth, featureWidth);
		rect.setAttribute("style", "fill:none");
		group.appendChild(rect);
		root.appendChild(group);

					
				
				
		

		//mutation
		
		Element circle = super.circle((featureWidth/2.0), (featureWidth/2.0), (featureWidth/2.0));
		circle.setAttribute("id","mut");
		circle.setAttribute("class","mutation");
		setTitle(circle,"mutation");
		root.appendChild(circle);
				
		return root;
		}
	
	Element writeConsensus(SamReader samReader) {
		double consensus_y=0;
		double featureHeight =10;
		Element root= group();
		root.appendChild(comment("BEGIN_CONSENSUS"));
		final AutoMap<SimplePosition,Counter<Character>,Counter<Character>> pos2count = AutoMap.make(()->new Counter<Character>());
		
		try(CloseableIterator<SAMRecord> iter=samReader.query(optimizeIntervals(0),false)) {
			while(iter.hasNext()) {
				SAMRecord rec= iter.next();
				byte[] bases = rec.getReadBases();
				if(bases.equals(SAMRecord.NULL_SEQUENCE)) continue;
				for(AlignmentBlock ab: rec.getAlignmentBlocks())
					{
					for(int i=0;i< ab.getLength();i++) {
						final char base = (char) Character.toUpperCase(bases[(ab.getReadStart()-1)+i]);
						pos2count.insert(new SimplePosition(rec.getContig(), ab.getReferenceStart()+i)).incr(base);
						}
					}
				}
			}
				
		/* write consensus */
		final int max_depth = (int) pos2count.entrySet().stream().mapToLong(KV->KV.getValue().getTotal()).max().orElse(0L);
			for(Region interval:this.all_regions) {
			for(int pos= interval.getStart();
				max_depth > 0L &&
				interval.getBaseWidth() >5 && //ignore if too small
				pos<= interval.getEnd();++pos)
				{
				final char rc= interval.getGenomicSequence().charAt(pos-1);
				final Counter<Character> cons = pos2count.get(new SimplePosition(interval.getContig(),pos));
				if(cons==null) continue;
				final char ac=cons.getMostFrequent();
				final double cov_height = (cons.getTotal()/(double)max_depth)*(featureHeight-1);
				if(cons.getCountCategories()==1 && ac==rc) {
					final double x0  = interval.baseToPixel(pos);
					final Element rect = super.rect(
							x0,
							(consensus_y+ featureHeight-cov_height),
							(interval.getBaseWidth()-1),
							(cov_height)
							);
					rect.setAttribute("style", "stroke:none;fill:gray;");
					setTitle(rect,String.valueOf(ac)+" "+cons.count(ac));
					}
				else
					{
					double x0  = interval.baseToPixel(pos);
					for(Character a : cons.keySetDecreasing()) {
						final double aw = (interval.getBaseWidth()-1)*((double)cons.count(a)/cons.getTotal());
						final Element rect = super.rect(
							x0,
							(consensus_y+ featureHeight-cov_height),
							aw,
							cov_height
							);
						rect.setAttribute("style",  "stroke:none;fill:"+AcidNucleics.cssColor(a));
						setTitle(rect, ""+a+" "+cons.count(a)+"/"+cons.getTotal()+" "+(int)(100.0*(double)cons.count(a)/cons.getTotal())+"%");
						x0+=aw;
						}
					}
				}
			}
		return root;
		}
	
	

	
	
	private void printDocument(
		final XMLStreamWriter w,
		final SAMSequenceDictionary dict
		) 
		throws XMLStreamException
		{
		final UrlSupplier urlSupplier=new UrlSupplier(dict);
		final String IDT_XNS="https://umr1087.univ-nantes.fr/";
		final Insets insets=new Insets(20, 20, 20, 20);
		final Dimension dim=new Dimension(insets.left+insets.right+this.drawinAreaWidth,0);
		dim.height += (int)context.getHeight();
		//dim.height+=100;
		
		
		

		
		
		appendStyle(
				"g.maing {stroke:black;stroke-width:0.5px;fill:none;}\n"+
				".frame {stroke:darkgray;fill:none;}\n"+
				".maintitle {stroke:blue;fill:none;}\n"+
				".nocall {stroke:orange;fill:none;stroke-dasharray:2;}\n"+
				".homref {stroke:blue;fill:none;}\n"+
				".homvar {stroke:blue;fill:blue;}\n"+
				".het {stroke:blue;fill:blue;}\n"+
				"line.insert {stroke:red;stroke-width:4px;}\n"+
				".rulerline {stroke:lightgray;stroke-width:1px;}\n"+
				".rulerlabel {fill:gray;font-size:7px;}\n"+
				"text.maintitle {fill:darkgray;stroke:black;text-anchor:middle;font-size:25px;}\n"+
				"path.samplename {stroke:black;stroke-width:3px;}\n"+
				"circle.mutation {stroke:red;fill:none;stroke-width:"+(context.featureWidth<5?0.5:3.0)+"px;}"+
				"line.gene {stroke:blue;stroke-width:4px;}"+
				"rect.exon {stroke:blue;fill:orange;stroke-width:4px;}"+
				""
				);
		for(char base : new char[] {'A','C','G','T','N'}) {
			appendStyle(".b"+base+" {stroke:"+AcidNucleics.cssColor(base)+";stroke-width:"+(context.featureWidth<5?0.5:1.0)+"px}\n");
			}
		w.writeEndElement();//style
		
		w.writeStartElement("defs");
		
		
		Element g = super.group();
		
		Element text = super.text(
				(this.drawinAreaWidth/2.0), 
				y, 
				context.sampleName+" "+ SequenceDictionaryUtils.getBuildName(dict).orElse("") +" "+this.interval.toNiceString()
				);
		text.setAttribute("class","maintitle");
		//TODO append text
		
		if(hyperlink.apply(this.interval).isPresent()) {
			w.writeEndElement();//a
			}
		y+=HEIGHT_MAIN_TITLE;
		
		
		
		
		
		
		}
	
	Element variants(VCFReader vcfReader)	{
		double y=0;
		final Element root= group();
		for(Region interval: this.all_regions) {
			final Element g= group();
			root.appendChild(g);
			try(CloseableIterator<VariantContext> iter=vcfReader.query(interval)) {
				while(iter.hasNext()) {
					VariantContext ctx=iter.next();
					final Genotype gt=ctx.getGenotype(context.sampleName);
					if(gt==null) continue;
					String url = null;
					if(ctx.hasID()) {
						url = urlSupplier.of(ctx.getID()).stream().map(G->G.getUrl()).findFirst().orElse(null);
						}
					if(url==null) {
						url = hyperlink.apply(ctx).orElse(null);
						}
				
				final String svType;
				switch(gt.getType()) {
					case HET:
					case HOM_REF:
					case HOM_VAR:
					case NO_CALL: svType = gt.getType().name().toLowerCase().replace("_", "");break;
					default: svType="";
					}
				final double x0  = interval.baseToPixel(ctx.getStart());
				Element use = super.use(
						x0,
						y,
						svType
						);
				use.setAttribute("x",format(x0));
				use.setAttribute("y",format(y));
				
				
				use.setAttribute("href", "#"+svType);
				final StringBuilder sb = new StringBuilder();
				sb.append(ctx.getContig()+":"+ctx.getStart()+" ");
				sb.append(ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/")));
				sb.append("\n");
				sb.append(gt.getType().name()).append("\n");
				if(gt.isFiltered()) sb.append("FILTERED:"+gt.getFilters()+"\n");
				if(gt.hasDP()) sb.append("DP:"+gt.getDP()).append("\n");
				if(gt.hasGQ()) sb.append("GQ:"+gt.getGQ()).append("\n");
				if(gt.hasAD()) sb.append("AD:"+IntStream.of(gt.getAD()).mapToObj(V->String.valueOf(V)).collect(Collectors.joining(","))).append("\n");
				if(gt.hasPL()) sb.append("PL:"+IntStream.of(gt.getPL()).mapToObj(V->String.valueOf(V)).collect(Collectors.joining(","))).append("\n");
				setTitle(use,sb.toString());
				
				use = super.anchor(use, url);
				
				g.appendChild(use);
				}
				}
			}
	
		return root;
		}
	
	
	Node makeGTF(
			double y,
			double featureHeight,
			TabixReader tbr
				) throws IOException
			{
			final Node  root = group();
			root.appendChild(comment("BEGIN GTF"));
			for(Region region: this.all_regions) {
				final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
				 final String ctg = cvt.apply(region.getContig());
				 if(StringUtils.isBlank(ctg)) continue;
				 final Node  g = group();
				 root.appendChild(g);
				 final List<GTFLine> genes  = new ArrayList<>();
				 final List<GTFLine> exons  = new ArrayList<>();
				 final Set<String> gene_names = new HashSet<>();
				 final GTFCodec codec = new GTFCodec();
				 final TabixReader.Iterator iter=tbr.query(ctg,region.getStart(),region.getEnd());
				 for(;;) {
					 final String line = iter.next();
					 if(line==null) break;
					 if(StringUtils.isBlank(line) ||  line.startsWith("#")) continue;
					 final GTFLine gtfline = codec.decode(line);
					 if(gtfline==null) continue;
					 if(gtfline.isGene()) {
						 genes.add(gtfline);
						 final String geneName = gtfline.getGeneName();
						 if(!(StringUtils.isBlank(geneName)) || geneName.equals(".")) {
							 gene_names.add(geneName);
						 	}
					 	 }
					 else if(gtfline.getType().equals("exon")) {
						 exons.add(gtfline);
					 	 }
				 	}/* end for */
			 for(GTFLine gtf:genes) {
				final Element e= line(
						region.baseToPixelTrim(gtf.getStart()),
						(y+ featureHeight*0.5),
						region.baseToPixelTrim(gtf.getEnd()+1),
						(y+ featureHeight*0.5),
						Maps.of("class", "gene")
						);
				setTitle(e, gtf.getGeneName());
				g.appendChild(e);
			 	}
			 for(GTFLine gtf:exons) {
				 final double x1 = region.baseToPixelTrim(gtf.getStart());
				 final double x2 = region.baseToPixelTrim(gtf.getEnd()+1);
	
				 final Element e= rect(
							x1,
							y,
							(x2-x1),
							(featureHeight-1),
							Maps.of("class", "exon")
							);
					
				 setTitle(e, gtf.getGeneName());
				 g.appendChild(e);
			 	 }
			 }/* end if contig!null */
			root.appendChild(comment("END GTF"));
			return root;
	 		}
	
	private static class PseudoByteArray {
		final byte[] array;
		final int len;
		final char c;
		PseudoByteArray(byte[] array) {
			this.array=array;
			this.len = array.length;
			this.c='\0';
			}
		PseudoByteArray(int len, char c) {
			this.array=null;
			this.len=len;
			this.c=c;
			}
		public int length() {
			return this.len;
			}
		public char at(int i) {
			return this.array==null?this.c:(char)this.array[i];
			}
		}
	
	Node makeSamRecord(
			final Region region,
			final SAMRecord record,
			final double featureHeight
			)	{
			double mid_y=  featureHeight/2.0;
			final double y_h95= featureHeight*0.90;
			final double y_top5= (featureHeight-y_h95)/2.0;
			final double y_bot5= y_top5+y_h95;
			final double arrow_w= featureHeight/3.0;
	
			
			/* print that sam record */
			final int unclipped_start= record.getUnclippedStart();
			final Cigar cigar = record.getCigar();
			if(cigar==null) return fragment();
			final PseudoByteArray bases = 
					record.getReadBases().equals(SAMRecord.NULL_SEQUENCE)?
					new PseudoByteArray(cigar.getReadLength(),'N'):
					new PseudoByteArray(record.getReadBases());
					;
			
			final PseudoByteArray qualities =
					record.getBaseQualities().equals(SAMRecord.NULL_QUALS)?
					new PseudoByteArray(cigar.getReadLength(),'#'):
					new PseudoByteArray(record.getBaseQualities());

			
			Element g = super.group();
			setTitle(g,getReadNameForPopup(record));
			

			
			int readPos=0;
			final Map<Integer,String> pos2insertions=new HashMap<Integer,String>();
			final List<CigarElement> cigarElements= cigar.getCigarElements();
			
			/* find position of arrow */
			int arrow_cigar_index=-1;
			for(int cidx=0; cidx< cigarElements.size(); cidx++ )
				{
				final CigarElement ce = cigarElements.get(cidx);
				final CigarOperator op=ce.getOperator();
				
				switch(op)
					{
					case H:case S: if(!this.showClipping) break;//threw
					case M:case EQ: case X:
						{
						arrow_cigar_index=cidx;
						}
					default:break;
					}
				if(record.getReadNegativeStrandFlag() && arrow_cigar_index!=-1)
					{
					break;
					}
				}
	
			
			int refPos=unclipped_start;
			
			/* loop over cigar string */
			for(int cidx=0; cidx< cigarElements.size(); cidx++ )
				{
				final CigarElement ce = cigarElements.get(cidx);
				final CigarOperator op=ce.getOperator();
				boolean in_clip=false;
				switch(ce.getOperator())
					{
					case D://cont
					case N:
						{
						final double x1 = region.baseToPixelTrim(refPos);
						final double x2 = region.baseToPixelTrim(refPos+ ce.getLength());
						if(x1<x2 && region.overlaps(refPos, refPos+ ce.getLength()) )
							{
							final Element line = super.line(
									x1,
									mid_y,
									x2,
									mid_y,
									Maps.of("class","indel")
									);
							setTitle(line, op.name()+ce.getLength());
							g.appendChild(line);
							}						
						refPos += ce.getLength();
						break;
						}
					case I: 
						{
						if(region.overlaps(refPos, refPos) ) {
							final StringBuilder sb=new StringBuilder(ce.getLength());
							for(int i=0;i< ce.getLength();++i)
								{
								sb.append((char)bases.at(readPos+i));
								}
							pos2insertions.put(refPos, sb.toString());
							}
						readPos+=ce.getLength();
						break;
						}
					case H:
						{
						if(!this.showClipping)
							{
							refPos+=ce.getLength();
							break;
							}
						in_clip=true;
						//NO break;
						}
					case S:
						{
						if(!this.showClipping)
							{
							readPos+=ce.getLength();
							refPos+=ce.getLength();
							break;
							}
						in_clip=true;
						//NO break;
						}
					case X:
					case EQ:
					case M:
						{
						int match_start = refPos;
						int match_end = refPos + ce.getLength();
						final Element path;
						//print sam background
						final StringBuilder sb=new StringBuilder();
						if(record.getReadNegativeStrandFlag() &&
							match_start >= region.getStart() && 
							match_start <= region.getEnd() &&
							cidx==arrow_cigar_index)
							{
							path = this.createShape()
								.moveTo( region.baseToPixel(match_start)+arrow_w, y_top5)
								.horizonalR( region.baseToPixel(trim(match_end))-region.baseToPixel(match_start)-arrow_w)
								.verticalR(y_h95)
								.make()
								;
							
							sb.append(" M ").append(format(baseToPixel(match_start)+arrow_w)).append(',').append(y_top5);
							sb.append(" h ").append(format());
							sb.append(" v ").append(format(y_h95));
							sb.append(" h ").append(format(-(baseToPixel(trim(match_end))-baseToPixel(match_start)-arrow_w)));
							sb.append(" L ").append(format(baseToPixel(match_start))).append(',').append(mid_y);
							sb.append(" Z");
							}
						else if(!record.getReadNegativeStrandFlag() &&
								match_end >= region.getStart() && 
								match_end <= region.getEnd() &&
								cidx==arrow_cigar_index
								)
							{
							sb.append(" M ").append(format(baseToPixel(match_end)-arrow_w)).append(',').append(y_top5);
							sb.append(" h ").append(format(-(baseToPixel(match_end)-baseToPixel(trim(match_start))-arrow_w)));
							sb.append(" v ").append(format(y_h95));
							sb.append(" h ").append(format(baseToPixel(match_end)-baseToPixel(trim(match_start))-arrow_w));
							sb.append(" L ").append(format(baseToPixel(match_end))).append(',').append(mid_y);
							sb.append(" Z");
							}
						else
							{
							sb.append(" M ").append(format(baseToPixel(trim(match_start)))).append(',').append(y_top5);
							sb.append(" h ").append(format(baseToPixel(trim(match_end))-baseToPixel(trim(match_start))));
							sb.append(" v ").append(format(y_h95));
							sb.append(" h ").append(format(-(baseToPixel(trim(match_end))-baseToPixel(trim(match_start)))));
							sb.append(" Z");
							}
						
						path.setAttribute("d",sb.toString().trim());
						if( ce.getOperator()==CigarOperator.H || ce.getOperator()==CigarOperator.S)
							{
							path.setAttribute("style","fill:yellow;");
							}
						else
							{
							final int[] rgb = new int[3];
							if(record.getReadFailsVendorQualityCheckFlag()) {
								rgb[0]=255; rgb[1]=102; rgb[2]=255;
								}
							else if(record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !record.getMateReferenceName().equals(this.interval.getContig())) {
								rgb[0]=102; rgb[1]=178; rgb[2]=255;
								}
							else if(record.getReadPairedFlag() && record.getMateUnmappedFlag()) {
								rgb[0]=51; rgb[1]=255; rgb[2]=51;
								}
							else if(record.getReadPairedFlag() && !record.getProperPairFlag()) {
								rgb[0]=153; rgb[1]=76; rgb[2]=0;
								}
							else
								{
								rgb[0]=192; rgb[1]=192; rgb[2]=192;
								}
							final double f= Math.max(record.getMappingQuality(),60)/60.0;
							for(int i=0;i< rgb.length;i++) rgb[i] = (int)(rgb[i]*f);
							path.setAttribute("style","fill:rgb("+rgb[0]+","+rgb[1]+","+rgb[2]+")");
							}
						
						if(op.consumesReadBases())
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final int pos_i = refPos+i;
								boolean pos_in_interval =  (region.getStart() <= pos_i && pos_i<=region.getEnd() );
								char ca=(char)bases.at(readPos);
								if(!in_clip)
									{
									Counter<Character> count= consensus.get(pos_i) ;
									if(count==null)
										{
										count=new Counter<>();
										consensus.put(refPos+i,count) ;
										}
									count.incr(ca);
									}
								char cb= genomicSequence.charAt(refPos+i-1);
								
								/* mutation ? */
								if(this.print_bases && cb!='N' && ca!='N' && cb!=ca && !in_clip && pos_in_interval)
									{
									final Element use = super.element("use");
									use.setAttribute("x",format( baseToPixel(refPos+i)));
									use.setAttribute("href", "#mut");
									
									}
								
								if( pos_in_interval && 
									((this.print_bases&& context.featureWidth>5) || (cb!=ca ))
									 )
									{
									final Element use = super.element("use");
									use.setAttribute("x",format( baseToPixel(refPos+i)));
									//w.setAttribute("y",format(y_top));never needed
									use.setAttribute("href", "#b"+ca);
									}
								readPos++;
								}
							}
						
	
						
						refPos+=ce.getLength();
						break;
						}
					default:
						{
						throw new RuntimeException("Unknown SAM op: "+ce.getOperator());
						}
					}
				}
			
			for(Integer pos:pos2insertions.keySet())
				{
				if(pos < region.getStart())  continue;
				if(pos > region.getEnd())  continue;
				final String insertion=pos2insertions.get(pos);
				
				
				final double x= region.baseToPixel(pos);
				final Element line = super.line(x, x, y_top5, y_bot5);
				setTitle(line, "Insertion "+insertion);
				line.setAttribute("class", "insert");
				g.appendChild(line);
				}

			
			return g;
			}
		
	
			
			
					
		
		
	public void setDrawinAreaWidth(int drawinAreaWidth) {
		this.drawinAreaWidth = drawinAreaWidth;
		}
		
	}
