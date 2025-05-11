package com.github.lindenb.jvarkit.tools.manhattan;

import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC

## Node

this program is very unstable. I often change everything...


END_DOC

*/
@Program(
		name="manhattan",
		description="Manhattan plot SVG picture from different sources.",
		keywords={"chromosome","reference","chart","visualization","svg"},
		creationDate="20220525",
		modificationDate="20240525",
		generate_doc=true,
		jvarkit_amalgamion = true
		)
public class Manhattan extends Launcher {
	
	private static final Logger LOG = Logger.of(Manhattan.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxFile = null;
	@Parameter(names={"--bed","-L"},description="BED file containing one or more region to observe. Use the whole chromosome if undefined")
	private Path bedPath = null;
	@DynamicParameter(names={"-D","--define"},description="Dynamic parameters")
	private Map<String,String> __dynaParams = new HashMap<>();
	@Parameter(names={"--dimension"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,300);
	@Parameter(names="--include",description="include chromosomes matching the following expression")
	private String includeContigsRegex = "(chr)?[0-9XY]+";
	@Parameter(names="-l",description="List available handlers and exit")
	private boolean list_handlers = false;
	@Parameter(names="-n",description="handler name")
	private String handlerName="default";


	private AttributeMap dynaParams = AttributeMap.wrap(Collections.emptyMap());
	
	private class Item extends SimpleInterval {
		final double value;
		final String label;
		Item(String ctg,int start,int end,double value,final String label) {
			super(ctg,start,end);
			this.value = value;
			this.label = label;
			}
		}
	
	private abstract class Handler {
		protected List<Track> tracks = new ArrayList<>();
		protected Map<String,Track> name2track = new HashMap<>();
		protected XMLStreamWriter w;
		protected List<Region> regions = new ArrayList<>();
		private Dimension dimension=null;
		protected String filename="undefined";
		protected Function<Locatable, String> loc2url = L->null;
		protected Function<Locatable, String> loc2title = L->null;
		protected Function<Locatable, String> loc2clazz = L->null;
		protected Function<Locatable, String> loc2style = L->null;
		protected Function<Locatable, String> loc2shape = L->null;
		
		class RectPainter {
			String _style = null;
			String _title = null;
			String _href = null;
			String _clazz = null;
			String _shape = null;
			private final Locatable _container;
			private final double _x0;
			private final double _width0;
			RectPainter style(String s) { this._style = s; return this;}
			RectPainter clazz(String s) { this._clazz = s; return this;}
			RectPainter title(String s) { this._title = s; return this;}
			RectPainter href(String s) { this._href = s; return this;}
			RectPainter shape(String s) { this._shape = s; return this;}
			RectPainter reset() {
				return style(null).clazz(null).title(null).href(null).shape(null);
				}
			RectPainter(final Locatable loc,double x,double width) {
				this._container = loc;
				this._x0 = x;
				this._width0 = width;
				}
			RectPainter(final Region rgn) {
				this(rgn,rgn.x,rgn.width);
				}
			RectPainter(final Pane pane) {
				this(pane.region);
				}
			double toPixel(int loc) {
				loc = Math.max(this._container.getStart(), loc);
				loc = Math.min(this._container.getEnd()+1, loc);
				return this._x0 + ((loc - this._container.getStart())/(double)this._container.getLengthOnReference())*this._width0;
				}
			
			void paintY1Y2(final Locatable loc,double y1,double y2) throws XMLStreamException {
				if(!_container.overlaps(loc)) return;
				double x0 = toPixel(loc.getStart());
				double x1 = toPixel(loc.getEnd()+1);
				final double minsize = dynaParams.getDoubleAttribute("min-width").orElse(0.1);
				if((x1-x0) < minsize) {
					double mid = (x1+x0)/2.0;
					x0 = mid - minsize/2.0;
					x1 = mid + minsize/2.0;
					}
				
				if(!StringUtils.isBlank(_href)) {
					w.writeStartElement("a");
					w.writeAttribute("href", _href);
					}
				if(this._shape!=null && this._shape.equalsIgnoreCase("circle")) {
					w.writeStartElement("circle");
					w.writeAttribute("cx",format((x0+x1)/2.0));
					w.writeAttribute("cy",format((y1+y2)/2.0));
					w.writeAttribute("r",format((x1-x0)/2.0));
					}
				else
					{
					w.writeStartElement("rect");
					w.writeAttribute("x",format(x0));
					w.writeAttribute("y",format(y1));
					w.writeAttribute("width",format(x1-x0));
					w.writeAttribute("height",format(y2-y1));
					}
				
				if(!StringUtils.isBlank(_style)) {
					w.writeAttribute("style", _style);
					}
				if(!StringUtils.isBlank(_clazz)) {
					w.writeAttribute("class", _clazz);
					}
				

				
				if(!StringUtils.isBlank(_title)) {
					w.writeStartElement("title");
					w.writeCharacters(_title);
					w.writeEndElement();
					}
				w.writeEndElement();
				if(!StringUtils.isBlank(_href)) {
					w.writeEndElement();
					}
				}
			}
		
		
		class Region  extends SimpleInterval {
			private String label = null;
			double x;
			double width;
			Region(final Locatable ssr) {
				super(ssr);
				}
			String getName() {
				return StringUtils.isBlank(label)?super.toNiceString():this.label;
				}
			
			}
	
		/** an horizontal track */
		class Track {
			protected List<Locatable> items = new ArrayList<>();
			double y =0;
			double height=0.0;
			final String label;
			Track(final String label) {
				this.label=label;
				}
			public String getName() {
				return label;
				}
		
			}
	
	
		/** intersection track & display range */
		class Pane
			{
			final List<Locatable> items = new ArrayList<>();
			final List<List<Locatable>> rows = new ArrayList<>();
			final Track track;
			final Region region;
			Pane(final Region region,final Track track) {
				this.region = region;
				this.track = track;
				}
			
			
			Pane pileup() {
				Collections.sort(this.items,comparator);
				RectPainter painter = new RectPainter(this);
				final double d = dynaParams.getDoubleAttribute("item-distance").orElse(1.0);
				while(!this.items.isEmpty()) {
					final Locatable curr  = this.items.remove(0);
					int y;
					for(y=0;y<this.rows.size();y++) {
						final List<Locatable> row = this.rows.get(y);
						Locatable last  = row.get(row.size()-1);
						if(painter.toPixel(last.getEnd()+1) + d <= painter.toPixel(curr.getStart())) {
							row.add(curr);
							break;
							}
						}
					if(y==this.rows.size()) {
						final List<Locatable> row = new ArrayList<>();
						row.add(curr);
						this.rows.add(row);
						}
					}
				return this;
				}
			}

		protected boolean overlaps(final Locatable loc) {
			return this.regions.stream().anyMatch(R->R.overlaps(loc));
			}
		
		protected Track createTrack(final String s) {
			return new Track(s);
			}
		
		
		protected Track getTrackByName(final String s) {
			Track t = this.name2track.get(s);
			if(t==null) {
				t = createTrack(s);
				this.tracks.add(t);
				this.name2track.put(s, t);
				}
			return t;
			}
		
		protected Region  newRegion(final Locatable loc) {
			final Region rgn = new Region(loc);
			this.regions.add(rgn);
			return rgn;
			}
		
		
		protected List<Pane> createPanes() {
			final List<Pane> panes = new ArrayList<>();
			for(Track track : this.tracks) {
				for(Region region: this.regions) {
					final List<Locatable> L = track.items.stream().
							filter(R->R.overlaps(region)).
							collect(Collectors.toList());
					if(L.isEmpty()) continue;
					final Pane pane = new Pane(region,track);
					pane.items.addAll(L);
					panes.add(pane);
					}
				track.items.clear();
				}
			return panes;
			}
		
		public String getName() { return this.getClass().getSimpleName();}
		public String getDescription() { return getName();}
		public double getWidth() {return this.dimension.getWidth();}
		public double getHeight() {return this.dimension.getHeight();}
		public int getLeftMargin() { return dynaParams.getIntAttribute("left-margin").orElse(100);}
		public int getRightMargin() { return dynaParams.getIntAttribute("right-margin").orElse(100);}
		public int getTopMargin() { return dynaParams.getIntAttribute("top-margin").orElse(100);}
		public int getBottomMargin() { return dynaParams.getIntAttribute("bottom-margin").orElse(100);}
		protected boolean isSex(final Locatable loc) {
			return loc.getContig().matches("(chr)?[XY]");
			}
		void adjustTracksWidth() {
			Collections.sort(this.regions,comparator);
			final long reference_length = this.regions.stream().mapToLong(R->R.getLengthOnReference()).sum();

		
			final double d1 = dynaParams.getDoubleAttribute("h-spacing-interval").orElse(1);
			final  double width = this.getWidth() - (this.regions.size() -1)* d1;

			double x = getLeftMargin();
			for(Region t:this.regions) {
				t.x = x;
				t.width = (t.getLengthOnReference()/(double)reference_length) * width;
				x += t.width;
				x += d1;
				}
			}
		
		
		void title(Object o) throws XMLStreamException {
			if(o==null) return;
			String s = String.valueOf(o);
			if(StringUtils.isBlank(s)) return;
			this.w.writeStartElement("title");
			this.w.writeCharacters(s);
			this.w.writeEndElement();
			}
		
		void beginDocument() throws XMLStreamException {
			this.w.writeStartDocument("UTF-8", "1.0");
			this.w.writeStartElement("svg");
			this.w.writeDefaultNamespace(SVG.NS);
			this.w.writeAttribute("width", format(getWidth()+getLeftMargin()+getRightMargin()));
			this.w.writeAttribute("height", format(getHeight()+getTopMargin()+getBottomMargin()));
			title(dynaParams.getAttribute("title",this.filename));
			this.w.writeStartElement("style");
			writeStyle();
			this.w.writeEndElement();
			}
		void writeStyle() throws XMLStreamException {
			w.writeCharacters(".contig {stroke:dimgray;stroke-width:0.8px;}");
			w.writeCharacters(".contig0 {fill:gainsboro;}");
			w.writeCharacters(".contig1 {fill:lightgrey;}");
			w.writeCharacters(".textContig {fill:black;stroke:none;font-size:6px;text-anchor:middle;}");
			}
		void drawContigs()  throws XMLStreamException {
			this.w.writeStartElement("g");
			for(int i=0;i<this.regions.size();i++) {
				final Region r = this.regions.get(i);
				this.w.writeStartElement("g");
				
				this.w.writeStartElement("rect");
				this.w.writeAttribute("class","contig contig"+(i%2));
				this.w.writeAttribute("x", format(r.x));
				this.w.writeAttribute("y", format(getTopMargin()));
				this.w.writeAttribute("width", format(r.width));
				this.w.writeAttribute("height", format(getHeight()));
				this.w.writeEndElement();
				
				this.w.writeStartElement("text");
				this.w.writeAttribute("class","textContig");
				this.w.writeAttribute("x", format(r.x+r.width/2.0));
				this.w.writeAttribute("y", format(getTopMargin()-10));
				w.writeCharacters(r.getName());
				w.writeEndElement();
				
				this.w.writeEndElement();
				}
			this.w.writeEndElement();
			}
		
		void drawTracks()  throws XMLStreamException {
			for(int i=0;i<this.tracks.size();i++) {
				final Track t = this.tracks.get(i);
				this.w.writeStartElement("text");
				this.w.writeAttribute("class","textTrack");
				this.w.writeAttribute("x", "0");
				this.w.writeAttribute("y", "0");
				this.w.writeAttribute("transform", "translate("+format(getLeftMargin()+getWidth())+","+format(t.y)+") rotate(90) ");
				w.writeCharacters(t.getName());
				w.writeEndElement();
				}
			}
		
		void finish() throws XMLStreamException {
		
			}
		void endDocument() throws XMLStreamException {
			this.w.writeEndElement();//svg
			this.w.writeEndDocument();
			}
		
		protected RectPainter createPainter(Pane pane) {
			return new RectPainter(pane);
			}
		
		
		protected void paint(RectPainter painter,final Locatable loc,double y0,double y1)  throws XMLStreamException {
			painter.
				href(loc2url.apply(loc)).
				title(loc2title.apply(loc)).
				style(loc2style.apply(loc)).
				clazz(loc2clazz.apply(loc)).
				shape(loc2shape.apply(loc)).
				paintY1Y2(loc, y0, y1);
			}
		protected void finish01() throws XMLStreamException {
			final List<Pane> panes = createPanes();
			
			for(final Pane pane:panes) {
				pane.pileup();
				}
			final int maxrows = panes.stream().mapToInt(P->P.rows.size()).max().orElse(0);
			if(maxrows==0) return;
			
			double d1 = dynaParams.getDoubleAttribute("h-spacing-interval").orElse(1);
			double height = this.getHeight() - (this.tracks.size() -1)* d1;
			double y1 = getTopMargin();
			for(Track t : this.tracks) {
				int nrows = panes.stream().filter(P->P.track == t).mapToInt(P->P.rows.size()).max().orElse(0);
				t.y = y1;
				t.height = (nrows/(double)maxrows) * height;
				y1 += t.height;
				y1 += d1;
				}
			/* plot each pane */
			double maxFeatureHeight =  dynaParams.getDoubleAttribute("max-feature-height").orElse(12);
			for(Pane pane:panes) {
				/* bottom y for this pane */
				y1 = pane.track.y + pane.track.height;
				final RectPainter painter  = createPainter(pane);
				/* rect height */
				final double dy= Math.min(maxFeatureHeight, pane.track.height / pane.rows.size());
				for(int y=0;y<pane.rows.size();++y) {
					y1-= dy;
					for(Locatable loc: pane.rows.get(y)) {
						paint(painter.reset(),loc, y1,y1+(dy*0.9));
						}
					}
				}
				
			}

		protected abstract CloseableIterator<? extends Locatable> open(final List<String> args)  throws IOException;
		
		/** find a track for this loc and add the item */
		abstract void visit(final Locatable loc);
		}
	
	
		private abstract class AbstractTsvHander extends Handler {
			/* header */
			protected final List<String> columns = new ArrayList<>();
			protected final Map<String,Integer> column2index = new HashMap<>();
			/** how to split lines */
			protected Pattern splitter= Pattern.compile("[\t]");
			
			/** and handler parsing TSV files */
			private class TsvIterator extends AbstractCloseableIterator<Item> {
				final BufferedReader br;
				TsvIterator(final BufferedReader br) {
					this.br = br;
					}
				@Override
				protected Item advance() {
					try {
						for(;;) {
							final String line = br.readLine();
							if(line==null) return null;
							if(StringUtils.isBlank(line)) continue;
							if(handleLine(line)) continue;
							final Item item = parseItem(line);
							if(item==null) {
								LOG.warn("cannot parse "+line+". Skipping. (handler:"+getName()+")");
								continue;
								}
							return item;
							}
						}
					catch(final IOException err ) {
						throw new RuntimeIOException(err);
						}
					}
				@Override
				public void close() {
					try{br.close();}
					catch(IOException err) {}
					}
				}
			/** register columns names and index if needed */
			protected void registerHeader(final List<String> columns) {
				this.columns.clear();
				this.column2index.clear();
				this.columns.addAll(columns);
				for(int i=0;i< columns.size();i++) this.column2index.put(columns.get(i), i);
				}
			
			/** are we able to handle this line ? */
			protected boolean handleLine(String line) {
				return false;
				}
			/** convert line to Item */
			protected Item parseItem(final String line) {
				final String[] tokens = this.splitter.split(line);
				return parseTokens(tokens);
				}
			/** convert tokens[] to Item */
			protected abstract Item parseTokens(final String[] tokens);
			
			@Override
			protected CloseableIterator<? extends Locatable> open( final List<String> args) throws IOException {
				final String input = oneFileOrNull(args);
				if(!StringUtils.isBlank(input)) {
					this.filename = input;
					}
				return new TsvIterator(openBufferedReader(input));
				}
			}
	
	
	
		private abstract class AbstractVcfHander extends Handler {
			/* vcf header*/
			protected VCFHeader vcfHeader;
			
			protected String toTitle(final VariantContext ctx) {
				return new SimpleInterval(ctx).toNiceString();
				}
			
			protected void registerHeader(final VCFHeader h) {
				this.vcfHeader = h;
				}
			
			protected VariantContext simplify(VariantContext vc) {
				return vc;
				}
			
			@Override
			protected CloseableIterator<? extends Locatable> open( final List<String> args) throws IOException {
				final String input = oneFileOrNull(args);
				if(!StringUtils.isBlank(input)) {
					this.filename = input;
					}
				final VCFIterator iter = openVCFIterator(input);
				registerHeader(iter.getHeader());
				return iter;
				}
			}
		
		private class ContigPosValueHandler1 extends AbstractTsvHander{
			@Override
			protected Item parseItem(String line) {
				return null;
				}
			@Override
			void visit(Locatable loc) {
				
				}
			@Override
			protected Item parseTokens(String[] tokens) {
				// TODO Auto-generated method stub
				return null;
			}
			}
		
		private class DellyCNVHander1 extends AbstractVcfHander {
			int CN_max=0;
			DellyCNVHander1() {
				super.loc2shape = L -> "circle";
				super.loc2title = L -> toTitle(VariantContext.class.cast(L));
				super.loc2style = L -> toStyle(VariantContext.class.cast(L));
				}

			private String toStyle(final VariantContext ctx) {
				int cn = ctx.getAttributeAsInt("CNMAX",0);
				int f=(int)((cn/(double)this.CN_max)*255);
				return "stroke:none;fill:rgba("+f+","+(255-f)+",255,0.9)";	
				}
			@Override
			protected String toTitle(final VariantContext ctx) {
				return super.toTitle(ctx)+ " CN:"+ ctx.getAttributeAsInt("CNMAX",0);	
				}
			@Override
			public String getName() {
				return "delly.cnv.1";
				}
			@Override
			protected VariantContext simplify(VariantContext vc)
				{
				int cn = vc.getGenotypes().stream().
					map(G->G.getExtendedAttribute("CN")).
					filter(S->S!=null).
					mapToInt(S->Integer.parseInt(String.valueOf(S))).
					max().orElse(0);
				if(cn==2 && !isSex(vc)) return null;
				if(cn==0) return null;
				this.CN_max = Math.max(this.CN_max , cn);
				return new VariantContextBuilder(vc).noGenotypes().
					rmAttribute("SVMETHOD").
					rmAttribute("IMPRECISE").
					attribute("CNMAX",cn).
					make();
				}
			@Override
			void visit(final Locatable loc) {
				if(!overlaps(loc)) return;
				
				VariantContext ctx = VariantContext.class.cast(loc);
				final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE,"");
				if(StringUtils.isBlank(svType)) return;
				
				final Track track = getTrackByName(svType);
				if(track==null) return;
				ctx = simplify(ctx);
				if(ctx==null) return;
				track.items.add(ctx);
				}
			@Override
			void finish() throws XMLStreamException {
				adjustTracksWidth();
				drawContigs();
				finish01();
				drawTracks();
				}
			}
	
	/** SAM sequence dict associated to the data, will be used to get the whole contigs */
	private SAMSequenceDictionary dict=null;
	private Comparator<Locatable> comparator = null;

	
	private String format(double s) {
		return String.valueOf(s);
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try	{
			this.dynaParams = AttributeMap.wrap(this.__dynaParams);
			final List<Handler> handlers = Arrays.asList(
					new DellyCNVHander1()
					);
			if(this.list_handlers) {
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					for(Handler h:handlers) {
						out.print(h.getName());
						out.print("\t");
						out.println(h.getDescription());
						}
					}
				return 0;
				}
			
			final Pattern regex = Pattern.compile(this.includeContigsRegex);
			final SAMSequenceDictionary dict0  = SequenceDictionaryUtils.extractRequired(this.faidxFile);
			this.dict = new SAMSequenceDictionary(dict0.getSequences().stream().
				filter(SR->regex.matcher(SR.getSequenceName()).matches()).
				collect(Collectors.toList()));
			this.comparator = new ContigDictComparator(this.dict).createLocatableComparator();
			final Handler handler  = handlers.
					stream().
					filter(S->S.getName().equals(this.handlerName)).
					findFirst().orElseThrow(()->new IllegalArgumentException("cannot find handler for "+this.handlerName));
			handler.dimension = this.dimension;
			
			if(this.bedPath!=null) {
				try(BedLineReader br = new BedLineReader(this.bedPath)) {
					br.setContigNameConverter(ContigNameConverter.fromOneDictionary(this.dict));
					br.setValidationStringency(ValidationStringency.LENIENT);
					while(br.hasNext()) {
						handler.newRegion(br.next());
					}
				}
			} else {
				for(SAMSequenceRecord ssr : this.dict.getSequences()) {
					final Handler.Region rgn  = handler.newRegion(ssr);
					rgn.label = ssr.getSequenceName();
					}
				}
			if(handler.regions.isEmpty()) {
				LOG.error("no region found or defined");
				return -1;
				}
		
				
			try(CloseableIterator<? extends Locatable> iter=handler.open(args)) {
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					final XMLOutputFactory xof = XMLOutputFactory.newFactory();
					handler.w = xof.createXMLStreamWriter(pw);
					handler.beginDocument();
					while(iter.hasNext()) {
						handler.visit(iter.next());
						}
					handler.finish();
					handler.endDocument();
					handler.w.flush();
					handler.w.close();
					pw.flush();
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	public static void main(final String[] args) {
		new Manhattan().instanceMainWithExit(args);
		}
	}
