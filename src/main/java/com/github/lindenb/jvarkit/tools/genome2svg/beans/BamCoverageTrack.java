package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;

public class BamCoverageTrack extends AbstractBamTrack {
	private int minCoverage = -1;
	public BamCoverageTrack() {
		setFeatureHeight(100);
		}
	
	public void setMinCoverage(int minCoverage) {
		this.minCoverage = minCoverage;
		}
	public int getMinCoverage() {
		return minCoverage;
		}
	
	@Override
	public String getShortDesc() {
		String s = super.getShortDesc();
		if(StringUtils.isBlank(s)) {
			if(getPaths().size()==1) {
				s = getPaths().get(0).getShortDesc();
				}
			else
				{
				s = String.valueOf(getPaths().size())+" Bams";
				}
			}
		return s;
		}
	
	private  class Fragment implements Locatable, Comparable<Fragment >{
		private SAMRecord rec1 = null;
		private SAMRecord rec2 = null;
		Fragment(final SAMRecord rec1) {
			this.rec1 = rec1;
			}
		
		@Override
		public int compareTo(Fragment o) {
			return Integer.compare(getStart(), o.getStart());
			}
		
		SAMRecord getFirst() {
			return this.rec1;
			}
		SAMRecord getSecond() {
			return this.rec1==null?getFirst():this.rec2;
			}
		
		@Override
		public String getContig() {
			return getFirst().getContig();
			}
		@Override
		public int getStart() {
			return Math.min(getFirst().getStart(), getSecond().getStart());
			}
		
		String getTitle() {
			return rec1.getReadName()+" "+rec1.getContig()+":"+getStart()+"-"+getEnd();
			}
		
		private int getMateEnd(final SAMRecord rec) {
			if(!rec.getReadPairedFlag()) return -1;
			if(rec.getReadUnmappedFlag()) return -1;
			int n = rec.getMateAlignmentStart();
			if(rec.hasAttribute(SAMTag.MC)) n = Math.max(n, SAMUtils.getMateAlignmentEnd(rec));
			return n;
			}
		
		@Override
		public int getEnd() {
			int n = this.rec1.getEnd();
			if(this.rec2!=null) {
				Math.max(rec2.getEnd(),n);
				}
			else
				{
				n  = Math.max(getMateEnd(this.rec1),n);
				}
			return n;
			}
		
		void put(SAMRecord rec) {
			if(this.rec2!=null) return;
			if(!rec1.getReadPairedFlag()) return;
			if(rec1.getMateUnmappedFlag()) return;
			if(!this.rec1.getReadName().equals(rec.getReadName())) return;
			if(!rec.getReadPairedFlag()) return;
			if(rec.getFirstOfPairFlag()==this.rec1.getFirstOfPairFlag()) return;
			if(rec.getSecondOfPairFlag()==this.rec1.getSecondOfPairFlag()) return;
			if(this.rec1.getStart() < rec2.getStart()) {
				this.rec2=rec;
				}
			else
				{
				this.rec2=this.rec1;
				this.rec1=rec;
				}
			}
		
		
		Element display(SVGContext ctx,double height) {
			final Element g= ctx.element("g");
			boolean require_line=true;
			
			if(rec2!=null) {
				if(rec1.getEnd() >= rec2.getStart() || (int)ctx.pos2pixel(rec1.getEnd()) >= (int)ctx.pos2pixel(rec2.getStart()) ) {
					require_line = false;
					Element rect12 = ctx.rect(
							Math.min(rec1.getStart(),rec2.getStart()),
							Math.max(rec1.getEnd(),rec2.getEnd()),
							ctx.y,
							height
							);
					g.appendChild(rect12);
					}
				else
					{
					Element rect = ctx.rect(rec1.getStart(),rec1.getEnd(),ctx.y,height);
					g.appendChild(rect);
					
					rect = ctx.rect(rec2.getStart(),rec2.getEnd(),ctx.y,height);
					g.appendChild(rect);
					}
				}
			else
				{
				require_line = false;
				Element rect1 = ctx.rect(rec1.getStart(),rec1.getEnd(),ctx.y,height);
				g.appendChild(rect1);
				}
			
			
			
			if(require_line) {
				Element line = ctx.line(this,ctx.y);
				g.insertBefore(line, g.getFirstChild());
				}
			g.appendChild(ctx.title(getTitle()));
			return ctx.anchor(g, getUrlForLocatable(this));
			}
		}
	private boolean paintLowResolution(SVGContext ctx,double featureHeight) {
		for(PathBean bamBean:getPaths()) {
			paintLowResolution( BamBean.class.cast(bamBean), ctx,featureHeight);
			}
		return true;
		}
	
	private Element paintLowResolution(BamBean bamBean,SVGContext ctx, double featureHeight) {
		final long max_reads = 1000L;
		final Map<String,Fragment> readName2frag  = new HashMap<>();
		try(SamReader sr= bamBean.openSamReader()) {
			final Locatable loc = bamBean.resolveContig(ctx.loc);
			if(loc!=null) {
				try(CloseableIterator<SAMRecord> iter = sr.queryOverlapping(loc.getContig(), loc.getStart(), loc.getEnd())) {
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						Fragment frag = readName2frag.get(rec.getReadName());
						if(frag==null) {
							frag = new Fragment(rec);
							readName2frag.put(rec.getReadName(), frag);
							if(readName2frag.size() > max_reads) return null;
							}
						frag.put(rec);
						}
					}
				}
		} catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		final Pileup<Fragment> pileup = new Pileup<>(ctx.createCollisionPredicate());
		pileup.addAll(readName2frag.
				values().
				stream().
				filter(F->(getSamFilter().test(F.getFirst())  || getSamFilter().test(F.getSecond()) )).
				sorted().
				collect(Collectors.toList())
				);
		Element g0 = ctx.element("g");
		for(List<Fragment> row: pileup.getRows()) {
			ctx.y+=1;
			for(Fragment frag : row) {
				g0.appendChild(frag.display(ctx, featureHeight));
				}
			ctx.y+= featureHeight;
			ctx.y+=1;
			}
		return g0;
		}
	
	private int maxCoverage(BamBean bamBean, Locatable interval) {
		final Locatable loc = bamBean.resolveContig(interval);
		if(loc==null) return 0;
		long maxCov=0L;
		try(SamReader sr= bamBean.openSamReader()) {
			final IntervalList il = new IntervalList(sr.getFileHeader());
			il.add(new Interval(loc));
			try(SamLocusIterator iter = new SamLocusIterator(sr, il,true)) {
				while(iter.hasNext()) {
					final SamLocusIterator.LocusInfo locusInfo = iter.next();
					maxCov = Math.max(maxCov,
							locusInfo.getRecordAndOffsets().
							stream().
							filter(LI->getSamFilter().test((LI.getRecord()))).
							count()
							)
							;
					}
				}
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		return (int)maxCov;
		}
	
		private boolean paintCoverage(SVGContext ctx) {
			final int width = (int)ctx.image_width;
			final double featureHeight =  Math.max(10,getFeatureHeight());
			String sampleName;
			
			int userMaxCov=0;
			for(PathBean bamBean: getPaths()) {
				userMaxCov = Math.max(maxCoverage(BamBean.class.cast(bamBean), ctx.loc), userMaxCov);
				}
			
			final double[] array = new double[width];
	
			for(PathBean bamBean: getPaths()) {
				try(SamReader sr= BamBean.class.cast(bamBean).openSamReader()) {
					sampleName = BamBean.class.cast(bamBean).getSampleName();
					Arrays.fill(array, -999);
					try(CloseableIterator<SAMRecord> iter = sr.queryOverlapping(ctx.loc.getContig(), ctx.loc.getStart(), ctx.loc.getEnd())) {
						while(iter.hasNext()) {
							final SAMRecord rec = iter.next();
							if(!getSamFilter().test(rec) || rec.getReadUnmappedFlag()) continue;
							int prev=-1;
							for(AlignmentBlock block : rec.getAlignmentBlocks()) {
								for(int x=0;x< block.getLength();++x) {
									int p1 = block.getReferenceStart()+x;
									if(p1 < ctx.loc.getStart()) continue;
									if(p1 > ctx.loc.getEnd()) break;
									int p2 = (int)ctx.pos2pixel(p1);
									if(prev==p2) continue;
									prev=p2;
									if(p2<0 || p2>array.length) continue;
									if(array[p2]<0) array[p2]=0;
									array[p2]++;
									}
								}
							}
						}
					int x=0;
					while(x<array.length) {
						if(array[x]<0) {
							array[x] = x>0?array[x-1]:0;
							}
						x++;
						}
					double maxcov = Math.max(1, Arrays.stream(array).max().orElse(1));
					if(getMinCoverage()>0 && maxcov < getMinCoverage()) maxcov = getMinCoverage();
					//defs
					{
					Element grad= ctx.element("linearGradient");
					ctx.defsNode.appendChild(grad);
					grad.setAttribute("id",getId()+"grad");
					grad.setAttribute("x1","0%");
					grad.setAttribute("y1","0%");
					grad.setAttribute("x2","0%");
					grad.setAttribute("y2","100%");
					//grad.setAttribute("gradientTransform","rotate(90)");
					
					Element stop = ctx.element("stop");
					grad.appendChild(stop);
					stop.setAttribute("offset", "0%");
					stop.setAttribute("stop-color",(maxcov>50?"red":maxcov>20?"green":"blue"));
					
					stop = ctx.element("stop");
					grad.appendChild(stop);
					stop.setAttribute("offset", "100%");
					stop.setAttribute("stop-color", "darkblue");
					}
					
					
					final Element g0 = ctx.element("g");
					ctx.tracksNode.appendChild(g0);
					insertTitle(g0,ctx);
					
					final Element legend = ctx.element("text",sampleName);
					legend.setAttribute("x", format(-ctx.left_margin));
					legend.setAttribute("y",String.valueOf(ctx.y+10));
					g0.appendChild(legend);
					
					final Element g = ctx.element("g");
					g0.appendChild(g);
					g.setAttribute("transform", "translate(0,"+ctx.y+")");
					ctx.tracksNode.appendChild(g);
					
					final Element path = ctx.element("polyline");
					path.setAttribute("style", "stroke:red;fill:url('#"+getId()+"grad')");
					g.appendChild(path);
					final List<Point2D> points = new ArrayList<>(array.length+2);
					points.add(new Point2D.Double(0,featureHeight));
					for(int i=0;i< array.length;i++) {
						points.add(new Point2D.Double(i, featureHeight - (array[i]/maxcov)*featureHeight));
						}
					points.add(new Point2D.Double(array.length-1,featureHeight));
					points.add(new Point2D.Double(0,featureHeight));
					x=1;
					while(x<points.size()) {
						if(x>0 && x+1 < points.size() &&  points.get(x-1).getY()==points.get(x).getY() && points.get(x).getY()==points.get(x+1).getY()) {
							points.remove(x);
							}
						else
							{
							++x;
							}
						}
									
					path.setAttribute("points", points.stream().
							map(P->format(P.getX())+","+format(P.getY())).
							collect(Collectors.joining(" "))
							);
					// yaxis
					ctx.styleNode.appendChild(ctx.text("."+getId()+"ticks{stroke:darkgray;stroke-width:1px;}\n"));
					ctx.styleNode.appendChild(ctx.text("."+getId()+"lbl{text-anchor:end;font-size:10px;}\n"));
					final int nTicks = 10;
					for(int i=0;i< nTicks;i++) {
						final double cov = (maxcov/nTicks)*i;
						final double v = featureHeight - (cov/maxcov) * featureHeight;
						final Element line = ctx.element("line");
						line.setAttribute("x1", format(0));
						line.setAttribute("y1", format(v));
						line.setAttribute("x2", format(-5));
						line.setAttribute("y2", format(v));
						line.setAttribute("class",getId()+"ticks");
	
						g.appendChild(line);
						
						final Element text = ctx.element("text",format(cov));
						text.setAttribute("class",getId()+"lbl");
						text.setAttribute("x", format(-5));
						text.setAttribute("y", format(v));
						g.appendChild(text);
						
					}
					
				}
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
				ctx.y+=featureHeight;
			}
				
		return true;
		}
	
	@Override
	public void paint(SVGContext ctx) {
		
		
		double save_y = ctx.y;
		if(paintLowResolution(ctx,10)) {
			return;
			}
		ctx.y = save_y;
		paintCoverage(ctx);
		}
}