package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;

public class BamCoverageTrack extends Track {
	private String reference = null;
	private String bamPath = null;
	private int featureHeight = 100;
	private Predicate<SAMRecord> samFilter = R->true;
	private int minCoverage = -1;
	public BamCoverageTrack() {
		}
	
	public void setReference(String reference) {
		this.reference = reference;
		}
	public String getReference() {
		return reference;
		}

	
	public void setBam(String path) {
		this.bamPath = path;
		}
	public String getBam() {
		return bamPath;
		}
	public void setMinCoverage(int minCoverage) {
		this.minCoverage = minCoverage;
		}
	public int getMinCoverage() {
		return minCoverage;
		}
	
	@Override
	public String getShortDesc() {
		return StringUtils.isBlank(super.shortDesc)?this.getBam():super.shortDesc;
		}
	
	public void setFeatureHeight(int featureHeight) {
		this.featureHeight = featureHeight;
		}
	public int getFeatureHeight() {
		return Math.max(10,featureHeight);
		}
	
	@Override
	public void paint(SVGContext ctx) {
		final int width = (int)ctx.image_width;
		final double featureHeight = getFeatureHeight();
		String sampleName = getBam();
		final SamReaderFactory srf = SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.LENIENT);
		if(!StringUtils.isBlank(getReference())) srf.referenceSequence(Paths.get(getReference()));
		
		try(SamReader sr=srf.open(SamInputResource.of(getBam()))) {
			final SAMFileHeader header = sr.getFileHeader();
			sampleName = header.getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtils.isBlank(S)).findFirst().orElse(getBam());
			final double[] array = new double[width];
			Arrays.fill(array, -999);
			try(CloseableIterator<SAMRecord> iter = sr.queryOverlapping(ctx.loc.getContig(), ctx.loc.getStart(), ctx.loc.getEnd())) {
				while(iter.hasNext()) {
					final SAMRecord rec = iter.next();
					if(!samFilter.test(rec) || rec.getReadUnmappedFlag()) continue;
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
			
			ctx.y+=featureHeight;
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	

}
