package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.MinMaxDouble;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;
import com.github.lindenb.jvarkit.wig.BigWigReader;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;

public class WiggleTrack extends Track {
	private Double minValue = null;
	private Double maxValue = null;
	private boolean normalize  = false;
	private boolean packed = false;
	public WiggleTrack() {
		setFeatureHeight(100);
		}
	
	public void setMinValue(double v) {
		this.minValue = v;
		}
	public void setMaxValue(double v) {
		this.maxValue = v;
		}
	
	public boolean isNormalize() {
		return normalize;
		}

	public void setNormalize(boolean normalize) {
		this.normalize = normalize;
		}
	
	public void setPacked(boolean packed) {
		this.packed = packed;
		}
	public boolean isPacked() {
		return packed && getPaths().size()>1;
		}
	
	@Override
	public void paint(SVGContext ctx) {
		if(getPaths().isEmpty()) {
			getLogger().warn("No Path for "+getId());
			return;
			}
		final double featureHeight = getFeatureHeight();
		double all_bigwig_max_value = 0.0;
		final List<? extends PathBean> wigglePaths = getPaths();
		final List<Double> max_values  = new ArrayList<>(wigglePaths.size());
		final Element g0 = ctx.element("g");
		for(int side=0; side < 2; ++side ) {
			/* plot each file */
			for(int wigIdx=0;wigIdx< wigglePaths.size();++wigIdx) {
				final PathBean path = wigglePaths.get(wigIdx);
				
				getLogger().debug("scanning "+path.getPath());
				final double[] array = new double[(int)ctx.image_width];
				final int[] count = new int[array.length];
				try(BigWigReader bbr = new BigWigReader(path.asPath())) {
					
					try(CloseableIterator<BigWigReader.WigItem> iter=bbr.query(ctx.loc)) {
						while(iter.hasNext()) {
							final BigWigReader.WigItem item=iter.next();
							if(!item.overlaps(ctx.loc)) continue;//paranoid
							final float v=item.getValue();
							int min_array_x = Math.max(0, (int)ctx.pos2pixel(item.getStart())-1);
							int max_array_x = Math.min(array.length, (int)ctx.pos2pixel(item.getEnd())+1);
							
							for(int x=min_array_x; x < max_array_x;++x) {
								int pos  = (int)ctx.pixel2genomic(x);
								if(pos< item.getStart() || pos > item.getEnd()) continue;
								count[x]++;
								array[x] += v;
								}
							}
						}
					
					for(int i=0; i < count.length;i++) {
						if(count[i]>1) array[i]/=count[i];
						}
					
					
					if(side==0) {
						double wig_max_value = Arrays.stream(array).max().orElse(1.0);
						if(this.minValue!=null && this.minValue.doubleValue() > wig_max_value) {
							wig_max_value = this.minValue.doubleValue();
							}
						if(this.maxValue!=null && this.maxValue.doubleValue() < wig_max_value) {
							wig_max_value = this.maxValue.doubleValue();
							}
						getLogger().debug("max value:"+wig_max_value);

						max_values.add(wig_max_value);
						}
					else
						{
						final double wig_max_value = (isNormalize()?all_bigwig_max_value:max_values.get(wigIdx));
						
						final Element g = ctx.element("g");
						g0.appendChild(g);
						
						final String style = "p"+path.getId();
						ctx.appendStyle(
							"."+style+"{stroke:"+(
										isPacked()?getColorPalette(wigIdx/(float)wigglePaths.size()):StringUtils.ifBlank(path.getStroke(),"gray"))+";fill:" +
										(
										isPacked()?
										"none":
										StringUtils.ifBlank(path.getFill(),getColorPalette(wigIdx/(float)wigglePaths.size()))
										) +
									";stroke-width:1px;}\n" +
									"."+style+"ti {stroke:none;fill:darkgray;text-anchor:middle;font-size:10px;}\n" +
									"."+style+"lbl {stroke:none;fill:darkgray;text-anchor:end;font-size:10px;}\n"
							);
						if(!isPacked() || wigIdx==0) {
							final Element title = ctx.element("text",path.asPath().getFileName().toString());
							title.setAttribute("class",style+"ti");
							title.appendChild(ctx.title(path.asPath().toString()));
							title.setAttribute("x", format(ctx.image_width/2.0));
							title.setAttribute("y", format(ctx.y+12));
							g.appendChild(title);
							
							final Element text = ctx.element("text",path.getShortDesc());
							text.setAttribute("class",style+"lbl");
							text.setAttribute("x",format(-2));
							text.setAttribute("y",format(ctx.y+12));
							g.appendChild(text);
								
							ctx.y+=12;
							}
						
						final List<Point2D.Double> points= new ArrayList<>(count.length+2);
						points.add(new Point2D.Double(0, ctx.y+featureHeight));
						for(int i=0;i< count.length;++i) {
							double v = array[i];
							if(this.maxValue!=null && v>this.maxValue.doubleValue()) v= this.maxValue.doubleValue();
							final double h = ctx.y + featureHeight - (v/wig_max_value)*featureHeight;
							points.add(new Point2D.Double(i, h));
							}
						points.add(new Point2D.Double(count.length, ctx.y+featureHeight));
						
						int i=1;
						while(i < points.size()) {
							if(i+2  < points.size() &&
									points.get(i  ).getY()==points.get(i+1).getY() &&
									points.get(i+1).getY()==points.get(i+2).getY()) {
								points.remove(i+1);
								}
							else
								{
								i++;
								}
							}
						
						g.appendChild(ctx.comment("BEGIN BIGWIG "+ path));
		
						
						
						
						final Element polygon = ctx.element("polygon");
						polygon.setAttribute("class", style);
						polygon.setAttribute("points",
								points.stream().
								map(P->format(P.getX())+","+format(P.getY()))
								.collect(Collectors.joining(" ")));
						polygon.appendChild(ctx.title(path.getLongDesc()));
						g.appendChild(polygon);
		
						if(!isPacked() || wigIdx==0) {
							g.appendChild(ctx.createYAxis(
									0,
									ctx.y,
									0,
									wig_max_value,
									featureHeight
									));
							}
						
						if(!isPacked() || wigIdx+1== wigglePaths.size()) {
							ctx.y += featureHeight;
							ctx.y += 1;
							}
						}
					}
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
				}//end for each bigwig
			
			if(side==0) {
				all_bigwig_max_value = max_values.stream().mapToDouble(V->V.doubleValue()).max().orElse(1.0);
				}
			
			} //end for side
		ctx.tracksNode.appendChild(g0);
		}
		
	

	
	
}
