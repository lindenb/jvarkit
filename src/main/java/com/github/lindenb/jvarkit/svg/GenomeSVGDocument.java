/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.svg;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.util.LocatableDelegate;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

public class GenomeSVGDocument extends SVGDocument {

public class IntervalInfo extends LocatableDelegate<Locatable> {
	final String id = nextId();
	int index;//index in intervals, useful to colorize using %2==0
	double yTop = 0;
	private double width=0;
	final Element group;//group associated to this interval used for clipping
	
	public IntervalInfo(final Locatable delegate) {
		super(delegate);
		this.group = GenomeSVGDocument.this.group();
		this.group.appendChild( comment("interval "+ LocatableUtils.toNiceString(delegate)));
		}
	
	public double getWidth() {
		return this.width;
		}

	
	public double getHeight() {
		return 0;
		}
	
	GenomeSVGDocument owner() {
		return GenomeSVGDocument.this;
		}
	
	public double length2pixel(int w) {
		return ((w)/(double)getLengthOnReference())* this.getWidth();
		}

	
	public double pos2pixel(int pos) {
		return ((pos-getStart())/(double)getLengthOnReference())* this.getWidth();
		}
	public int trimPos(int pos) {
		return Math.max(getStart(), Math.min(pos, getEnd()));
		}
	public double trimPixel(double pix) {
		return Math.max(0, Math.min(getWidth(), pix));
		}
	
	Element rect(Locatable r,double y,double height,Map<String,Object> atts) {
		final double x0 = pos2pixel(trimPos(r.getStart()));
		final double x1 = pos2pixel(trimPos(r.getEnd()+1));
		return owner().rect( x0, y, (x1-x0), height, atts);
		}
	Element rect(Locatable r,double y,double height) {
		return rect(r,y,height,Collections.emptyMap());
		}
	
	private int bestTicks(final int max) {
		if(max<=10) return 1;
		if(max<=100) return 10;
		if(max<=1000) return 100;
		final int ndigit=(int)Math.ceil(Math.log10(max-1));
		return Math.max(1,(int)Math.pow(10, ndigit-2));
		}

	void finish() {
		//background-color
		final Element background = rect(this,0,this.getHeight(),Maps.of("class", style2class("stroke:none;fill:"+(this.index==0?"white":"gray")+";")));
		this.group.insertBefore(background, this.group.getFirstChild());
		
		// add clip
		final String id = nextId();
		
		
		final Element clipPath = element("clipPath",Maps.of("id",id));
		owner().defsElement.appendChild(clipPath);
		clipPath.appendChild(
				createShape().
					moveToR(0, -1).
					horizonal(getWidth()).
					make()
				);
		final Element t = element("text");
		this.group.appendChild(t);
		final Element tp = element("textPath",toNiceString(),Maps.of("href","#"+id));
		setTitle(tp,toNiceString());
		t.appendChild(tp);
		
		this.group.insertBefore(anchor(tp,this),  this.group.getFirstChild());
		}
	
	public void plotHightlights() {
		if(owner().highlights.isEmpty()) return;
		final Element highgroup = group(Maps.of("class", style2class("highlight","stroke:none;fill:pink;opacity:0.4")));
		for(Interval rgn: owner().highlights) {
			if(!rgn.overlaps(this)) continue;
			final Element rect = rect(rgn,0,this.getHeight());
			setTitle(rect,rgn.getName());
			highgroup.appendChild(rect);
			}
		this.group.insertBefore(highgroup, this.group.getFirstChild());
		}
	
	public Element insertRuler() {
		/** vertical ruler */
		final Element ruler_gh = group();
		
		final int sep= bestTicks(this.getLengthOnReference());
		int pos = this.getStart() + this.getStart()%sep;
		
		while(pos<= this.getEnd())
			{
			double x= pos2pixel(pos);
			final Element line = line(x,x,this.yTop,this.yTop,Maps.of("class", style2class("ruler","stroke:lightgray;stroke-width:0.5px")));
			callbacks_on_finish.add(()->line.setAttribute("y2", String.valueOf(GenomeSVGDocument.this.lastY)));
			setTitle(line,StringUtils.niceInt(pos));
			ruler_gh.appendChild(line);
			
			final Element label = text(
					0,0,
					StringUtils.niceInt(pos),
					Maps.of(
						"class", style2class("rulerlabel","stroke:gray;stroke-width:0.5px;font-size:7px;"),
						"transform", translate(x,this.yTop)+" rotate(90) "
					));
			ruler_gh.appendChild(label);
			
			pos+=sep;
			}
		return ruler_gh;
		}
	}

private final SAMSequenceDictionary dict;
private final AttributeMap properties;
private final List<IntervalInfo> intervals;
private final IntervalTreeMap<IntervalInfo> intervalsTreeMap;
private final double image_width;
private final double margin_left;
private final List<Runnable> callbacks_on_finish = new ArrayList<>();
private double lastY = 0;
private final List<Interval> highlights = new ArrayList<>();

public GenomeSVGDocument(
		final SAMSequenceDictionary dict,
		final List<Interval> intervals,
		final AttributeMap properties
		) {
	this.dict = Objects.requireNonNull(dict);
	if(Objects.requireNonNull(intervals).isEmpty()) throw new IllegalArgumentException();
	this.properties=properties;
	
	final Comparator<Locatable> sorter = new ContigDictComparator(this.dict).createLocatableComparator();
	this.intervals = intervals.stream().
		sorted(sorter).
		map(R->new IntervalInfo(R)).
		collect(Collectors.toList());
	this.intervalsTreeMap = new IntervalTreeMap<>();
	for(IntervalInfo ii: this.intervals) {
		this.intervalsTreeMap.put(ii.toInterval(),ii);
		}
	this.image_width =  this.properties.getDoubleAttribute("image-width").orElse(700);
	this.margin_left =  this.properties.getDoubleAttribute("margin-left").orElse(this.image_width/10.0);
	double spaceBetweenInterval = this.properties.getDoubleAttribute("space-between-regions").orElse(1);
	final double margin_right =  this.properties.getDoubleAttribute("margin-right").orElse(10);
	
	double adjust_image_width = this.image_width - Math.max(0,(intervals.size()-1)*spaceBetweenInterval);
	final long sum_length_on_ref = this.intervals.stream().mapToLong(R->R.getLengthOnReference()).sum();
	double x = this.margin_left;
	for(int i=0;i< this.intervals.size();++i) {
		if(i>0) x+=spaceBetweenInterval;
		final IntervalInfo ii = this.intervals.get(i);
		ii.index=i;
		ii.group.setAttribute("transform",translate(x,this.lastY));//TODO shit top
		ii.width = (ii.getLengthOnReference()/(double)sum_length_on_ref)* adjust_image_width;
		x+= ii.width;
		}
	
	setWidth(this.margin_left+this.image_width+margin_right+1);
	
	
	this.callbacks_on_finish.add(()->{
	
		});
	}

public void finish() {
	for(IntervalInfo ii: this.intervals) {
		ii.finish();
		}
	
	// draw black box at the end
	double w = Double.parseDouble(svgElement.getAttribute("width")) -1;
	Element g = rect(0,0,w,this.lastY,Maps.of("class",style2class("frame","fill:none;stroke:black;")));
	svgElement.appendChild(g);
	g = rect(0,0,this.margin_left,this.lastY,Maps.of("class",style2class("frame","fill:none;stroke:black;")));
	svgElement.appendChild(g);
	setHeight(lastY+1);
	
	
	callbacks_on_finish.stream().forEach(R->R.run());
	}

public String getUrl(final Locatable loc) {
	return Hyperlink.compile(this.dict).apply(loc).orElse("");
	}

public Element anchor(Element wrapped, Locatable loc) {
	return this.anchor(wrapped, getUrl(loc));
	}

public void insertRuler() {
	/** vertical ruler */
	final Element ruler_gv = group();
	for(IntervalInfo ii:this.intervals) {
		final Element g = ii.insertRuler();
		ruler_gv.appendChild(g);
		}

	}

/** return true if any IntervalInfo overlaps loc */
public boolean overlaps(final Locatable loc) {
	return this.intervals.stream().anyMatch(R->R.overlaps(loc));
	}

}
