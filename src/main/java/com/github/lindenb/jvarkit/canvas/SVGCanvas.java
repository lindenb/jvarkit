package com.github.lindenb.jvarkit.canvas;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.PathIterator;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.Map;
import java.util.Objects;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.util.RuntimeIOException;

public class SVGCanvas extends Canvas {
	private final int width;
	private final int height;
	private final OutputStream ostream;
	private final XMLStreamWriter w;
	public SVGCanvas(Path out,int width,int height,boolean compressed,final FunctionalMap<String,Object> params) throws IOException,XMLStreamException {
		super();
		this.width=width;
		this.height=height;
				
		OutputStream base;
		
		if(out==null) {
			base = System.out;
			}
		else
			{
			base = Files.newOutputStream(out);
			}
		if(compressed) {
			base = new GZIPOutputStream(base);
			}
		this.ostream=base;
		w = XMLOutputFactory.newFactory().
				createXMLStreamWriter(ostream, "UTF-8");
		w.writeStartDocument("UTF-8", "1.0");
		w.writeStartElement("svg");
		w.writeDefaultNamespace(SVG.NS);
		w.writeAttribute("version", "1.1");
		w.writeAttribute("width", String.valueOf(width));
		w.writeAttribute("height", String.valueOf(height));
		
		
		w.writeStartElement("title");
		w.writeCharacters(String.valueOf(params.getOrDefault(KEY_TITLE,"")));
		w.writeEndElement();
		
		
		
		final Object meta= params.get(KEY_SVG_META);
		if(meta!=null && meta instanceof Consumer) {
			w.writeStartElement("meta");
			@SuppressWarnings("unchecked")
			Consumer<XMLStreamWriter> cons=(Consumer<XMLStreamWriter>)meta;
			cons.accept(this.w);
			w.writeEndElement();
			}
		
		
		
		w.writeStartElement("g");
		writeStyle();
		begin(params.minus(KEY_SVG_META).minus(KEY_TITLE));
		}
	
	private String toString(Object v) {
		if(v ==null) throw new IllegalArgumentException();
		if(v instanceof Color) return ColorUtils.toRGB(Color.class.cast(v));
		return String.valueOf(v);
		}
	
	private Map.Entry<String,String> toString(Map.Entry<String,Object> h) {
		String v;
		if((h.getKey().equals(KEY_FILL) || h.getKey().equals(KEY_STROKE)) && h.getValue()==null) {
			v = "none";
			}
		else if((h.getKey().equals(KEY_FONT_SIZE)  && (h.getValue() instanceof Number))) {
			v = toString(h.getValue())+"px";
			}
		else
			{
			v = toString(h.getValue());
			}
		return new AbstractMap.SimpleEntry<String,String>(h.getKey(),v);
		}
	
	private FunctionalMap<String, Object> whatsNew() {
		FunctionalMap<String, Object> fm1 = super.states.peek();
		if(super.states.size()<2) return fm1;
		FunctionalMap<String, Object> ret = FunctionalMap.make();
		final FunctionalMap<String, Object> fm0 = super.states.get(super.states.size()-2);
		for(final String key:fm1.keySet()) {
			Object v1= fm1.get(key);
			if(!fm0.containsKey(key)) {
				ret=ret.plus(key,v1);
				continue;
				}
			Object v0 = fm0.get(key);
			if(Objects.equals(v1, v0)) continue;
			ret=ret.plus(key,v1);
			}
		return ret;
		}
		
	
	@Override
	public Canvas begin(final FunctionalMap<String, Object> fm) {
		super.begin(fm);
		try {
			this.w.writeStartElement("g");
			}
		catch(final XMLStreamException err) {
			throw new RuntimeIOException(err);
			}
		return this;
		}
	
	@Override
	public Canvas end() {
		try {
			this.w.writeEndElement();//g
			}
		catch(final XMLStreamException err) {
			throw new RuntimeIOException(err);
			}
		return super.end();
		}
	
	@Override
	public Canvas circle(double cx, double cy, double r, FunctionalMap<String, Object> fm) {
		this.states.push(this.states.peek().plus(fm));
		try {
			beginAnchor();
			this.w.writeStartElement("circle");
			this.w.writeAttribute("cx", round(cx));
			this.w.writeAttribute("cy", round(cy));
			this.w.writeAttribute("r", round(r));
			writeStyle();
			inner();
			this.w.writeEndElement();
			endAnchor();
			}
		catch(final XMLStreamException err) {
			throw new RuntimeIOException(err);
			}
		this.states.pop();
		return this;
		}
	
	@Override
	public Canvas rect(double x, double y, double width, double height, FunctionalMap<String, Object> fm) {
		this.states.push(this.states.peek().plus(fm));
		try {
			beginAnchor();
			this.w.writeStartElement("rect");
			this.w.writeAttribute("x", round(x));
			this.w.writeAttribute("y", round(y));
			this.w.writeAttribute("width", round(width));
			this.w.writeAttribute("height", round(height));
			writeStyle();
			inner();
			this.w.writeEndElement();
			endAnchor();
			}
		catch(final XMLStreamException err) {
			throw new RuntimeIOException(err);
			}
		this.states.pop();
		return this;
		}
	@Override
	public Canvas text(String s, double x, double y, FunctionalMap<String, Object> fm) {
		if(StringUtils.isBlank(s)) return this;
		this.states.push(this.states.peek().plus(fm));
		try {
			beginAnchor();
			this.w.writeStartElement("text");
			this.w.writeAttribute("x",round(x));
			this.w.writeAttribute("y",round(y));
			writeStyle();
			w.writeCharacters(s);
			inner();
			this.w.writeEndElement();
			endAnchor();
			
			this.w.flush();
			}
		catch(final XMLStreamException err) {
			throw new RuntimeIOException(err);
			}
		this.states.pop();
		return this;
		}
	
	
	private String round(double x) {
		String s= String.format("%.2f",x);
		if(s.endsWith(".00")) s=s.substring(0,s.length()-3);
		return s;
	}
	
	@Override
	public Canvas shape(Shape shape, final FunctionalMap<String, Object> fm) {
		this.states.push(this.states.peek().plus(fm));

		final StringBuilder path = new StringBuilder();
		final double tab[] = new double[6];
		final PathIterator pathiterator = shape.getPathIterator(null);

		while(!pathiterator.isDone())
			{
				int currSegmentType= pathiterator.currentSegment(tab);
				switch(currSegmentType) {
				case PathIterator.SEG_MOVETO: {
					path.append( "M " + round(tab[0]) + " " + round(tab[1]) + " ");
					break;
				}
				case PathIterator.SEG_LINETO: {
					path.append( "L " + round(tab[0]) + " " + round(tab[1]) + " ");
					break;
				}
				case PathIterator.SEG_CLOSE: {
					path.append( "Z ");
					break;
				}
				case PathIterator.SEG_QUADTO: {
					path.append( "Q " + round(tab[0]) + " " + round(tab[1]));
					path.append( " "  + round(tab[2]) + " " + round(tab[3]));
					path.append( " ");
					break;
				}
				case PathIterator.SEG_CUBICTO: {
					path.append( "C " + round(tab[0]) + " " + round(tab[1]));
					path.append( " "  + round(tab[2]) + " " + round(tab[3]));
					path.append( " "  + round(tab[4]) + " " + round(tab[5]));
					path.append( " ");
					break;
				}
				default:
				{
					throw new IllegalStateException("Cannot handled "+currSegmentType);
				}
				}
			pathiterator.next();
			}
		
		try {
			beginAnchor();
			this.w.writeStartElement("path");
			this.w.writeAttribute("d",path.toString());
			writeStyle();
			inner();
			this.w.writeEndElement();
			endAnchor();
			this.w.flush();
			}
		catch(final XMLStreamException err) {
			throw new RuntimeIOException(err);
			}
		this.states.pop();
		return this;
		}
	
	@Override
	public Canvas comment(String s) {
		if(!StringUtils.isBlank(s)) {
			try {
				w.writeComment(s);
				}
			catch(XMLStreamException err) {
				throw new RuntimeIOException(err);
				}
			}
		return this;
		}
		
	private void writeStyle() throws XMLStreamException {
		final FunctionalMap<String, Object> fm = whatsNew();
		final String css= fm.
			stream().
			filter(KV->!(
					KV.getKey().equals(KEY_HREF) ||
					KV.getKey().equals(KEY_TITLE) || 
					KV.getKey().equals(KEY_SVG_META))).
			map(KV->toString(KV)).
			map(KV->KV.getKey()+":"+toString(KV.getValue())).
			collect(Collectors.joining(";"));
		if(css.isEmpty()) return;
		w.writeAttribute("style", css.toString());
		}
	
	private void beginAnchor() throws XMLStreamException {
		final String href = whatsNew().getOrDefault(KEY_HREF,"").toString();
		if(!StringUtils.isBlank(href)) {
			w.writeStartElement( "a");
			w.writeAttribute("href", href);
			w.writeAttribute("target", "_blank");
			}
	}
	private void endAnchor()  throws XMLStreamException {
		final String href = whatsNew().getOrDefault(KEY_HREF,"").toString();
		if(!StringUtils.isBlank(href)) {
			w.writeEndElement();
			}
		}
	
	private void inner()  throws XMLStreamException  {
		FunctionalMap<String, Object> fm = whatsNew();
		Object o = fm.get(KEY_TITLE);
		if(o!=null) {
			String s=String.valueOf(o);
			if(!StringUtils.isBlank(s)) {
				w.writeStartElement("title");
				w.writeCharacters(s);
				w.writeEndElement();
				}
			}
		}
	
	@Override
	public void close() throws IOException {
		try {
			this.w.writeComment("Generated with jvarkit. Author: Pierre Lindenbaum Phd.");
			end();
			this.w.writeEndElement();//g
			this.w.writeEndElement();//svg
			this.w.writeEndDocument();
			this.w.close();
			this.ostream.flush();
			this.ostream.close();
			}
		catch(final XMLStreamException err) {
			throw new IOException(err);
			}
		}
	@Override
	public int getWidth() {
		return width;
		}
	@Override
	public int getHeight() {
		return height;
		}
	}
