package com.github.lindenb.jvarkit.canvas;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.nio.file.Path;

import javax.xml.stream.XMLStreamException;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.FunctionalMap;

public class CanvasFactory {
	private int width=700;
	private int height=400;
	public CanvasFactory() {
		}
	
	public CanvasFactory setWidth(int width) {
		this.width = width;
		return this;
		}
	
	public CanvasFactory setHeight(int height) {
		this.height = height;
		return this;
		}
	
	public CanvasFactory setDimension(int width,int height) {
		return setWidth(width).setHeight(height);
		}
	public Canvas open(Path output) throws IOException {
		if(output==null) throw new IllegalArgumentException("path cannot be null");
		String fname=output.getFileName().toString().toLowerCase();
		if(fname.endsWith(".png")) return open(output,"png");
		if(fname.endsWith(".jpg")) return open(output,"jpg");
		if(fname.endsWith(".jpeg")) return open(output,"jpg");
		if(fname.endsWith(".svg")) return open(output,"svg");
		if(fname.endsWith(".svg.gz")) return open(output,"svg.gz");
		if(fname.endsWith(".ps") || fname.endsWith(".eps")) return open(output,"ps");
		if(fname.endsWith(".ps.gz") || fname.endsWith(".eps.gz")) return open(output,"ps.gz");
		throw new IllegalArgumentException("Cannot find output format for "+output);
		}
	
	public Canvas open(Path output, String format) throws IOException {
		if(StringUtils.isBlank(format)) {
			return this.open(output);
			}
		if(format.startsWith(".")) return open(output,format.substring(1));
		format=format.toLowerCase();
		if(format.equals("png")) {
			return new Graphics2DCanvas(output, width, height,format,FunctionalMap.make());
			}
		if(format.equals("svg") || format.equals("svg.gz")) {
			try {
				return new SVGCanvas(output, width, height,format.endsWith(".gz"),FunctionalMap.make());
				}
			catch(XMLStreamException err) {
				throw new IOException(err);
				}
			}
		
		if(format.equals("ps") || format.equals("ps.gz")) {
			return new PSCanvas(output, width, height,format.endsWith(".gz"),FunctionalMap.make());
			}
		throw new IllegalArgumentException("Cannot find output format for "+output);
		}
	}
