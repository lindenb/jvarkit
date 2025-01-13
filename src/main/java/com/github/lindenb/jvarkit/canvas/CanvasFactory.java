/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.canvas;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Optional;

import javax.xml.stream.XMLStreamException;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.FunctionalMap;

public class CanvasFactory {
	public enum Format {
		SVG,SVG_GZ,PNG,JPG,PS,PS_GZ;
		/** return file suffix, WITHOUT the first dot */
		public String getSuffix() {
			return this.name().toLowerCase().replace('_', '.');
			}
		};
	public static final String OPT_FORMAT_DESC="output format";
	private int width=700;
	private int height=400;
	private Format format=Format.SVG;
	public CanvasFactory() {
		}
	
	public CanvasFactory setFormat(Format format) {
		this.format = format;
		return this;
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
	
	private Optional<Format> formatFromSuffix( String fname) {
		if(StringUtils.isBlank(fname)) return Optional.empty();
		fname=fname.toLowerCase();
		if(fname.endsWith(".png")) return Optional.of(Format.PNG);
		if(fname.endsWith(".jpg")) return Optional.of(Format.JPG);
		if(fname.endsWith(".jpeg")) return  Optional.of(Format.JPG);
		if(fname.endsWith(".svg")) return Optional.of(Format.SVG);
		if(fname.endsWith(".svg.gz")) return Optional.of(Format.SVG_GZ);
		if(fname.endsWith(".ps") || fname.endsWith(".eps")) return Optional.of(Format.PS);
		if(fname.endsWith(".ps.gz") || fname.endsWith(".eps.gz")) return Optional.of(Format.PS_GZ);
		return Optional.empty();
		}
	
	public Canvas open(final Path outputOrNull,final FunctionalMap<String,Object> fmap) throws IOException {
		Format fmt = outputOrNull==null?null:formatFromSuffix(outputOrNull.getFileName().toString()).orElse(null);
		if(fmt==null) fmt=this.format;
		if(fmt==null) throw new IllegalArgumentException("cannot find output format");
		switch(fmt) {
			case PNG: case JPG: {
				return new Graphics2DCanvas(outputOrNull, width, height,fmt,fmap);
				}
			case SVG:
			case SVG_GZ:
				try {
					return new SVGCanvas(outputOrNull, width, height,fmt.equals(Format.SVG_GZ),FunctionalMap.make());
					}
				catch(XMLStreamException err) {
					throw new IOException(err);
					}
			case PS:
			case PS_GZ:
				return new PSCanvas(outputOrNull, width, height,fmt.equals(Format.PS_GZ),FunctionalMap.make());
			default:
				throw new IllegalArgumentException("Cannot find output format for "+outputOrNull+ " / "+fmt);
			}
		}
	}
