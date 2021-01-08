/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.jcommander.converter;
import java.awt.Dimension;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Iterator;
import java.util.function.Function;
import java.util.zip.GZIPInputStream;

import javax.imageio.ImageIO;
import javax.imageio.ImageReader;
import javax.imageio.stream.ImageInputStream;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;


public class DimensionConverter
	implements Function<String,Dimension>
	{
	public static final String OPT_DESC="a dimension can be specified as '[integer]x[integer]' or it can be the path to an existing png,jpg,xcf,svg file. ";
	
	public static class StringConverter implements IStringConverter<Dimension> {
		@Override
		public Dimension convert(final String dimStr) {
			return new DimensionConverter().apply(dimStr);
			}
	}
	
	
	@Override
	public Dimension apply(final String dimStr) {
		if(StringUtils.isBlank(dimStr)) throw new IllegalArgumentException("empty string");
		if(dimStr.toLowerCase().matches("\\d+x\\d+")) {
			final int x_symbol = dimStr.toLowerCase().indexOf("x");
			return new Dimension(
					Integer.parseInt(dimStr.substring(0,x_symbol)),
					Integer.parseInt(dimStr.substring(1+x_symbol))
					);
			}
	
		final Path f= Paths.get(dimStr);
		if(!Files.exists(f) || !Files.isRegularFile(f)) {
			throw new IllegalArgumentException("not an existing file: "+f);
			}
		if(f.getFileName().toString().endsWith(".svg") || f.getFileName().toString().endsWith(".svg.gz")) {
			try(InputStream fis=Files.newInputStream(f))
				{
				InputStream in=fis;
				if( f.getFileName().toString().endsWith(".gz")) {
					in=new GZIPInputStream(fis);
					}
				Integer w = null;
				Integer h = null;
				XMLEventReader r = XMLInputFactory.newInstance().createXMLEventReader(in);
				while(r.hasNext()) {
					final XMLEvent evt = r.nextEvent();
					if(!evt.isStartElement()) continue;
					final StartElement E = evt.asStartElement();
					if("svg".equals(E.getName().getLocalPart())) {
						Attribute att = E.getAttributeByName(new QName("width"));
						if(att!=null) {
							w= (int)Double.parseDouble(att.getValue());
							}
						 att = E.getAttributeByName(new QName("height"));
						if(att!=null) {
							h= (int)Double.parseDouble(att.getValue());
							}
						}
					break;
					}
				r.close();
				if(w!=null && h!=null) return new Dimension(w,h);
				}
			catch(final Exception err) {
				throw new IllegalArgumentException(err);
				}
			throw new IllegalArgumentException("cannot get dimension from SVG file "+f);
			}
		
		
		if(f.getFileName().toString().endsWith(".xcf"))
				{
				try(InputStream fis=Files.newInputStream(f))
					{
					byte array[]=new byte[9];
					fis.read(array);
					if(!Arrays.equals(array, "gimp xcf ".getBytes()))
						{
						throw new IOException("bad gimp xcf header");
						}
					
					array=new byte[5];
					if(fis.read(array)!=array.length) {
						throw new IOException("bad gimp xcf header");
						}
					//LOG.info("version "+new String(array));
					array=new byte[8];
					if(fis.read(array)!=array.length) {
						throw new IOException("bad gimp xcf header");
						}
					final ByteBuffer buf = ByteBuffer.wrap(array); // big endian by default
				    buf.put(array);
				    buf.position(0);
				    final int w= buf.getInt();
				    final int h= buf.getInt();
				    //LOG.info("width of "+f+" is "+w+"x"+h);
				    return new Dimension(w,h);
					}
				catch(final IOException err) {
					throw new IllegalArgumentException(err);
					}
				}
			
			try(ImageInputStream in = ImageIO.createImageInputStream(f.toFile())){
				if(in!=null) {
				    final Iterator<ImageReader> readers = ImageIO.getImageReaders(in);
				    if (readers.hasNext()) {
				        final ImageReader reader = readers.next();
				        try {
				            reader.setInput(in);
				           return new Dimension(reader.getWidth(0), reader.getHeight(0));
				        } finally {
				            reader.dispose();
					        }
					    }
					}
				} 			
			catch(final IOException err) {
				throw new IllegalArgumentException(err);
				}
			throw new IllegalArgumentException("cannot get dimension from "+dimStr);
			}
		}
