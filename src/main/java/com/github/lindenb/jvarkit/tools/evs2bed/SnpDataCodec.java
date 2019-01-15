/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.evs2bed;

import java.io.StringReader;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.stream.StreamSource;
import edu.washington.gs.evs.SnpData;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;

public class SnpDataCodec
extends AsciiFeatureCodec<SnpDataFeature>
	{
	public static final String XML_FILE_SUFFIX=".xml";
	@SuppressWarnings("unused")
	private static edu.washington.gs.evs.ObjectFactory _fooljavac=null;
	private Unmarshaller unmarshaller;
	public static final String ROOT_ELEMENT="evsData";
	
	public SnpDataCodec()
		{
		super(SnpDataFeature.class);
		try {
			JAXBContext jc = JAXBContext.newInstance(SnpData.class);
	        this.unmarshaller = jc.createUnmarshaller();			}
		catch (Exception e)
			{
			throw new RuntimeException(e);
			}
		}
	
	@Override
	public boolean canDecode(String path)
		{
		return path.endsWith(XML_FILE_SUFFIX);	
		}
	@Override
	public Object readActualHeader(LineIterator reader) {
		StringBuilder sb=new StringBuilder();
		String line;
		while((line=reader.peek())!=null && !line.startsWith("<snpList"))
			{
			sb.append(reader.next());
			}
		return sb.toString();
		}
	
	
	@Override
	public SnpDataFeature decode(String line)
		{
		try {
			//last line of xml file
			if(line.startsWith("</")) return null;
			//System.err.println(line);
			JAXBElement<SnpData> jaxbElement = (JAXBElement<SnpData>)this.unmarshaller.unmarshal(
					new StreamSource(new StringReader(line)),
					SnpData.class
					);
			return new SnpDataFeature(jaxbElement.getValue());
		} catch (Exception e) {
			throw new RuntimeException(e);
			}
		}
	
	}
