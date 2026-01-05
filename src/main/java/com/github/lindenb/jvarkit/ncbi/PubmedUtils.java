/*

The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.ncbi;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Map;
import java.util.Optional;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.util.Maps;

public class PubmedUtils extends NcbiConstants {
	private static final Map<String,String> MONTH2MONTH = Maps.of(
			"Jan","01","Feb","02","Mar","03","Apr","04","May","05","Jun","06",
			"Jul","07","Aug","08","Sep","09","Oct","10","Nov","11","Dec","12"
			);
	/** returns URL for given PMID */
	public static String pmidToURL(final String pmid) {
		return "https://pubmed.ncbi.nlm.nih.gov/" + pmid;
		}
	
	private static void skip(final XMLEventReader r) throws XMLStreamException
		{
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isEndElement()) break;
			else if(evt.isStartElement())
				{
				skip(r);
				}
			}
		}
	
	public static Optional<Date> scanDate(final XMLEventReader r) throws XMLStreamException {
		String year = null;
		String month="01";
		String day="01";
		
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isStartElement())
				{
				String name =evt.asStartElement().getName().getLocalPart();
				if(name.equals("Year")) {
					year = r.getElementText();
					}
				else if(name.equals("Month")) {
					month =   r.getElementText();
					month = MONTH2MONTH.getOrDefault(month, month);
					}
				else if(name.equals("Day")) {
					day = r.getElementText();
					}
				else
					{
					skip(r);
					}
				}
			else if(evt.isEndElement()) {
				if(year!=null) {
					final String pattern = "yyyy-MM-dd";
					final SimpleDateFormat simpleDateFormat = new SimpleDateFormat(pattern);
					
					try {
						final Date date = simpleDateFormat.parse(String.join("-", year,month,day));
						return Optional.of(date);
						}
					catch(Throwable err) {
						
						}
					return Optional.empty();
					}
				return null;
				}
			}
		
		 return null;
		}
}
