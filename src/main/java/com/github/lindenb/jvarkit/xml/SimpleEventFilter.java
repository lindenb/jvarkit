/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.xml;

import java.util.Objects;
import java.util.function.Predicate;

import javax.xml.namespace.QName;
import javax.xml.stream.EventFilter;
import javax.xml.stream.events.XMLEvent;



/**
Example:

```
final XMLEventReader reader0 = inputFactory.createXMLEventReader(r);
final XMLEventReader reader = inputFactory.createFilteredReader(reader0, new OWLFilter());
StAX2Model.read(reader,this.ontModel,"file://"+path.toString()); 
```
*/
public class SimpleEventFilter implements EventFilter {
	private int in_reject_depth = 0;
	private final Predicate<QName> acceptQName;
	public SimpleEventFilter(final Predicate<QName> acceptQName) {
		this.acceptQName = Objects.requireNonNull(acceptQName);
		}
	
	private boolean acceptName(final QName qName) {
		return acceptQName.test(qName);
		}
	
	@Override
	public boolean accept(final XMLEvent event) {
		if(event.isProcessingInstruction()) return false;
		if(event.isStartElement()) {
			if(!acceptName( event.asStartElement().getName()) ) {
				in_reject_depth++;
				}
			//System.err.println("<"+event.asStartElement().getName().getLocalPart()+"> Returning "+(in_reject_depth==0)+" and after will be "+in_reject_depth);
			return in_reject_depth==0;
			}
		else if(event.isEndElement()) {
			boolean curr  = in_reject_depth==0;
			if( !acceptName(event.asEndElement().getName())) in_reject_depth--;
			//System.err.println("</"+event.asEndElement().getName().getLocalPart()+"> Returning "+curr+" and after will be "+in_reject_depth);
			return curr;
			}
		//System.err.println("("+event.getEventType()+") in_reject_depth "+in_reject_depth);
		return in_reject_depth==0;
		}
	

}
