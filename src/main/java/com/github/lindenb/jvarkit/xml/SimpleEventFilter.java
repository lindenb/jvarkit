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
package com.github.lindenb.jvarkit.xml;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
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
	private int reject_depth = 0;
	private final Predicate<List<QName>> acceptQName;
	private final LinkedList<QName> stack = new LinkedList<>();
	public SimpleEventFilter(final Predicate<List<QName>> acceptQName) {
		this.acceptQName = Objects.requireNonNull(acceptQName);
		}
	
	private boolean testStack() {
		return acceptQName.test(Collections.unmodifiableList(this.stack));
		}
	
	@Override
	public boolean accept(final XMLEvent event) {
		if(event.isProcessingInstruction()) return false;
		if(event.isStartElement()) {
			final QName qName = event.asStartElement().getName();
			this.stack.add(qName);
			if(!testStack() ) {
				reject_depth++;
				}
			return reject_depth==0;
			}
		else if(event.isEndElement()) {
			final boolean curr  = reject_depth==0;
			if(!testStack())  reject_depth--;
			final QName qName = event.asEndElement().getName();
			final QName pop = stack.removeLast();
			if(!qName.equals(pop)) throw new IllegalStateException("expected "+pop+" but got "+qName);
			return curr;
			}
		return reject_depth==0;
		}
	

}
