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
