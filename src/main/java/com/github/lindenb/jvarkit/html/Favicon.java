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
package com.github.lindenb.jvarkit.html;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
/**
 * URL for a favicon
 */
public interface Favicon {
public static Favicon make() {
	return new Favicon() {
		@Override
		public String getURL() {
			return "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAIAAACQd1PeAAAAAXNSR0IB2cksfwAAAARnQU1BAACxjwv8YQUAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAAAlwSFlzAAAuIwAALiMBeKU/dgAAAAxJREFUCNdj8NttCgACmQE/IckioQAAAABJRU5ErkJggg==";
			}
		@Override
		public String toString() {
			return getURL();
			}
		};
	}
public String getURL();
	
/** write link element */
public default void write(final XMLStreamWriter w) throws XMLStreamException {
    w.writeEmptyElement("link");
    w.writeAttribute("rel", "icon");
    w.writeAttribute("href",getURL());
	}

/** append link element to head element 
 * 
 * @return the link element
 * */
public default Element appendToHead(final Element headElement) {
	final Document dom = headElement.getOwnerDocument();
	if(dom==null) return null;
	final Element link = dom.createElement("link");
	headElement.appendChild(link);
	link.setAttribute("rel", "icon");
	link.setAttribute("href", getURL());
	return link;	
	}

}
