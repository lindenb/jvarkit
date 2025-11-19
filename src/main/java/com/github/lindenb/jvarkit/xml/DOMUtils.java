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

import org.w3c.dom.Attr;
import org.w3c.dom.CharacterData;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.lang.StringUtils;

public class DOMUtils {

	
	/** returns human XPATH-like path for this node 
	 * @param node the node
	 * @return the path
	 */
	public static String getNodePath(final org.w3c.dom.Node node) {
		if(node==null) return "null";
		String s;
		switch(node.getNodeType()) {
			case Node.CDATA_SECTION_NODE: s= "#cdata"; break;
			case Node.COMMENT_NODE: s= "#comment"; break;
			case Node.TEXT_NODE: s= "#text"; break;
			case Node.DOCUMENT_NODE: s= "<doc>"; break;
			case Node.ATTRIBUTE_NODE: 
				s= "@"+Attr.class.cast(node).getNodeName();
				Element E=Attr.class.cast(node).getOwnerElement();
				if(E!=null) s = getNodePath(E)+"/"+s;
				break;
			case Node.ELEMENT_NODE: 
				String suffix="";
				if(node.getParentNode()!=null) {
					int i=0;
					for(Node c=node.getParentNode().getFirstChild();c!=null;c=c.getNextSibling()) {
						if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
						if(Element.class.cast(node).getNodeName().equals(Element.class.cast(c).getNodeName())) {
							i++;
							if(c==node) {
								suffix="["+i+"]";
								break;
								}
							}
						}
					}
				s= Element.class.cast(node).getNodeName()+suffix;
				break;
			default: s="TODO";break;
			}
		return (node.getParentNode()!=null?getNodePath(node.getParentNode())+"/":"")+s;
		}

	/** retturn true if node is null or element without children (or blank text), or blank text */
	public static boolean isBlank(final org.w3c.dom.Node node) {
		if(node==null) return true;
		switch(node.getNodeType()) {
			case Node.COMMENT_NODE: return true;
			case Node.TEXT_NODE:
			case Node.CDATA_SECTION_NODE:
				return StringUtils.isBlank(CharacterData.class.cast(node).getData());
			case Node.DOCUMENT_NODE: //cont
			case Node.ELEMENT_NODE: //cont
			case Node.DOCUMENT_FRAGMENT_NODE:
				for(Node c=node.getFirstChild();c!=null;c=c.getNextSibling()) {
					if(c.getNodeType()==Node.ELEMENT_NODE) return false;
					if(!isBlank(c)) return false;
					}
				return true;
			default:throw new IllegalArgumentException("TODO");
			}
		}


}
