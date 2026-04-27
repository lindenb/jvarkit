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
package com.github.lindenb.jvarkit.gtf;

import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.samtools.util.AbstractLocatable;

public class GTFTreeNode extends AbstractLocatable {
	/* package */ GTFTreeNode parent = null;
	private GTFTreeNode firstChild = null;
	private GTFTreeNode nextSibling = null;
	/* package */ GTFLine line = null;
	/* package */ GTFTreeNode() {
		}
	
	public boolean hasParentNode() {
		return this.parent!=null;
		}
	public GTFTreeNode getParent() {
		return parent;
		}
	
	public GTFLine getDelegate() {
		return this.line;
		}
	
	void addChild(GTFTreeNode child) {
		if(child.parent!=null && child.parent!=this) throw new IllegalStateException();
		if(child.line==null) throw new IllegalStateException();
		child.parent=this;
		if(this.firstChild==null) {
			this.firstChild= child;
			}
		else
			{
			GTFTreeNode c = this.firstChild;
			while(c.nextSibling!=null) {
				c=c.nextSibling;
				}
			c.nextSibling=child;
			}
		}
	
	void cleanup() {
		final List<GTFTreeNode> L=getChildren();
		int n1= L.size();
		L.stream().forEach(C->C.cleanup());
		L.removeIf(C->C.line==null);
		int n2  =L.size();
		if(n1!=n2) {
			this.firstChild=null;
			if(!L.isEmpty()) {
				this.firstChild= L.get(0);
				for(int i=0;i+1< L.size();i++) {
					L.get(i).nextSibling = L.get(i+1);
					}
				L.get(L.size()-1).nextSibling=null;
				}
			}
		}
	
	
	@Override
	public String getContig() {
		return getDelegate().getContig();
		}
	@Override
	public int getStart() {
		return getDelegate().getStart();
		}
	@Override
	public int getEnd() {
		return getDelegate().getEnd();
		}
	
	public String getType() {
		return getDelegate().getType();
		}
	public String getAttribute(String key,String def) {
		return getDelegate().hasAttribute(key)?getDelegate().getAttribute(key):def;
		}
	
	
	@Override
	public int hashCode() {
		return line.hashCode();
		}
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof GTFTreeNode)) return false;
		return this.line.equals(GTFTreeNode.class.cast(obj).line);
		}
	
	public List<GTFTreeNode> getChildren() {
		final List<GTFTreeNode> L= new ArrayList<>();
		GTFTreeNode c=this.firstChild;
		while(c!=null) {
			L.add(c);
			c=c.nextSibling;
			}
		return L;
		}
	
	public boolean isGene() { return getDelegate().isGene();}
	public boolean isTranscript() { return getDelegate().isTranscript();}
	public boolean isExon() { return getDelegate().isExon();}
	public boolean isCDS() { return getDelegate().isCDS();}
	
	
}
