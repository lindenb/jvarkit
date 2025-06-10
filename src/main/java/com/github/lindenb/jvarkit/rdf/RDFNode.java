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
package com.github.lindenb.jvarkit.rdf;

public interface RDFNode extends Comparable<RDFNode> {
	public boolean isResource();
	public boolean isLiteral();
	public default Resource asResource() {
		return Resource.class.cast(this);
		}
	public default Literal asLiteral() {
		return Literal.class.cast(this);
		}
	@Override
	default int compareTo(RDFNode o) {
		if(this.isResource()) {
			if(o.isLiteral()) {
				return -1;
				}
			else if(o.isResource()) {
				return asResource()._compareTo(o.asResource());
				}
			else
				{
				throw new IllegalStateException();
				}
			}
		else if(this.isLiteral()) {
			if(o.isLiteral()) {
				return asLiteral()._compareTo(o.asLiteral());
				}
			else if(o.isResource()) {
				return 1;
				}
			else
				{
				throw new IllegalStateException();
				}
			}
		else
			{
			throw new IllegalStateException();
			}
		}
	}
