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

*/
package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.function.Supplier;

import htsjdk.variant.vcf.VCFHeader;

public abstract class AbstractPredictionParserFactory<T extends PredictionParser,FACTORYTYPE> 
	implements Supplier<T>
	{
	private String _tag=null;
	private VCFHeader _header=null;
	protected AbstractPredictionParserFactory(final String tag) {
		this._tag=tag;
		}
	
	protected AbstractPredictionParserFactory(final String tag,final VCFHeader header) {
		this(tag);
		this._header=header;
		}
	
	public void setTag(final String tag) {
		this._tag=tag;
		}
	public String getTag() {
		return this._tag;
		}
	public void setHeader(final VCFHeader header) {
		this._header=header;
		}
	public VCFHeader getHeader() {
		return this._header;
		}
	@SuppressWarnings("unchecked")
	public FACTORYTYPE header(final VCFHeader vcfHeader) {
		this.setHeader(vcfHeader);
		return (FACTORYTYPE)this;
	}
		
	@SuppressWarnings("unchecked")
	public FACTORYTYPE tag(final String tag) {
		this.setTag(tag);
		return (FACTORYTYPE)this;
		}
	}
