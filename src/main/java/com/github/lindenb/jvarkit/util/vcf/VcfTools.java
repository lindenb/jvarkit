/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2017 creation

*/
package com.github.lindenb.jvarkit.util.vcf;

import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.variant.vcf.VCFHeader;

/** bean to be injected in various contexts like javascript context, etc... */
public class VcfTools {
private VCFHeader header=null;
private SnpEffPredictionParser snpEffPredictionParser=null;
private VepPredictionParser vepPredictionParser=null;
private AnnPredictionParser annPredictionParser=null;

public VcfTools() {
	
	}
public VcfTools(final VCFHeader header) {
	init(header);
	}
public void init(final VCFHeader header) {
	this.header=header;
	if(header!=null)
		{
		this.snpEffPredictionParser=new SnpEffPredictionParserFactory().header(header).get();
		this.vepPredictionParser=new VepPredictionParserFactory().header(header).get();
		this.annPredictionParser=new AnnPredictionParserFactory().header(header).get();
		}
	else
		{
		this.snpEffPredictionParser = null;
		this.vepPredictionParser= null;
		this.annPredictionParser= null;
		}
	
	}
}
