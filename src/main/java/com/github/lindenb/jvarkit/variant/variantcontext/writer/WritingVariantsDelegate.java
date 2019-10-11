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
package com.github.lindenb.jvarkit.variant.variantcontext.writer;

import java.nio.file.Path;

import com.beust.jcommander.Parameter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class WritingVariantsDelegate {
private SAMSequenceDictionary dict;
@Parameter(names= {"--generate-vcf-md5"},description="Generate MD5 checksum for VCF output.")
private boolean generate_md5 = false;

public WritingVariantsDelegate dictionary(final VCFHeader header) {
	return dictionary(header==null?null:header.getSequenceDictionary());
	}

public WritingVariantsDelegate dictionary(final SAMSequenceDictionary dict) {
	this.dict = dict;
	return this;
	}
	
public VariantContextWriter open(final Path pathOrNull) {
	final VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
	
	
	if(this.dict!=null) vcwb.setReferenceDictionary(this.dict);
	vcwb.clearOptions();
	
	if(pathOrNull!=null) {
		vcwb.setCreateMD5(this.generate_md5);
		vcwb.setOutputPath(pathOrNull);
		}
	else
		{
		vcwb.setCreateMD5(false);
		vcwb.setOutputStream(System.out);
		}
	return vcwb.build();
	}
}
