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
package com.github.lindenb.jvarkit.tools.vcfanonymize;

import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.VCFWriterAnonymizer;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
@Program(name="vcfanonym",
	description="Anonymize samples from VCF",
	keywords={"vcf","annotation"},
	modificationDate="20260724",
	creationDate="20260724",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfAnonymize extends OnePassVcfLauncher {

	public VcfAnonymize() {
		}

	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		final VCFWriterAnonymizer delegate = new VCFWriterAnonymizer(out);
		final VCFHeader h2 = iterin.getHeader();
		JVarkitVersion.getInstance().addMetaData(this, h2);
		delegate.writeHeader(h2);
		VCFUtils.copyVariantsTo(iterin, delegate);
		return 0;
		}
	
	public static void main(final String[] args) {
		new VcfAnonymize().instanceMainWithExit(args);
	}

}
