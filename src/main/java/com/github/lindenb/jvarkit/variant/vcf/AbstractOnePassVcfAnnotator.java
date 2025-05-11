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
package com.github.lindenb.jvarkit.variant.vcf;



import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/** a VCF filtering variant using a set of VariantAnnotator */
public abstract class AbstractOnePassVcfAnnotator extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(AbstractOnePassVcfAnnotator.class);

	private final List<VariantAnnotator> annotators =  new ArrayList<>();
	
	/** return the {@link VariantAnnotator} . will be called in beforeVcf */
	protected abstract List<VariantAnnotator> createVariantAnnotators();
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	/** createVariantAnnotators is called in init */
	@Override
	protected int beforeVcf() {
		try {
			this.annotators.addAll(createVariantAnnotators());
			if(this.annotators.isEmpty()) {
				LOG.warn("NO Variant Annotator was loaded !");
				}
			}
		catch(Throwable err) {
			getLogger().error(err);
			return -1;
			}
		return super.beforeVcf();
		}
	
	private void recursive(VariantContextWriter w, final int annotator_idx, final VariantContext ctx) throws IOException {
		if(annotator_idx==this.annotators.size()) {
			w.add(ctx);
			}
		else
			{
			for(VariantContext ctx2: this.annotators.get(annotator_idx).annotate(ctx)) {
				recursive(w,annotator_idx+1,ctx2);
				}
			}
		}
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator r, final VariantContextWriter w) {
		try {
			final VCFHeader header= r.getHeader();
			final VCFHeader h2=new VCFHeader(header);
			for(VariantAnnotator ann:this.annotators) {
				ann.fillHeader(h2);
				}
			
			JVarkitVersion.getInstance().addMetaData(this, h2);			
			w.writeHeader(h2);
	
			while(r.hasNext())
				{
				recursive(w,0,r.next());
				}
			return 0;
			} catch(final Throwable err ) {
			LOG.error(err);
			return -1;
			} 
		}
	
	@Override
	protected void afterVcf() {
		for(VariantAnnotator ann: this.annotators) {
			ann.close();
			}
		super.afterVcf();
		}
	}
