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


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

public abstract class AbstractVariantProcessor
extends RodWalker<Long, Long> implements TreeReducible<Long>
{
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter writer = null;

    
    protected AbstractVariantProcessor(){
    }
    
    @Override
    public void initialize() {
    	super.initialize();
    	}
    
    /** get the VCF header for 'variants' */
    protected VCFHeader getVcfHeader() {
    	if(this.variantCollection==null || this.variantCollection.variants==null) throw new IllegalStateException("variants was not set");
    	final Map<String,VCFHeader> map =GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Collections.singleton(this.variantCollection.variants.getName()));
    	final VCFHeader header= map.get(this.variantCollection.variants.getName());
    	if( header == null) throw new UserException.MalformedVCF("Cannot get VCF header from "+this.variantCollection.variants);
    	return header;
    	}
    
    /** this method is called my mapToMany */
    protected VariantContext mapVariant(final VariantContext ctx,
    		final RefMetaDataTracker tracker,
    		final ReferenceContext ref,
    		final AlignmentContext context)
    	{
    	logger.warn("method map not implemented.");
    	return null;
    	}
    
    /** this method is called my map(tracker,ref,context) */
    protected List<VariantContext> mapVariantToMany(
    		final VariantContext ctx,
    		final RefMetaDataTracker tracker,
    		final ReferenceContext ref,
    		final AlignmentContext context)
    	{
    	final VariantContext ctx2 = mapVariant(ctx,tracker,ref,context);
    	return ctx2==null?Collections.emptyList():Collections.singletonList(ctx2);
    	}
    
    
    
    @Override
    public Long map(
    		final RefMetaDataTracker tracker,
    		final ReferenceContext ref,
    		final AlignmentContext context) {
        if ( tracker == null ) return 0L;
        long processed=0L;
		for(final VariantContext ctx: tracker.getValues(
				this.variantCollection.variants,
				context.getLocation()
				))
			{
			for(final VariantContext ctx2:mapVariantToMany(ctx,tracker,ref,context))
				{
				this.writer.add(ctx2);
				
				}
			processed++;
			}
		return processed;
	    }
    @Override
    public Long reduceInit() { return 0L; }

    @Override
    public Long reduce(final Long value,final Long sum) { return value + sum; }

    @Override
    public Long treeReduce(final Long lhs,final Long rhs) {
        return lhs + rhs;
    }

    @Override
    public void onTraversalDone(final Long result) {
        logger.info("Processed " + result + " variants.\n");
    }

}
