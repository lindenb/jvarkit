/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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

* 2016 creation

*/
package com.github.lindenb.jvarkit.tools.vcfeigen;


import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;



/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */
public class VcfEigen01
	extends AbstractVcfEigen01
	{	
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfEigen01.class);
	
	public VcfEigen01()
		{
		
		}
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator r,
			final VariantContextWriter w
			) throws IOException
		{
		EigenInfoAnnotator annotator = null;
		try
			{
			if(eigenDirStr==null || eigenDirStr.trim().isEmpty()) {
				throw new IOException("Eigen directory is undefined : option -"+OPTION_EIGENDIRSTR);
			}
			final File eigenDirectory = new File(super.eigenDirStr);
			LOG.info("loading eigen directory "+eigenDirectory);
			annotator = new EigenInfoAnnotator(eigenDirectory);
			final VCFHeader header = r.getHeader();
			
			final  VCFHeader h2 = new VCFHeader(header);
			addMetaData(h2);
			for(final VCFInfoHeaderLine vihl: annotator.getInfoHeaderLines()) {
				if(h2.getInfoHeaderLine(vihl.getID())!=null) {
					throw new IOException("VCF INFO "+vihl.getID()+" already defined in input VCF.");
				}
				h2.addMetaDataLine(vihl);
			}
			
			
			final  SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
			w.writeHeader(h2);
			while (r.hasNext())
				{
				final  VariantContext variation = progress.watch(r.next());
				final Map<String,Object> m  = annotator.getAnnotations(variation);
				if(m==null || m.isEmpty())
					{
					w.add(variation);
					}
				else
					{
					final VariantContextBuilder vcb = new VariantContextBuilder(variation);
					for(final String key: m.keySet()) {
						vcb.attribute(key, m.get(key));
					}
					w.add(vcb.make());
					}
				}
			progress.finish();
			return RETURN_OK;
			}
		finally
			{
			CloserUtil.close(annotator);
			}
		}
	
	@Override
	public Collection<Throwable> initializeKnime()
		{
		if(super.eigenDirStr==null || super.eigenDirStr.trim().isEmpty())
			{
			return wrapException("Eigen directory is undefined : option -"+OPTION_EIGENDIRSTR);
			}
		return super.initializeKnime();
		}
	
	
	
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	
	public static void main(String[] args) throws Exception
		{
		new VcfEigen01().instanceMainWithExit(args);
		}

	}
