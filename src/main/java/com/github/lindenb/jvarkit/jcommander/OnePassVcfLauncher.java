/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.jcommander;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

public abstract class OnePassVcfLauncher extends Launcher {
private static final Logger LOG = Logger.build(OnePassVcfLauncher.class).make();
@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
protected Path outputFile=null;
@ParametersDelegate
protected WritingVariantsDelegate writingVariantsDelegate= new WritingVariantsDelegate();

private static class VCFIter implements VCFIterator {
	final VCFIterator delegate;
	final ProgressFactory.Watcher<VariantContext> progess;
	VCFIter(final VCFIterator delegate,final Logger logger) {
		this.delegate = delegate;
		this.progess =ProgressFactory.newInstance().
				dictionary(delegate.getHeader()).
				logger(logger).
				build();
		}
	@Override
	public VCFHeader getHeader() {
		return this.delegate.getHeader();
		}
	@Override
	public boolean hasNext() {
		return this.delegate.hasNext();
		}
	@Override
	public VariantContext peek() {
		return this.delegate.peek();
		}
	@Override
	public VariantContext next() {
		return progess.apply(this.delegate.next());
		}
	@Override
	public void close() {
		this.progess.close();
		this.delegate.close();
		}
	}

protected Logger getLogger() {
	return null;
	}

private void deleteOutputOnError() {
	if(this.outputFile==null) return;
	try {
		if(Files.deleteIfExists(this.outputFile)) {
			LOG.warning("The following file was deleted: "+this.outputFile);
			}
	} catch(final IOException err) {
		//ignore
	}
}
/** initialize things before opening the vcf */
protected int beforeVcf() {
	return 0;
	}

/** initialize things after closing the vcf */
protected void afterVcf() {
	}

@Override
public int doWork(final List<String> args) {
	VCFIterator in = null;
	VariantContextWriter vcw = null;
	final String input = super.oneFileOrNull(args);
	if(input!=null &&
		this.outputFile!=null &&
		!IOUtil.isUrl(input) && Paths.get(input).equals(this.outputFile)) {
		LOG.error("Input == output : "+ input);
		return -1;
		}
	
	try {
		if(beforeVcf()!=0) {
			LOG.error("initialization failed");
			return -1;
			}
		}
	catch (final Throwable err) {
		LOG.error(err);
		return -1;
		}
	
	try {
		final BcfIteratorBuilder bcb = new BcfIteratorBuilder();
		if(input==null) {
			in = bcb.open(stdin());
			}
		else {
			in = bcb.open(input);
			}
		if(getLogger()!=null) {
			in = new VCFIter(in, getLogger());
			}
		vcw = this.writingVariantsDelegate.dictionary(in.getHeader()).open(this.outputFile); 
		final int err = doVcfToVcf(input==null?"<stdin>":input, in,vcw);
		vcw.close();
		vcw=null;
		in.close();
		in = null;
		if(err!=0) deleteOutputOnError();
		return err;
		}
	catch (final Throwable err) {
		LOG.error(err);
		if(vcw!=null) vcw.close();
		deleteOutputOnError();
		return -1;
		}
	finally
		{
		try {
			afterVcf();
			}
		catch (final Throwable err) {
			LOG.error(err);
			}	
		if(in!=null) in.close();
		if(vcw!=null) vcw.close();
		}
	}
}
