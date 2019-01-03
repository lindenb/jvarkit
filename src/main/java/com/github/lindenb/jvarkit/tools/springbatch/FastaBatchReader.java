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
package com.github.lindenb.jvarkit.tools.springbatch;

import java.io.InputStreamReader;
import java.util.Objects;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.item.file.ResourceAwareItemReaderItemStream;
import org.springframework.batch.item.support.AbstractItemCountingItemStreamItemReader;
import org.springframework.core.io.Resource;

import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequence;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequenceReader;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

public class FastaBatchReader
			extends AbstractItemCountingItemStreamItemReader<FastaSequence> implements
			ResourceAwareItemReaderItemStream<FastaSequence>
	{
	private static final Log LOG = LogFactory.getLog(FastaBatchReader.class);

	private CloseableIterator<FastaSequence> fastaIterator = null;
	private Resource resource = null;
	private int sequenceCapacity = 1000;
	@Override
	public void setResource(final Resource rsrc) {
		this.resource = rsrc;
		}
	public Resource getResource() {
		return resource;
		}
	
	public void setSequenceCapacity(int sequenceCapacity) {
		this.sequenceCapacity = sequenceCapacity;
		}
	
	public int getSequenceCapacity() {
		return sequenceCapacity;
		}
	
	@Override
	protected void doOpen() throws Exception {
		Objects.requireNonNull(getResource(), "Input resource must be set");
		if(LOG.isInfoEnabled()) LOG.info("opening "+getResource());
		if(!this.getResource().isReadable()) throw new IllegalStateException("Resource is not readable.");
		final FastaSequenceReader fastaReader = new FastaSequenceReader();
		fastaReader.setSequenceCapacity(getSequenceCapacity());
		this.fastaIterator = fastaReader.iterator(
			new InputStreamReader(this.getResource().getInputStream()));
		}
	@Override
	protected void doClose() throws Exception {
		CloserUtil.close(this.fastaIterator);
		this.fastaIterator = null;
		}

	@Override
	protected FastaSequence doRead() throws Exception {
		return  this.fastaIterator==null || !this.fastaIterator.hasNext() ?
				null :
				this.fastaIterator.next()
				;
		}
	
}
