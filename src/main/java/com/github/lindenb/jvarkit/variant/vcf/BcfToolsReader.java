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
package com.github.lindenb.jvarkit.variant.vcf;

import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

public class BcfToolsReader implements VCFReader {
	private static final Logger LOG = Logger.build(BcfToolsReader.class).make();

	private static int WARNING = 0;
	private final VCFCodec codec = new VCFCodec();
	private final VCFHeader header;
	private final String path;
	private List<String> extraParams;
	public BcfToolsReader(final String path) {
		this(path,null);
		}
	public BcfToolsReader(final String path,final List<String> extraParams) {
		this.path = path;
		if(extraParams==null || extraParams.isEmpty()) {
			this.extraParams= Collections.emptyList();
		} else {
			this.extraParams = new ArrayList<>(this.extraParams);
		}
		
		Process proc = null;
		if(WARNING==0) {
			WARNING=1;
			LOG.warn("the htsjdk library is current not able to read bcf file version > " + BCF2Codec.ALLOWED_BCF_VERSION + 
					". This tool will try to read bcf file using bcftools.");
			}
		final String cmd[] = {
				BcfToolsUtils.getBcftoolsExecutable() ,
				"view",
				"--header-only",
				path
				};
		try {
			proc =new ProcessBuilder(cmd).
				redirectError(Redirect.INHERIT).
				start();
			try(htsjdk.tribble.readers.LineReader  in = AsciiLineReader.from(new PositionalBufferedStream(proc.getInputStream()))) {
				final LineIterator r= new LineIteratorImpl(in);
				this.header = (VCFHeader)codec.readActualHeader(r);
				}
			}
		catch(final Throwable err) {
			proc.destroy();
			throw new RuntimeIOException(err);
			}
		finally
			{
			}
		}
	@Override
	public VCFHeader getHeader() {
		return header;
		}
	
	public boolean isQueryable () {
		return Files.exists(Paths.get(getPath()+FileExtensions.CSI)) || Files.exists(Paths.get(getPath()+FileExtensions.TABIX_INDEX));
		}
	@Override
	public CloseableIterator<VariantContext> iterator() {
		return create(null);
	}
	
	public CloseableIterator<VariantContext> query(final Locatable loc) {
		return query(loc.getContig(),loc.getStart(),loc.getEnd());
	}
	@Override
	public CloseableIterator<VariantContext> query(String s,int start,int end) {
		return create(s+":"+start+"-"+end);
	}
	
	public String getPath() {
		return path;
		}
	
	private CloseableIterator<VariantContext> create(final String location) {
		final List<String> cmd =  new ArrayList<>();
		cmd.add(BcfToolsUtils.getBcftoolsExecutable());
		cmd.add("view");
		cmd.add("--no-header");
		cmd.addAll(this.extraParams);
		cmd.add(getPath());
		if(location!=null) cmd.add(location);
		try {
			return new BcfIterator(this.codec,new ProcessBuilder(cmd).
					redirectError(Redirect.INHERIT).
					start());
		} catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	@Override
	public void close() throws IOException {
		//nothing
		}
	
	private class BcfIterator extends AbstractCloseableIterator<VariantContext> {
		final Process proc;
		final VCFCodec codec;
		final htsjdk.tribble.readers.LineReader r;
		BcfIterator(final VCFCodec codec,final Process proc) {
			this.proc = proc;
			this.codec = codec;
			this.r = AsciiLineReader.from(new PositionalBufferedStream(proc.getInputStream()));
			}
		@Override
		protected VariantContext advance() {
			try {
				final String line = r.readLine();
				if(line==null) {
					this.proc.waitFor();
					return null;
					}
				return this.codec.decode(line);
				}
			catch(final Throwable err) {
				throw new RuntimeIOException(err);
				}
			
			}
		@Override
		public void close() {
			r.close();
			proc.destroy();
			}
		}
	

	}
