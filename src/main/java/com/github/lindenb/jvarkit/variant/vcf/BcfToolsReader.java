package com.github.lindenb.jvarkit.variant.vcf;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.LineReader;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class BcfToolsReader implements Closeable {
	private final VCFCodec codec = new VCFCodec();
	private final VCFHeader header;
	private final String path;
	public BcfToolsReader(final String path) {
		this.path = path;
		Process proc = null;
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
			try(htsjdk.tribble.readers.LineReader  in = AsciiLineReader.from(proc.getInputStream())) {
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
	
	public VCFHeader getHeader() {
		return header;
		}
	
	public boolean isQueryAble() {
		return Files.exists(Paths.get(getPath()+FileExtensions.CSI)) || Files.exists(Paths.get(getPath()+FileExtensions.TABIX_INDEX));
	}
	
	public CloseableIterator<VariantContext> iterator() {
		return create(null);
	}
	
	public CloseableIterator<VariantContext> query(Locatable loc) {
		return query(loc.getContig(),loc.getStart(),loc.getEnd());
	}
	
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
		cmd.add(getPath());
		if(location!=null) cmd.add(location);
		try {
		return new BcfIterator(this.codec,new ProcessBuilder(cmd).
				redirectError(Redirect.INHERIT).
				start());
		} catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	@Override
	public void close() throws IOException {
		//nothing
		}
	
	private class BcfIterator extends AbstractIterator<VariantContext>
		implements CloseableIterator<VariantContext> {
		final Process proc;
		final VCFCodec codec;
		final htsjdk.tribble.readers.LineReader r;
		BcfIterator(final VCFCodec codec,final Process proc) {
			this.proc = proc;
			this.codec = codec;
			r = AsciiLineReader.from(proc.getInputStream());
			
			}
		@Override
		protected VariantContext advance() {
			try {
				final String line = r.readLine();
				if(line==null) {
					this.proc.waitFor();
					return null;
					}
				return codec.decode(line);
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
