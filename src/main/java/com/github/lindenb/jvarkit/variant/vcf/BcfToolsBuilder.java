/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.vcf.readers.DelegateVcfIterator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class BcfToolsBuilder {
	
public interface BcfToolsView {
	public String getFile();
	public VCFHeader getHeader();
	public VCFIterator open();
	public VCFIterator query(final Locatable locatable);
}

private static class HtsJdkView implements BcfToolsView {
	private VCFHeader header = null;
	private final String file;
	HtsJdkView(final String file) {
		this.file= file;
	}
	@Override
	public String getFile() {
		return this.file;
		}
	@Override
	public VCFHeader getHeader() {
		if(this.header==null) {
			try(VCFIterator r=open()) {
				r.getHeader();//do nothing but hide warning from javac
				}
			}
		return this.header;
		}
	@Override
	public VCFIterator open() {
		try {
			return new VCFIteratorBuilder().open(this.getFile());
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	@Override
	public VCFIterator query(final Locatable locatable) {
		final VCFFileReader r = new VCFFileReader(Paths.get(getFile()), true);
		final CloseableIterator<VariantContext> iter = r.query(locatable);
		final PeekableIterator<VariantContext> peeker = new PeekableIterator<>(iter);
		return new VCFIterator() {
			@Override
			public VariantContext next() {
				return peeker.next();
			}
			
			@Override
			public boolean hasNext() {
				return peeker.hasNext();
			}
			
			@Override
			public void close() {
				peeker.close();
				iter.close();
				r.close();
			}
			
			@Override
			public VariantContext peek() {
				return peeker.next();
			}
			
			@Override
			public VCFHeader getHeader() {
				return r.getFileHeader();
			}
		};
		}
	}

private class BcfToolsViewImpl implements BcfToolsView {
	private final String file;
	private VCFHeader header = null;
	BcfToolsViewImpl(final String file) {
		this.file= file;
	}
	
	@Override
	public String getFile() {
		return this.file;
	}
	@Override
	public VCFHeader getHeader() {
		if(this.header==null) {
			try(VCFIterator r=open()) {
				r.getHeader();//do nothing but hide warning from javac
				}
			}
		return this.header;
		}
	@Override
	public VCFIterator open() {
		return _query(null);
		}
	@Override
	public VCFIterator query(final Locatable locatable) {
		return _query(locatable);
		}
	private VCFIterator _query(final Locatable locatable) {
		final List<String> cmd = new ArrayList<>();
		cmd.add( BcfToolsBuilder.this.bcftools_exe);
		cmd.add("view");
		cmd.add(getFile());
		if(locatable!=null) {
			cmd.add(locatable.getContig()+":"+locatable.getStart()+"-"+locatable.getEnd());
			}
		return _create(cmd);
		}
	
	private VCFIterator _create(final List<String> cmd)  {
		Process proc = null;
		try {
			proc =new ProcessBuilder(cmd).
				redirectError(Redirect.INHERIT).
				start();
			}
		catch(final Throwable err) {
			throw new RuntimeIOException(err);
			}
		VCFIterator delegate = null;
		InputStream in = null;
		try {
			in  = proc.getInputStream();
			delegate = new  VCFIteratorBuilder().open(in);
			}
		catch(final Throwable err) {
			CloserUtil.close(in);
			proc.destroy();
			throw new RuntimeIOException(err);
			}
		final BcfToolsVcfIterator biter = new BcfToolsVcfIterator(delegate,proc);
		if(this.header==null) this.header=biter.getHeader();
		return biter;
		}

	}


public BcfToolsView getView(final String file) {
	if(!IOUtil.isUrl(file)) {
		final BCF2Codec codec = new BCF2Codec();
		if(codec.canDecode(file)) {
			return new HtsJdkView(file);
			}
		}
	return new BcfToolsViewImpl(file);
	}
	
private String bcftools_exe;
public BcfToolsBuilder() {
	this.bcftools_exe =  System.getProperty("bcftools.exe");
	
	if(StringUtils.isBlank(this.bcftools_exe)) {
		for(String dir:StringUtils.ifBlank(System.getenv("PATH"),"").split(Pattern.quote(File.pathSeparator))) {
			if(StringUtils.isBlank(dir)) continue;
			final File exe = new File(dir,"bcftools");
			if(!exe.exists()) continue;
			if(exe.isDirectory()) continue;
			if(!exe.canExecute()) continue;
			this.bcftools_exe = exe.getPath();
			break;
			}
		}
	
	if(StringUtils.isBlank(this.bcftools_exe)) {
		this.bcftools_exe = "bcftools";
		}
	
	}




private class BcfToolsVcfIterator extends DelegateVcfIterator  {
	final Process proc;
	public BcfToolsVcfIterator(VCFIterator delegate,final Process proc) {
		super(delegate);
		this.proc = proc;
		}
	@Override
	public void close() {
		super.close();
		this.proc.destroy();
		}
	}

}