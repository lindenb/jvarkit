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
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.eclipse.jetty.io.RuntimeIOException;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.MayBeGzipInputStream;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.bcf2.BCFVersion;
import htsjdk.variant.vcf.VCFFileReader;

/** Factory creating VCFReader */
public abstract class VCFReaderFactory {
	private boolean requireIndex = true; /* because default is true in VCFFileReader */
	
	/** initialize with requireIndex = true */
	protected VCFReaderFactory() {	
		}
	
	public static VCFReaderFactory makeDefault() {
		return new VCFReaderFactory() {
			};
		}
	
	public VCFReaderFactory setRequireIndex(boolean b) {
		this.requireIndex = b;
		return this;
		}
	
	public boolean isRequireIndex() {
		return requireIndex;
		}
	
	private void tabixedURL(String url,boolean requireIndex) throws IOException
		{
		if(requireIndex) {
			final String tbiUrl= url + FileExtensions.TABIX_INDEX;
			try(InputStream is=new URL(tbiUrl).openStream()) {
				Path tmpTbi = Files.createTempFile(IOUtils.getDefaultTempDir(),"tabix", FileExtensions.TABIX_INDEX);
				IOUtils.copyTo(is, tmpTbi);
				}
			}
		} 
	
	/** open new VCFReader with default {@link #isRequireIndex()} */
	public VCFFileReader open(final String pathOrUrl) {
		return open(pathOrUrl,isRequireIndex());
	}
	
	/** open new VCFReader */
	public VCFFileReader open(final String pathOrUrl,boolean requireIndex) {
		if(IOUtil.isUrl(pathOrUrl)) {
			
		}

		return open(Paths.get(pathOrUrl), requireIndex);
		}

	
	/** open new VCFReader with default {@link #isRequireIndex()} */
	public VCFFileReader open(final Path p) {
		return open(p,isRequireIndex());
		}
	
	/** open new VCFReader */
	public VCFFileReader open(final Path path,boolean requireIndex) {
		if(path.getFileName().toString().endsWith(FileExtensions.BCF)) {
			final BCFVersion version;
			try(InputStream is = new MayBeGzipInputStream(path)) {
				version = BCFVersion.readBCFVersion(is);
				if(version==null) throw new IllegalArgumentException("File "+path+" doesn't look like a BCF file.");
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			if(!BCF2Codec.ALLOWED_BCF_VERSION.equals(version)) {
				throw new IllegalArgumentException("path is bcf version "+version +" but support is for "+BCF2Codec.ALLOWED_BCF_VERSION);
				}
			}
		
		return new VCFFileReader(path, requireIndex);
		}
	
	/** open new VCFReader */
	public final VCFFileReader open(final File path,boolean requireIndex) {
		return open(path.toPath(),requireIndex);
		}
	
	/** open new VCFReader with default {@link #isRequireIndex()} */
	public final VCFFileReader open(final File p) {
		return open(p,isRequireIndex());
		}

}
