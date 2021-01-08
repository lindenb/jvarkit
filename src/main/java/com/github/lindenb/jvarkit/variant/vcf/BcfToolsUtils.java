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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.regex.Pattern;


import com.github.lindenb.jvarkit.io.MayBeGzipInputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.bcf2.BCFVersion;

public class BcfToolsUtils {
	/* path to bcftools OR BLANK if not found */
	
	private static final String BCFTOOLS_EXE =get_bcftools_exe() ;
	
	private static String get_bcftools_exe() {
		String execpath =  System.getProperty("bcftools.exe");
		
		if(StringUtils.isBlank(execpath)) {
			for(String dir:StringUtils.ifBlank(System.getenv("PATH"),"").split(Pattern.quote(File.pathSeparator))) {
				if(StringUtils.isBlank(dir)) continue;
				final File exe = new File(dir,"bcftools");
				if(!exe.exists()) continue;
				if(exe.isDirectory()) continue;
				if(!exe.canExecute()) continue;
				execpath = exe.getPath();
				break;
				}
			}
		if(StringUtils.isBlank(execpath)) {
			execpath = null;
			}
		return execpath;
		}

	/** tell wheter bcftools is available  */
	public static boolean isBcfToolsAvailable() {
		return !StringUtils.isBlank(BCFTOOLS_EXE);
		}
	
	/** return path to bcftools executable path. Never null. Default is 'bcftools' */
	public static String getBcftoolsExecutable() {
		if(!isBcfToolsAvailable()) {
			return "bcftools";
			}
		return BCFTOOLS_EXE;
		}
	/** tell wheter file should be opened with bcftool */
	public static boolean isBcfToolsRequired(final Path path) {
		if(path.getFileName().toString().endsWith(FileExtensions.BCF)) {
			final BCFVersion version;
			try(InputStream is = new MayBeGzipInputStream(path)) {
				version = BCFVersion.readBCFVersion(is);
				if(version==null) throw new IllegalArgumentException("File "+path+" doesn't look like a BCF file.");
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			if(BCF2Codec.ALLOWED_BCF_VERSION.getMajorVersion()!=(version.getMajorVersion()) ||
				(BCF2Codec.ALLOWED_BCF_VERSION.getMinorVersion()< version.getMinorVersion())) {
				return true;
				}
			}
		return false;
		}

	}
