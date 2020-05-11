package com.github.lindenb.jvarkit.variant.vcf;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.regex.Pattern;

import org.eclipse.jetty.io.RuntimeIOException;

import com.github.lindenb.jvarkit.io.MayBeGzipInputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.bcf2.BCFVersion;

public class BcfToolsUtils {
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
			execpath = "bcftools";
			}
		return execpath;
		}
	/** return path to bcftools executable path. Never null. Default is 'bcftools' */
	public static String getBcftoolsExecutable() {
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
			if(BCF2Codec.ALLOWED_BCF_VERSION.getMajorVersion()!=(version.getMajorVersion())) {
				return true;
				}
			}
		return true;
		}
	
	}
