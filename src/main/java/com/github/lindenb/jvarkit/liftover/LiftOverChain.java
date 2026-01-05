/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.liftover;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.function.UnaryOperator;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.IOUtil;

public class LiftOverChain {
	public static final String OPT_DESC="LiftOver chain file. Can be a local chain file, a URL 'https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToCriGri1.over.chain.gz', or a chain identifier like 'hg19ToHg38'.";
	/** convert lifover basename chain e.g. hg19ToCriGri1 to ucsc URL ( https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToCriGri1.over.chain.gz ) 
	 * or return the original String
	 * @param src like hg19ToCriGri1
	 * @return
	 */
	public static String convert(final String src) {
		if(src==null) return null; // yes, used in Launcher.oneFileOrNull
		if(src.contains("/") || src.contains(".")) return src;
		int to = src.indexOf("To");
		if(to<=0 || to+4>=src.length()) return src;
		final Pattern pat = Pattern.compile("[a-z][a-zA-Z]*[0-9]+");
		final String build1 = src.substring(0,to);
		if(!pat.matcher(build1).matches()) return src;
		final String build2 = src.substring(to+2);
		if(!pat.matcher(build2.toLowerCase()).matches()) return src;
		if(!Character.isUpperCase(build2.charAt(0))) return src;
		return "https://hgdownload.soe.ucsc.edu/goldenpath/"+build1+"/liftOver/"+src+".over.chain.gz";
		}
	/** 
	 * load lifover . can be a url, a path, or a identifier like hg19ToHg38
	 */
	public static LiftOver load(final String src0) throws IOException {
		final String src1 = LiftOverChain.convert(src0);
		if(IOUtil.isUrl(src1)) {
			try(InputStream in = IOUtils.mayBeGzippedInputStream(IOUtils.toURL(src1).openStream())) {
				return new LiftOver(in,src1);
				}
			}
		else
			{
			final Path p = Paths.get(src1);
			IOUtil.assertFileIsReadable(p);
			try(InputStream in = IOUtils.mayBeGzippedInputStream(Files.newInputStream(p))) {
				return new LiftOver(in,src1);
				}
			}
		}
	
	public static LiftOver load(final String src,SAMSequenceDictionary dictSrc, SAMSequenceDictionary dictDest) throws IOException {
		return load(src, ContigNameConverter.fromOneDictionary(dictSrc), ContigNameConverter.fromOneDictionary(dictDest));
		}
	
	public static LiftOver load(final String src,UnaryOperator<String> convertSrc, UnaryOperator<String> convertDest) throws IOException {
		try(InputStream in = new LiftOverChainInputStream(src, convertSrc, convertDest)) {
			return new LiftOver(in,convert(src));
			}
		}
}
