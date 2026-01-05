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

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

public class LiftOverLoader {
	private static final Logger LOG = Logger.of(LiftOverLoader.class);
	public static final String OPT_DICT_R1="Source of chromosome names to convert chromosome names ('chr1'->'1') for the source assembly. Could be a dict, a fai, etc...";
	public static final String OPT_DICT_R2="Source of chromosome names to convert chromosome names ('chr1'->'1') for thetarget assembly. Could be a dict, a fai, etc...";
	private ContigNameConverter sourceConvert = ContigNameConverter.getIdentity();
	private ContigNameConverter targetConvert = ContigNameConverter.getIdentity();

	
	public LiftOverLoader setSourceDictionary(final SAMSequenceDictionary dict) {
		return setSourceMapper(dict==null?null:ContigNameConverter.fromOneDictionary(dict));
	}
	public LiftOverLoader setSourceMapper(final ContigNameConverter convert) {
		this.sourceConvert = convert==null?ContigNameConverter.getIdentity():convert;
		return this;
	}
	
	public LiftOverLoader setTargetDictionary(final SAMSequenceDictionary dict) {
		return setTargetMapper(dict==null?null:ContigNameConverter.fromOneDictionary(dict));
	}
	
	public LiftOverLoader setTargetMapper(final ContigNameConverter convert) {
		this.targetConvert = convert==null?ContigNameConverter.getIdentity():convert;
		return this;
	}

	public LiftOver load(final String urlOrFile) {
		if(IOUtil.isUrl(urlOrFile)) {
			try(BufferedReader in= IOUtils.openURIForBufferedReading(urlOrFile)) {
				return load(in,urlOrFile);
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		else
			{
			return load(Paths.get(urlOrFile));
			}
		}

	
	public LiftOver load(final Path path) {
		try(BufferedReader br = IOUtil.openFileForBufferedReading(path)) {
			return load(br,path.toString());
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	public LiftOver load(final BufferedReader br,final String sourceName) {
		try(InputStream in= new LiftOverChainInputStream(br, sourceConvert, targetConvert)) {
			return new LiftOver(in,sourceName);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}

}
