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
package com.github.lindenb.jvarkit.bed;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

/**
 * A simple BedLine reader/Iterator
 * @author lindenb
 */
public class BedLineReader extends AbstractCloseableIterator<BedLine>  {
private static final Logger LOG = Logger.build(BedLineReader.class).make();
private final BedLineCodec codec = new BedLineCodec();
private final BufferedReader br;
private final String sourceName;
public BedLineReader(final Path path) {
	this(IOUtil.openFileForBufferedReading(path),path.toString());
	}
public BedLineReader(final File path) {
	this(IOUtil.openFileForBufferedReading(path),path.toString());
	}
public BedLineReader(final InputStream in,final String sourceName) throws IOException{
	this(new BufferedReader(new InputStreamReader(in,"UTF-8")),sourceName);
	}

public BedLineReader(final BufferedReader br,final String sourceName) {
	this.br = br;
	this.sourceName = sourceName==null?"undefined":sourceName;
	}

@Override
protected BedLine advance() {
	try {
		for(;;) {
			final String line = br.readLine();
			if(line==null) break;
			if(StringUtil.isBlank(line) || BedLine.isBedHeader(line)) {
				continue;
				}
			final BedLine bed = this.codec.decode(line);
			if (bed==null) {
				LOG.warn("Cannot read bed entry \""+line.replace("\t", "\\t")+"\" in "+sourceName);
				continue;
				}
			return bed;
			}
	} catch(final IOException err) {
		throw new RuntimeIOException(err);
		}
	return null;
	}

public void close() {
	try{if(br!=null) br.close();} catch(final IOException err) {}
	}

@Override
public String toString() {
	return "BedLineReader("+sourceName+")";
	}
}	
