/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.ucsc;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.UnaryOperator;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;

/**
 * Convert, on the fly, the content of an Chain File, converting the chromosomes
 */
public class LiftOverChainInputStream extends InputStream {
	private final Pattern splitter = Pattern.compile("\\s");
	private BufferedReader delegate;
	private ByteArrayInputStream buffer=null;
	private final UnaryOperator<String> convertSource;
	private final UnaryOperator<String> convertDest;
	private final Set<String> notFound1 = new TreeSet<>();
	private final Set<String> notFound2 = new TreeSet<>();
	private int num_ignored=0;
	private boolean valid = false;
	private final Charset charset = Charset.forName("UTF-8");
	private String sourceName="undefined";

	
	public LiftOverChainInputStream(BufferedReader delegate, 
			UnaryOperator<String> convertSource,
			UnaryOperator<String> convertDest) throws IOException  {
		this.delegate=delegate;
		this.convertSource=convertSource==null?A->A:convertSource;
		this.convertDest=convertDest==null?A->A:convertDest;
		}
	
	public LiftOverChainInputStream(InputStream delegate, 
			UnaryOperator<String> convertSource,
			UnaryOperator<String> convertDest) throws IOException  {
		this(new BufferedReader(new InputStreamReader(delegate, Charset.forName("UTF-8"))),convertSource,convertDest);
		}
	
	private static InputStream _openPath(Path path) throws IOException {
		IOUtil.assertFileIsReadable(path);
		return  IOUtils.mayBeGzippedInputStream(Files.newInputStream(path));
		}
	private static InputStream _open(String path) throws IOException {
		String src1 = LiftOverChain.convert(path);
		final InputStream is;
		if(IOUtil.isUrl(src1)) {
			is= IOUtils.mayBeGzippedInputStream(new URL(src1).openStream());
			}
		else
			{
			is =  _openPath(Paths.get(src1));
			}
		return is;
		}
	
	public LiftOverChainInputStream(String path, UnaryOperator<String> convertSource,UnaryOperator<String> convertDest) throws IOException {
		this(_open(path),convertSource,convertDest);
		this.sourceName=LiftOverChain.convert(path);
		}
	public LiftOverChainInputStream(Path path, UnaryOperator<String> convertSource,UnaryOperator<String> convertDest) throws IOException {
		this(_openPath(path),convertSource,convertDest);
		this.sourceName=path.toString();
		}
	
	private boolean fillBuffer()  throws IOException {
		if(this.delegate==null)return false;
		this.buffer=null;
		for(;;) {
			final String line = this.delegate.readLine();
			if(line==null) {
				this.delegate.close();
				this.delegate=null;
				return false;
				}
			if(line.startsWith("chain")) {
				valid = false;
				final String chainFields[]=splitter.split(line);
				 if (chainFields.length != 13) {
					 throw new JvarkitException.TokenErrors(13, chainFields);
				 	}
				 final String fromSequenceName = chainFields[2];
				 String ctg = this.convertSource.apply(fromSequenceName);
				 if(StringUtils.isBlank(ctg)) {
					 this.notFound1.add(fromSequenceName);
					 ++num_ignored;
					 continue;
				 	}
				 chainFields[2] = ctg;
				 final String toSequenceName = chainFields[7];
				 ctg = this.convertDest.apply(toSequenceName);
				 if(StringUtils.isBlank(ctg)) {
					 this.notFound2.add(toSequenceName);
					 ++num_ignored;
					 continue;
				 	}
				 chainFields[7] = ctg;
				valid=true;
				final String join =  String.join(" ", chainFields) + "\n";
				this.buffer = new ByteArrayInputStream(join.getBytes(charset));
				return true;
				}
			else if(valid==false) {
				++this.num_ignored;
				}
			else
				{
				final String join =  line + "\n";
				this.buffer = new ByteArrayInputStream(join.getBytes(charset));
				return true;
				}
			}
		}
	
	@Override
	public int read() throws IOException {
		for(;;) {
			if(buffer==null) {
				fillBuffer();
				if(buffer==null) return -1;
				}
			int c= buffer.read();
			if(c!=-1) return c;
			buffer=null;
			}
		}
	@Override
	public void close() throws IOException {
		if(this.delegate!=null) {
			delegate.close();
			this.delegate=null;
			}
		this.buffer=null;
		}
	public void log(final Logger LOG) {
		LOG.info("source: "+this.sourceName);
		LOG.info("number of lines skipped "+num_ignored);
		LOG.info("unmatched source contigs "+String.join("; ",notFound1));
		LOG.info("unmatched dest contigs "+String.join("; ",notFound2));
		}
	}
