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
package com.github.lindenb.jvarkit.liftover;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.derby.iapi.services.io.ArrayInputStream;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

public class LiftOverLoader {
	private static final Logger LOG = Logger.build(LiftOverLoader.class).make();

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
		final StringBuilder sb = new StringBuilder();
		try(CloseableIterator<String> iter = iterator(br)) {
			while(iter.hasNext()) {
				sb.append(iter.next()).append('\n');
				}			
			}
		try(InputStream in= new ArrayInputStream(sb.toString().getBytes())) {
			return new LiftOver(in,sourceName);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}

	public CloseableIterator<String> iterator(final BufferedReader in) {
		return new MyIterator(Objects.requireNonNull(in));
		}
	
	private class MyIterator extends AbstractCloseableIterator<String> {
		private final BufferedReader br;
		private boolean valid = false;
		private final Pattern splitter = Pattern.compile("\\s");
		private final Set<String> notFound1 = new TreeSet<>();
		private final Set<String> notFound2 = new TreeSet<>();
		private int num_ignored=0;

		
		MyIterator(final BufferedReader br) {
			this.br = br;
			}
		@Override
		protected String advance() {
			try {
				for(;;) {
					String line = br.readLine();
					if(line==null) {
						close();
						return null;
						}
					if(line.startsWith("chain")) {
						valid = false;
						final String chainFields[]=splitter.split(line);
						 if (chainFields.length != 13) {
							 throw new JvarkitException.TokenErrors(13, chainFields);
						 	}
						 final String fromSequenceName = chainFields[2];
						 String ctg = sourceConvert.apply(fromSequenceName);
						 if(StringUtils.isBlank(ctg)) {
							 notFound1.add(fromSequenceName);
							 ++num_ignored;
							 continue;
						 	}
						 chainFields[2] = ctg;
						 final String toSequenceName = chainFields[7];
						 ctg = targetConvert.apply(toSequenceName);
						 if(StringUtils.isBlank(ctg)) {
							 notFound2.add(toSequenceName);
							 ++num_ignored;
							 continue;
						 	}
						 chainFields[7] = ctg;
						valid=true;
						return String.join(" ", chainFields);
						}
					else if(valid)
						{
						return line;
						}
					else
						{
						++num_ignored;
						}
					}
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			try {
				if(num_ignored>0) LOG.info("number of lines skipped "+num_ignored);
				if(!notFound1.isEmpty()) LOG.info("unmatched source contigs "+notFound1.stream().collect(Collectors.joining("; ")));
				if(!notFound2.isEmpty()) LOG.info("unmatched dest contigs "+notFound2.stream().collect(Collectors.joining("; ")));
				this.br.close();
				}
			catch(final IOException err) {
				LOG.warn(err);
				}
		}
	}
	
}
