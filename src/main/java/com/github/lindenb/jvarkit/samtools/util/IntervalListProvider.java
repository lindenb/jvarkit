/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.samtools.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.iterator.LineIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.IntervalList.IntervalListCodec;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public abstract class IntervalListProvider {
	public static final String OPT_DESC="A source of intervals. The following suffixes are recognized: "
			+ "vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz."
			+ "Otherwise it could be an empty string (no interval) or a list of plain interval separated by '[ \\t\\n;,]' ";
	public static final class StringConverter implements IStringConverter<IntervalListProvider> {
		@Override
		public IntervalListProvider convert(final  String s) {
			return IntervalListProvider.from(s);
			}
		}
	
	final protected String path;
	
	protected IntervalListProvider(final String path) {
		this.path =path;
	}
	
	
	public abstract Stream<? extends Locatable> stream();

	
	public static IntervalListProvider empty() {
		return new IntervalListProvider("(empty)") {
			@Override
			public Stream<? extends Locatable> stream() {
				return Stream.empty();
				}
			};
		}
	
	public static final IntervalListProvider from(final String path) {
		if(StringUtils.isBlank(path)) {
			return empty();
			}
		if(FileExtensions.VCF_LIST.stream().anyMatch(S->path.endsWith(S))) {
			return new ProviderIsVcf(path);
			}
		else if(path.endsWith(FileExtensions.BED) || path.endsWith(".bed.gz")) {
			return new ProviderIsBed(path);
			}
		else if(StringUtils.endsWith(path, ".gtf",".gtf.gz",".gff",".gff.gz")) {
			return new ProviderGtfGff(path);
			}
		else if(path.endsWith(FileExtensions.INTERVAL_LIST)) {
			return new ProviderIsIntervalList(path);
			}
		else
			{
			return new ProviderIsString(path);
			}
		}
	
	
	private static class ProviderIsString extends IntervalListProvider
		{
		private final List<? extends Locatable> array;
		ProviderIsString(final String path) {
			super(path);
			
			final Function<String, Optional<SimpleInterval>> parser = IntervalParserFactory.newInstance().make();
			this.array =Arrays.
					stream(path.split("[ ;,\t\n]+")).
					filter(S->!StringUtils.isBlank(S)).
					map(S->(Locatable)parser.apply(S).orElseThrow(()->new IllegalArgumentException("Cannot convert interval"))).
					collect(Collectors.toList());
			
			}
		@Override
		public Stream<? extends Locatable> stream() {
			return array.stream();
			}
		}	

	
	private static class ProviderIsVcf extends IntervalListProvider
		{
		ProviderIsVcf(final String path) {
			super(path);
			}
		@Override
		public Stream<? extends Locatable> stream() {
			try {
				return new VCFIteratorBuilder().
						open(path).
						stream();
						
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		}	
	
	private static class ProviderIsBed extends IntervalListProvider
		{
		ProviderIsBed(final String path) {
			super(path);
			}
		@Override
		public Stream<? extends Locatable> stream() {
			try {
				final BedLineCodec codec = new BedLineCodec();
				final BufferedReader r = IOUtils.openURIForBufferedReading(this.path);
				return r.lines().
						filter(L->!(StringUtils.isBlank(L) || L.startsWith("#"))).
						map(L->codec.decode(L)).
						filter(L->L!=null).
						onClose(()->CloserUtil.close(r));
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		}	

	private static class ProviderGtfGff extends IntervalListProvider
		{
		ProviderGtfGff(final String path) {
			super(path);
			}
		@Override
		public Stream<? extends Locatable> stream() {
			try {
				final BufferedReader r = IOUtils.openURIForBufferedReading(this.path);
				return r.lines().
						filter(L->!(StringUtils.isBlank(L) || L.startsWith("#"))).
						map(L->CharSplitter.TAB.split(L)).
						map(L->new SimpleInterval(L[0],Integer.parseInt(L[3]),Integer.parseInt(L[4]))).
						onClose(()->CloserUtil.close(r));
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		}	
	private static class ProviderIsIntervalList extends IntervalListProvider
		{
		ProviderIsIntervalList(final String path) {
			super(path);
			}
		
		@Override
		public Stream<? extends Locatable> stream() {
			try {
				final IntervalListCodec codec=new IntervalListCodec();
				final BufferedReader r = IOUtils.openURIForBufferedReading(this.path);
				final LineIterator lr = new LineIterator(r);
				codec.readActualHeader(lr);
				return lr.stream().
						map(L->codec.decode(L)).
						filter(L->L!=null).
						onClose(()->CloserUtil.close(r));
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		}	

	@Override
	public String toString() {
		return this.path;
		}
}


