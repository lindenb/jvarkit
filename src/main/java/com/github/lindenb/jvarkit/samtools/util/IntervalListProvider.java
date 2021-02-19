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
package com.github.lindenb.jvarkit.samtools.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Stream;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.LineIterator;
import com.github.lindenb.jvarkit.variant.variantcontext.Breakend;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.IntervalList.IntervalListCodec;
import htsjdk.variant.variantcontext.VariantContext;
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
	
	private boolean skip_unknown_contig = false;
	protected boolean enable_single_point = false;
	protected boolean enable_breakend_vcf = false;
	private SAMSequenceDictionary dictionary = null;
	private ContigNameConverter contigNameConverter = null;
	protected Predicate<VariantContext> variantPredicate = null;
	final protected String pathStr;
	
	protected IntervalListProvider(final String path) {
		this.pathStr =path;
	}
	
	protected String getRequiredPath() {
		if(StringUtils.isBlank(this.pathStr)) throw new IllegalStateException("path is unspecified");
		return this.pathStr;
	}
	
	/** use this predicate when input source is a VCF */
	public IntervalListProvider setVariantPredicate(final Predicate<VariantContext> variantPredicate) {
		this.variantPredicate = variantPredicate;
		return this;
		}
	
	/** set dictionary that can be of help to decode contig names . Can be null */
	public IntervalListProvider dictionary(final SAMSequenceDictionary dictionary) {
		this.dictionary = dictionary;
		if(this.dictionary==null) {
			this.contigNameConverter = null;
		} else {
			this.contigNameConverter = ContigNameConverter.fromOneDictionary(this.dictionary);
		}
		return this;
	}
	
	public abstract Stream<? extends Locatable> stream();

	/** convert stream to IntervalTreeMap */
	public IntervalTreeMap<Locatable> toIntervalTreeMap() {
		final IntervalTreeMap<Locatable> intervalTreeMap = new IntervalTreeMap<>();
		this.stream().map(T->new Interval(T)).forEach(I->intervalTreeMap.put(I, I));
		return intervalTreeMap;
		}
	
	
	/** convert stream to optimizedQueryIntervals for reading indexed bam */
	public QueryInterval[] optimizedQueryIntervals() {
		if(this.dictionary==null) throw new IllegalStateException("No dictionary available.");
		return this.stream().
			map(R->{
				if(R==null) throw new NullPointerException("Locatable is null at this point.");
				final SAMSequenceRecord ssr = this.dictionary.getSequence(R.getContig());
				if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(R.getContig(),this.dictionary);
				return new QueryInterval(ssr.getSequenceIndex(), R.getStart(), R.getEnd()); 
				}).
			collect(HtsCollectors.optimizedQueryIntervals());
		}
	
	public IntervalListProvider skipUnknownContigs(boolean skip_unknown_contig) {
		this.skip_unknown_contig = skip_unknown_contig;
		return this;
		}
	
	public IntervalListProvider skipUnknownContigs() {
		return this.skipUnknownContigs(true);
		}

	/** enable 'chr:pos' notation */
	public IntervalListProvider enableSinglePoint(boolean enable_single_point) {
		this.enable_single_point = enable_single_point;
		return this;
		}
	
	/** enable 'chr:pos' notation */
	public IntervalListProvider enableSinglePoint() {
		return this.enableSinglePoint(true);
		}
	
	/** in a VCF, for a BND, add the mate interval */
	public IntervalListProvider enableBreakEndInterval(boolean enable_breakend_vcf) {
		this.enable_breakend_vcf = enable_breakend_vcf;
		return this;
		}
	
	/** in a VCF, for a BND, add the mate interval */
	public IntervalListProvider enableBreakEndInterval() {
		return this.enableBreakEndInterval(true);
		}
	
	
	protected Locatable remapContig(final Locatable loc) {
		if(this.contigNameConverter==null) return loc;
		final Optional<SimpleInterval> r = this.contigNameConverter.convertToSimpleInterval(loc);
		if(!r.isPresent()) {
			if(this.skip_unknown_contig) return null;
			throw new JvarkitException.ContigNotFoundInDictionary(loc.getContig(), this.dictionary);
			}
		return r.get();
		}
	
	
	/** default provider throwing an IllegalStateException */
	public static IntervalListProvider unspecified() {
		return new IntervalListProvider("(unspecified)") {
			@Override
			public Stream<? extends Locatable> stream() {
				throw new IllegalStateException("No interval/Bed/Gtf was provided.");
				}
			};
		}

	
	/** default provider returning empty */
	public static IntervalListProvider empty() {
		return new IntervalListProvider("(empty)") {
			@Override
			public Stream<? extends Locatable> stream() {
				return Stream.empty();
				}
			};
		}
	
	/** creates a stream when reading from a java.io.Reader */
	public Stream<? extends Locatable> fromInputStream(final InputStream r,final String format) {
		return fromReader(new InputStreamReader(r),format);
	}
	
	/** creates a stream when reading from a java.io.Reader */
	public Stream<? extends Locatable> fromReader(final Reader r,final String format) {
		final  Stream<? extends Locatable> st;
		if(format.equalsIgnoreCase("bed")) {
			st = fromBedReader(r);
			} 
		else if(format.equalsIgnoreCase("gff") || format.equalsIgnoreCase("gtf")) {
			st =  fromGffGtfReader(r);
			}
		else
			{
			throw new IllegalArgumentException("unsupported format "+format+".");
			}
		return st.onClose(()->CloserUtil.close(r));
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
		else if(path.endsWith(FileExtensions.INTERVAL_LIST) || path.endsWith(FileExtensions.INTERVAL_LIST+".gz")) {
			return new ProviderIsIntervalList(path);
			}
		else
			{
			return new ProviderIsString(path);
			}
		}
	
	
	private static class ProviderIsString extends IntervalListProvider
		{
		ProviderIsString(final String path) {
			super(path);
			}
		@Override
		public Stream<? extends Locatable> stream() {
			final Function<String, Optional<SimpleInterval>> parser = IntervalParserFactory.newInstance().
					throwOnError().
					enableSinglePoint(this.enable_single_point).
					make();
			return Arrays.
					stream(getRequiredPath().split("[ ;,\t\n]+")).
					filter(S->!StringUtils.isBlank(S)).
					map(S->(Locatable)parser.apply(S).orElseThrow(()->new IllegalArgumentException("Cannot convert interval \""+S+"\""))).
					map(L->this.remapContig(L)).filter(L->L!=null)
					;
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
						open(getRequiredPath()).
						stream().
						filter(CTX->variantPredicate==null ||variantPredicate.test(CTX) ).
						flatMap(V->{
							/* enable break end ?*/
							final List<Locatable> L = new ArrayList<>();
							L.add(V);
							if(enable_breakend_vcf) L.addAll(Breakend.parse(V));
							return L.stream();
						}).
						map(L->this.remapContig(L)).
						filter(L->L!=null);
						
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		}	
	
	public Stream<? extends Locatable> fromBedReader(final Reader r) {
		final BedLineCodec codec = new BedLineCodec();
		final BufferedReader br = (r instanceof BufferedReader ? BufferedReader.class.cast(r):new BufferedReader(r));
			return br.lines().
				filter(L->!(StringUtils.isBlank(L) || L.startsWith("#"))).
				map(L->codec.decode(L)).
				filter(L->L!=null).
				map(L->this.remapContig(L)).filter(L->L!=null)
				;
		}
	
	private static class ProviderIsBed extends IntervalListProvider
		{
		ProviderIsBed(final String path) {
			super(path);
			}
		@Override
		public Stream<? extends Locatable> stream() {
			try {
				final BufferedReader r = IOUtils.openURIForBufferedReading(getRequiredPath());
				return fromBedReader(r).onClose(()->CloserUtil.close(r));
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		}	

	public Stream<? extends Locatable> fromGffGtfReader(final Reader r) {
		final BufferedReader br = (r instanceof BufferedReader ? BufferedReader.class.cast(r):new BufferedReader(r));
			return br.lines().
					filter(L->!(StringUtils.isBlank(L) || L.startsWith("#"))).
					map(L->CharSplitter.TAB.split(L)).
					map(L->new SimpleInterval(L[0],Integer.parseInt(L[3]),Integer.parseInt(L[4]))).
					map(L->this.remapContig(L)).
					filter(L->L!=null)
				;
		}
	
	private static class ProviderGtfGff extends IntervalListProvider
		{
		ProviderGtfGff(final String path) {
			super(path);
			}
		@Override
		public Stream<? extends Locatable> stream() {
			try {
				final BufferedReader r = IOUtils.openURIForBufferedReading(this.getRequiredPath());
				return fromGffGtfReader(r).onClose(()->CloserUtil.close(r));
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
				final BufferedReader r = IOUtils.openURIForBufferedReading(this.getRequiredPath());
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
		return String.valueOf(this.pathStr);
		}
}


