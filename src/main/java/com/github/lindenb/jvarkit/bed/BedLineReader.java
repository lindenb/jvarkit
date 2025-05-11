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
package com.github.lindenb.jvarkit.bed;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.function.UnaryOperator;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

/**
 * A simple BedLine reader/Iterator
 * @author lindenb
 */
public class BedLineReader extends AbstractCloseableIterator<BedLine>  {
private static final Logger LOG=Logger.of(BedLineReader.class);
public static final String OPT_DESC="A Bed file: (CHROM)<tab>(START 0-based)<tab>(END)[<tab>otherfields...].";
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

public BedLineReader setValidationStringency(final ValidationStringency stringency) {
	this.codec.setValidationStringency(stringency);
	return this;
	}

public BedLineReader setContigNameConverter(final UnaryOperator<String> chromosomeConverter) {
	this.codec.setContigNameConverter(chromosomeConverter);
	return this;
	}

public  QueryInterval[] optimizeIntervals(final SAMSequenceDictionary dict) {
	Objects.requireNonNull(dict);
	return QueryInterval.optimizeIntervals(stream().
		map(B->{
			int tid = dict.getSequenceIndex(B.getContig());
			if(tid<0) {
				final String msg = "Chromosome is not in dictionary. Skipping : "+B.getContig();
				switch (this.codec.getValidationStringency()) {
				case STRICT:
					throw new IllegalArgumentException(msg);
				case LENIENT:
					LOG.warn(msg);
					return null;
				default:
					return null;
					}
				}
			return new QueryInterval(tid, B.getStart(), B.getEnd());
			}).
		filter(X->X!=null).
		toArray(X->new QueryInterval[X]));
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
				/* no need to display a message because validation stringency set in codec
				 * LOG.warn("Cannot read bed entry \""+line.replace("\t", "\\t")+"\" in "+sourceName);
				 */
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

/** convert all items to an IntervalTreeMap<X> . Each Bedline is converted to X using transform.
 * Items are ignored if 'transform' returns null
 * @param transform convert BedLine to X. if x==null, item is ignored.
 * @return
 */
public <X> IntervalTreeMap<X> toIntervalTreeMap(final Function<BedLine,X> transform) {
	final IntervalTreeMap<X> map = new IntervalTreeMap<>();
	while(this.hasNext()) {
		final BedLine bl = this.next();
		final X x  = transform.apply(bl);
		if(x==null) continue;
		final Interval r = new Interval(bl);
		map.put(r, x);
		}
	return map;
	}

/** convert all items to an IntervalTreeMapBedLine>
 */
public IntervalTreeMap<BedLine> toIntervalTreeMap() {
return toIntervalTreeMap(X->X);
}


/** convert all items to an IntervalTreeMap<List<X>> . Each Bedline is converted to X using transform.
 * Items are ignored if 'transform' returns null
 * @param transform convert BedLine to X. if x==null, item is ignored.
 * @return
 */
public <X> IntervalTreeMap<List<X>> toIntervalListTreeMap(final Function<BedLine,X> transform) {
	return toIntervalCollectionTreeMap(transform,()->new ArrayList<X>());
	}

/** convert all items to an IntervalTreeMap<List<X>> . Each Bedline is converted to X using transform.
 * Items are ignored if 'transform' returns null
 * @param transform convert BedLine to X. if x==null, item is ignored.
 * @return
 */
public <X,C extends Collection<X>> IntervalTreeMap<C> toIntervalCollectionTreeMap(final Function<BedLine,X> transform,final Supplier<C> collectionSupplier) {
	final IntervalTreeMap<C> map = new IntervalTreeMap<>();
	while(this.hasNext()) {
		final BedLine bl = this.next();
		final X x  = transform.apply(bl);
		if(x==null) continue;
		final Interval r = new Interval(bl);
		C col = map.get(r);
		if(col==null) {
			col = collectionSupplier.get();
			map.put(r, col);
			}
		col.add(x);
		}
	return map;
	}


/** convert all items to an IntervalTreeMapBedLine>
 */
public IntervalTreeMap<List<BedLine>> toIntervalListTreeMap() {
return toIntervalListTreeMap(X->X);
}

@Override
public String toString() {
	return "BedLineReader("+sourceName+")";
	}
}	
