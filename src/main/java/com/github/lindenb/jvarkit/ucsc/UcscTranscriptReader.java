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


import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.sql.parser.SchemaParser;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;

public class UcscTranscriptReader implements FeatureReader<UcscTranscript> {
	
	public static final String OPT_DESC=
			"Transcrips as genpred format https://genome.ucsc.edu/FAQ/FAQformat.html#format9  ."
			+ " The genePred format is a compact alternative to GFF/GTF because one transcript is described using only one line."
			+ "	Beware chromosome names are formatted the same as your REFERENCE. A typical KnownGene file is http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz ."
			+ "If you only have a gff file, you can try to generate a knownGene file with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)";
	
	public static final String OPT_TABIX_DESC= OPT_DESC+ ". File must be sorted on chromosome/txStart, compressed with bgzip and indexed as 0-based with tabix";

	
	public static final String SQL_DESC="Each instance of transcript can be associated to a .sql schema to help the software to decode the semantics of the columns. Eg.: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV20.sql ";
	
private final FeatureReader<UcscTranscript> delegate;

public	UcscTranscriptReader(final Path path) {
	this(path.toString());
	}

public	UcscTranscriptReader(final String uri) {
	if(!uri.endsWith(UcscTranscriptCodec.FILE_SUFFIX)) {
		throw new RuntimeIOException("path '"+uri+"' should end with '"+UcscTranscriptCodec.FILE_SUFFIX+"'");
		}
	final String sqlpath = uri.substring(0,uri.length()-UcscTranscriptCodec.FILE_SUFFIX.length())+".sql";
	SchemaParser.Table table;
	try(InputStream sqlin = IOUtils.openURIForReading(sqlpath)) {
		table = SchemaParser.parseTable(sqlin);
		}
	catch(Throwable err) {
		throw new RuntimeIOException("Every genePred file must be compressed with bgzip, indexed with tabix and associated with a .sql file. "
				+ "Eg: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV20.sql and http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV20.txt.gz", err);
		}
	
	final UcscTranscriptCodec codec = new UcscTranscriptCodec(table);
	this.delegate= AbstractFeatureReader.getFeatureReader(uri, codec);
	}
	
@Override
public void close() {
	try {this.delegate.close();}
	catch(final Throwable err) {}
	}

public CloseableTribbleIterator<UcscTranscript> query(final Locatable loc) throws IOException {
	return this.query(loc.getContig(),loc.getStart(),loc.getEnd());
	}

@Override
public CloseableTribbleIterator<UcscTranscript> query(String chr, int start, int end) throws IOException {
	return this.delegate.query(chr,start,end);
	}

@Override
public CloseableTribbleIterator<UcscTranscript> iterator() throws IOException {
	return this.delegate.iterator();
	}

@Override
public List<String> getSequenceNames() {
	return this.delegate.getSequenceNames();
	}

@Override
public Object getHeader() {
	return this.delegate.getHeader();
	}
@Override
public boolean isQueryable() {
	return true;
	}

public Map<String,List<UcscTranscript>> queryAsGeneMap(final Locatable loc) throws IOException {
	Set<String> geneNames= null;
	int minPos=0;
	int maxPos=0;
	try(CloseableTribbleIterator<UcscTranscript> iter = this.query(loc)) {
		while(iter.hasNext()) {
			final UcscTranscript tr= iter.next();
			if(StringUtils.isBlank(tr.getName2())) continue;
			if(geneNames==null) {
				geneNames = new HashSet<>();
				minPos = tr.getStart();
				maxPos = tr.getEnd();
				}
			else
				{
				minPos = Math.min(minPos, tr.getStart());
				maxPos = Math.min(maxPos, tr.getEnd());
				}
			geneNames.add(tr.getName2());
			}
		}
	if(geneNames==null) return Collections.emptyMap();
	final Map<String,List<UcscTranscript>> hash = new HashMap<String, List<UcscTranscript>>(geneNames.size());
	try(CloseableTribbleIterator<UcscTranscript> iter = this.query(loc)) {
		while(iter.hasNext()) {
			final UcscTranscript tr= iter.next();
			if(StringUtils.isBlank(tr.getName2())) continue;
			if(!geneNames.contains(tr.getName2())) continue;
			List<UcscTranscript> L =hash.get(tr.getName2());
			if(L==null) {
				L = new ArrayList<>();
				hash.put(tr.getName2(), L);
				}
			L.add(tr);
			}
		}
	return hash;
	}


public IntervalTreeMap<List<UcscTranscript>> readAsIntervalTreeMap() throws IOException {
	return readAsIntervalTreeMap(T->true);
	}
public IntervalTreeMap<List<UcscTranscript>> readAsIntervalTreeMap(final Predicate<UcscTranscript> predicate) throws IOException {
	final IntervalTreeMap<List<UcscTranscript>> tm = new IntervalTreeMap<>();
	
	try(CloseableTribbleIterator<UcscTranscript> iter =iterator()) {
		while(iter.hasNext()) {
		final UcscTranscript tr= iter.next();
			if(tr==null ||!predicate.test(tr)) continue;
			final Interval r = new Interval(tr);
			List<UcscTranscript> L = tm.get(r);
			if(L==null) {
				L=new ArrayList<>();
				tm.put(r,L);
				}
			L.add(tr);
			}
		}
	
	return tm;
	}

}
