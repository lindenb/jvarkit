/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.*;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.*;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.*;
import htsjdk.samtools.util.IntervalList;

/**
BEGIN_DOC

## Example

```
```


END_DOC

 */
@Program(
		name="bedmergecnv",
		description="TODO.",
		keywords={"bed","chromosome","contig"},
		creationDate="20200330",
		modificationDate="20200330"
		)
public class BedMergeCnv
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BedMergeCnv.class).make();
	
	@Parameter(names={"-o","--out"},description="TODO")
	private Path outputFile= null;
	@Parameter(names={"-f","--f"},description="TODO")
	private double fraction = 0.9;

	
	private class Cluster implements Locatable{
		private final Locatable base;
		private final List<Locatable> intervals = new ArrayList<>();
		Cluster(final Locatable	base) {
			this.base = base;
			this.intervals.add(base);
			}		
		@Override
		public String getContig() {
			return base.getContig();
			}
		@Override
		public int getStart() {
			return this.intervals.stream().mapToInt(R->R.getStart()).min().orElse(-1);
			}
		@Override
		public int getEnd() {
			return this.intervals.stream().mapToInt(R->R.getEnd()).max().orElse(-1);
			}

		private boolean test(final Locatable R1,final Locatable R2) {
			if(!R1.overlaps(R2)) return false;
			final int x1 = Math.max(R1.getStart(),R2.getStart());
			final int x2 = Math.min(R1.getEnd(),R2.getEnd());
			double L = CoordMath.getLength(x1,x2);
			double L1 = R1.getLengthOnReference();
			double L2 = R2.getLengthOnReference();
			if(L/L1 < fraction) return false;
			if(L/L2 < fraction) return false;
			return true;
			}

		public boolean add(final Locatable si) {
			if(!this.base.contigsMatch(si)) return false;
			if(this.intervals.stream().anyMatch(R->!test(R,si))) return false;
			this.intervals.add(si);
			return true;
			}
		}
	

		
	@Override
	public int doWork(final List<String> args) {
		final BedLineCodec codec = new BedLineCodec();
		PrintWriter w=null;
		try
			{
			final IntervalTreeMap<Cluster> clusters = new IntervalTreeMap<>();
			try(BufferedReader br = new BufferedReader(new InputStreamReader(stdin()))) {
				String line;
				while((line=br.readLine())!=null) {
					if( BedLine.isBedHeader(line) ) continue;
					final BedLine rec = codec.decode(line);
					if(rec==null) continue;
					boolean found=false;
					Interval key = new Interval(rec);
					for(final Cluster bc: clusters.getOverlapping(key)) {
						if(bc.add(rec)) { found=true;}
						}
					if(!found) {
						clusters.put(key,new Cluster(rec));
						}
					}
				}
			w = new PrintWriter(stdout());

			for(final Cluster bc: clusters.values()) {
				w.print(bc.getContig());
				w.print("\t");
				w.print(bc.getStart()-1);
				w.print("\t");
				w.print(bc.getEnd());
				w.print("\t");
				w.print(bc.intervals.size());
				w.println();
				}
			w.flush();
			w.close();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			}
		}
	

	public static void main(final String[] args)
		{
		new BedMergeCnv().instanceMainWithExit(args);
		}
	}
