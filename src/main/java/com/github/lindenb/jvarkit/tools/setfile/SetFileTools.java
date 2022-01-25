/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.setfile;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.setfile.SetFileReaderFactory;
import com.github.lindenb.jvarkit.setfile.SetFileRecord;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC


END_DOC

**/
@Program(name="setfiletools",
description="Utilities for the setfile format",
creationDate="20210125",
modificationDate="20210125",
keywords={"setfile"}
)

public class SetFileTools extends Launcher {
	private static final Logger LOG = Logger.build(SetFileTools.class).make();
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxRef = null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@Parameter(names={"-t","--trim-chr"},description="Remove chr prefix in chromosome names")
	protected boolean trim_chr_prefix = false;
	
	protected int n_jobs = -1;

	private abstract class AbstractIntersectIterator extends AbstractCloseableIterator<SetFileRecord> {
		final CloseableIterator<SetFileRecord> delegate;
		AbstractIntersectIterator(final CloseableIterator<SetFileRecord> delegate) {
			this.delegate = delegate;
			}
		abstract boolean overlap(Locatable loc);
		@Override
		protected final SetFileRecord advance() {
			while(this.delegate.hasNext()) {
				SetFileRecord rec = this.delegate.next();
				boolean foundCtx = false;
				for(int i=0;!foundCtx && i< rec.size();i++) {
					if(overlap(rec.get(i))) {
						foundCtx=true;
						break;
						}
					}
				if(!foundCtx) continue;
				return rec;
				}
			return null;
			}
		@Override
		public void close() {
			delegate.close();
			}
		}
	
	private class IntersectVcfIterator extends AbstractIntersectIterator {
		VCFReader vcfReader;
		IntersectVcfIterator(final Path vcfPath,final CloseableIterator<SetFileRecord> delegate) {
			super(delegate);
			}
		
		@Override
		boolean overlap(final Locatable loc) {
			try(CloseableIterator<VariantContext> iter2 = vcfReader.query(loc)) {
				if(iter2.hasNext())  return true;
				}
			return false;
			}
		
		@Override
		public void close() {
			delegate.close();
			super.close();
			}
		}
	
	private static class Cluster {
		final List<SetFileRecord> records = new ArrayList<>();
		long sum_length = 0L;
		void add(SetFileRecord rec) {
			this.records.add(rec);
			this.sum_length += rec.getSumOfLengthOnReference();
		}
		double getMean(final SetFileRecord ifAddThisRed) {
			double l= this.sum_length + ifAddThisRed.getSumOfLengthOnReference();
			return l/(records.size()+1);
		}
	}

	private enum Action {
		tobed,
		frombed
		};
		
	private String noChr(final String contig) {
		if(trim_chr_prefix && contig.toLowerCase().startsWith("chr")) {
			return contig.substring(3);
		}
		return contig;
	}
	
	private CloseableIterator<SetFileRecord> openSetFileIterator(final List<String> args) throws IOException {
		final SetFileReaderFactory srf = new SetFileReaderFactory(SequenceDictionaryUtils.extractRequired(this.faidxRef));
		return srf.open(pathOrNull);
		}
	
	
	private int intersectBed() {
		return 0;
	}
	
	private int makeClusters(final List<String> args) throws IOException {
		final List<Cluster> clusters = new ArrayList<>();
	
		try(CloseableIterator<SetFileRecord> iter = openSetFileIterator(args)) {
			final List<SetFileRecord> records = iter.stream().
				sorted((A,B)->Integer.compare(B.getSumOfLengthOnReference(), A.getSumOfLengthOnReference())).
				collect(Collectors.toList());
			while(!records.isEmpty()) {
				if(n_jobs>0) {
					if(clusters.size() < this.n_jobs) {
						final Cluster c = new Cluster();
						c.add(records.remove(0));
						}
					else {
						Cluster best =null;
						for(int i=0;i< clusters.size();i++) {
							if(best==null ) {
								best= clusters.get(i);
							}
						}
						best.add(records.remove(0));
					}
				}
			}
		}
	}

	
	private int toBed(final List<String> args) throws IOException {
		try(CloseableIterator<SetFileRecord> iter = openSetFileIterator(args)) {
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				while(iter.hasNext()) {
					final SetFileRecord rec = iter.next();
					for(Locatable loc : rec) {
						pw.print(noChr(loc.getContig()));
						pw.print("\t");
						pw.print(loc.getStart()-1);
						pw.print("\t");
						pw.print(loc.getEnd());
						pw.print("\t");
						pw.print(rec.getName());
						pw.println();
						}
					}
				pw.flush();
				}
			}
		return 0;
		}
	
	private int fromBed(final List<String> args) throws IOException {
		final String input = oneFileOrNull(args);
		final Function<BedLine,String> bed2name = bed->{
			if(bed.getColumnCount()<4) throw new IllegalArgumentException("Expected 4 columns but got "+bed);
			return bed.get(3);
			};
		try(BufferedReader br = super.openBufferedReader(input)) {
			try(BedLineReader blr = new BedLineReader(br, input)) {
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					EqualIterator<BedLine> iter = new EqualIterator<BedLine>(blr.stream().iterator(),(A,B)->bed2name.apply(A).compareTo(bed2name.apply(B)));
					while(iter.hasNext()) {
						final List<BedLine> lines = iter.next();
						pw.print(lines.get(0).getContig());
						for(int i=0;i< lines.size();i++) {
							pw.print(i==0?"\t":",");
							final BedLine rec = lines.get(i);
							pw.print(noChr(rec.getContig()));
							pw.print(":");
							pw.print(rec.getStart()+1);
							pw.print("-");
							pw.print(rec.getEnd());
							}
						pw.println();
						}
					iter.close();
					pw.flush();
					}
				}
			}
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args0) {
		if(args0.isEmpty()) {
			LOG.error("action parameter is missing.");
			return -1;
			}
		try {
			final Action action = Action.valueOf(args0.get(0));
			final List<String> args = args0.subList(1, args0.size());
			switch(action) {
				case frombed: return fromBed(args);
				case tobed: return toBed(args);
				default: LOG.error("not implemented "+action);return -1;
				}
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new SetFileTools().instanceMainWithExit(args);
	}

}
