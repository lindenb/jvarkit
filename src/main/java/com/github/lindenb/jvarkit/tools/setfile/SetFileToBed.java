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
package com.github.lindenb.jvarkit.tools.setfile;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.setfile.SetFileRecord;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;

/**
BEGIN_DOC

TODO

END_DOC

**/
@Program(name="setfile2bed",
description="Convert setfile to bed",
creationDate="20210125",
modificationDate="20240724",
keywords={"setfile","bed"},
jvarkit_amalgamion = true,
menu="BED Manipulation"
)
public class SetFileToBed extends AbstractSetFileTool {
	private static final Logger LOG = Logger.of(SetFileToBed.class);

	
	@Parameter(names={"--bed12"},description="convert to BED12. Skip setfiles with multiple contigs")
	protected boolean to_bed12 = false;

	private int toBed12(final List<String> args) throws IOException {
		int skipped=0;
		try(CloseableIterator<SetFileRecord> iter = openSetFileIterator(args)) {
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				while(iter.hasNext()) {
					final SetFileRecord rec = iter.next();
					final Set<String> chroms = rec.getChromosomes();
					if(chroms.size()!=1) {
						skipped++;
						continue;
						}
					final List<Locatable> merged = sortAndMerge(rec);
					final String contig = chroms.iterator().next();
					final int chromStart;
					final int chromEnd;
					pw.print(noChr(contig));
					pw.print("\t");
					pw.print(chromStart = merged.stream().mapToInt(X->X.getStart()).min().getAsInt()-1);
					pw.print("\t");
					pw.print(chromEnd = merged.stream().mapToInt(X->X.getEnd()).max().getAsInt());
					pw.print("\t");
					pw.print(rec.getName());
					pw.print("\t");
					pw.print((int)((merged.stream().mapToInt(L->L.getLengthOnReference()).sum()/(double)(chromEnd-chromStart))*1000));
					pw.print("\t");
					pw.print(".");//no strand
					pw.print("\t");
					pw.print(chromStart);
					pw.print("\t");
					pw.print(chromEnd);
					pw.print("\t");
					pw.print("0,0,255");
					pw.print("\t");
					pw.print(rec.size());
					pw.print("\t");
					pw.print(merged.stream().map(T->String.valueOf(T.getLengthOnReference())).collect(Collectors.joining(",")));
					pw.print("\t");
					pw.print(merged.stream().map(T->String.valueOf((T.getStart()-1)-chromStart)).collect(Collectors.joining(",")));
					pw.println();
					}
				pw.flush();
				}
			}
		if(skipped>0) {
			LOG.warn("skipped "+skipped+" records with multiple chromosomes.");
			}
		return 0;
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
	
	@Override
	public int doWork(final List<String> args) {
		
		try {
			return this.to_bed12?toBed12(args):toBed(args);
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	public static void main(final String[] args) {
		new SetFileToBed().instanceMainWithExit(args);
	}

}
