/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
import java.io.PrintWriter;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Locatable;

/**
BEGIN_DOC

```
$ echo -e "RF01\t150\t200\tA\nRF01\t190\t300\tA\nRF01\t350\t400\tA\nRF01\t150\t200\tB" |\
	java -jar dist/jvarkit.jar setfilefrombed -R src/test/resources/rotavirus_rf.fa frombed

A	RF01:151-300,RF01:351-400
B	RF01:151-200

```


END_DOC

**/
@Program(name="setfilefrombed",
description="Convert bed chrom/start/end/name sorted on 4th column to set file",
creationDate="20210125",
modificationDate="20240724",
keywords={"setfile","bed"},
jvarkit_amalgamion = true,
menu="BED Manipulation"
)
public class SetFileFromBed extends AbstractSetFileTool {
	private static final Logger LOG = Logger.build(SetFileFromBed.class).make();

	
	@Parameter(names={"--disable-interval-merge"},description="Do not merge overlapping intervals in a setFile record")
	protected boolean disable_interval_merge = false;

	
	private String bed2name(final BedLine bed) {
		if(bed.getColumnCount()<4) throw new IllegalArgumentException("Expected 4 columns but got "+bed);
		final String name= bed.get(3);
		if(StringUtils.isBlank(name)) throw new IllegalArgumentException("empty 4th columns but got "+bed);
		return name;
		};
	
	@Override
	public int doWork(final List<String> args) {
		
		try {
			final String input = oneFileOrNull(args);
			
			try(BufferedReader br = super.openBufferedReader(input)) {
				try(BedLineReader blr = new BedLineReader(br, input)) {
					blr.setValidationStringency(ValidationStringency.LENIENT);
					blr.setContigNameConverter(getContigConverter());

					try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
						final EqualIterator<BedLine> iter = new EqualIterator<BedLine>(blr.stream().iterator(),(A,B)->bed2name(A).compareTo(bed2name(B)));
						while(iter.hasNext()) {
							final List<BedLine> lines = iter.next();
							final List<Locatable> L = disable_interval_merge?
									lines.stream().map(B->Locatable.class.cast(B)).sorted(getSorter()).collect(Collectors.toList()):
									sortAndMerge(lines)
									;
							pw.print(bed2name(lines.get(0)));
							for(int i=0;i< L.size();i++) {
								pw.print(i==0?"\t":",");
								final Locatable rec = L.get(i);
								pw.print(noChr(rec.getContig()));
								pw.print(":");
								pw.print(rec.getStart());
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
		new SetFileFromBed().instanceMainWithExit(args);
	}

}
