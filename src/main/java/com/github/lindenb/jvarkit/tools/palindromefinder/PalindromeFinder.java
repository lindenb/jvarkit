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
package com.github.lindenb.jvarkit.tools.palindromefinder;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequence;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequenceReader;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;

@Program(
		name="palindromefinder",
		description="naive palindrome finder",
		keywords={"palindrome","fasta"},
		creationDate = "20240424",
		modificationDate = "20240424",
		generate_doc = false
		)
public class PalindromeFinder  extends Launcher {
private static final Logger LOG = Logger.build(PalindromeFinder.class).make();

@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private Path outfile=null;
@Parameter(names={"-m"},description="min size")
private int min_palindrome_size = 10;


@Parameter(names={"-g"},description="min GC%." + FractionConverter.OPT_DESC ,splitter = NoSplitter.class,converter = FractionConverter.class)
private double min_gc=0.0;
@Parameter(names={"-G"},description="max GC%." + FractionConverter.OPT_DESC ,splitter = NoSplitter.class,converter = FractionConverter.class)
private double max_gc=1.0;
@Parameter(names={"-p"},description="disable overlaping palindrom (eg: regions with poly AT)")
private boolean disable_overlapping=false;


static boolean is_reverse_complement(char c1,char c2) {
	return Character.toUpperCase(AcidNucleics.complement(c1)) == Character.toUpperCase(c2);
	}



private void scan(final PrintWriter out,final FastaSequence seq) {
	int i=0;
	while(i< seq.length()) {
		int next_i=i+1;
		for(int skip_bases=0;skip_bases<1;skip_bases++) {
			for(int odd_even=0;odd_even<2;++odd_even) {
				final int pos = i;
				int best_len=0;
				int best_start=0;
				int best_end=0;
				int length=1;
				int mismatch=0;
				for(;;) {
					int px = pos-((skip_bases+length)- (odd_even==0?0:1));
					if(px<0) break;
					int py = pos+(skip_bases+length);
					if(py>=seq.length()) break;
					char cx = seq.charAt(px);
					char cy = seq.charAt(py);
					if(!is_reverse_complement(cx,cy)) {
						mismatch++;
						break;
						}
					if(best_len<length-mismatch) {
							best_len=length;
							best_start = px;
							best_end = py+1;
							}
					length++;
					}
				
				if(best_len>=this.min_palindrome_size) {
					float gc=0;
					for(int x=best_start;x< best_end;++x) {
						switch(seq.charAt(x)) {
							case 'c': case 'C':
							case 'G': case 'g':
							case 's': case 'S': gc++;break;
							default:break;
							}
						}
					double percent = gc/(best_end-best_start);
					if(this.min_gc<=percent && percent<= this.max_gc) {
						out.print(seq.getName());
						out.print('\t');
						out.print(best_start);
						out.print('\t');
						out.print(best_end);
						out.print('\t');
						out.print(best_end-best_start);
						out.print('\t');
						out.print(seq.subSequence(best_start, best_end));
						out.print('\t');
						out.print(percent);
						out.println();
						if(disable_overlapping) {
							next_i = Math.max(next_i, best_end);
							break;
							}
						}
					}
				}
			}
		i=next_i;
		}
	}

private void scan(PrintWriter out,InputStream in) throws IOException {
	try(CloseableIterator<FastaSequence> r=new FastaSequenceReader().iterator(in)) {
		while(r.hasNext()) {
			scan(out,r.next());
			}
		}
	}

@Override
public int doWork(List<String> args) {
	if(min_palindrome_size<4) {
		LOG.error("min size is too low");
		return -1;
		}
	try {
		final List<Path> paths = IOUtils.unrollPaths(args);
		try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outfile)) {
			if(paths.isEmpty()) {
				scan(out,System.in);
				}
			else
				{
				for(Path p:paths) {
					try(InputStream in = IOUtils.openPathForReading(p)) {
						scan(out,in);
						}
					}
				}
			out.flush();
			}
		return 0;
		}
	catch(Throwable err) {
		LOG.error(err);
		return -1;
		}
	}

public static void main(String[] args) {
	new PalindromeFinder().instanceMainWithExit(args);
	}
}
