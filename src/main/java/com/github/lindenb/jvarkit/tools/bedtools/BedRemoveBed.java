/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bedtools;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.BitSet;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
/**
BEGIN_DOC

## Motivation

remove bed from each record of bed

## EXAMPLES

### Example


END_DOC
*/
@Program(
		name="bedremovebed",
		description="Remove bed file from each record of input bed file. Output is a SETFILE",
		keywords={"bed"},
		creationDate="20221210",
		modificationDate="20221210",
		jvarkit_amalgamion = true
		)
public class BedRemoveBed extends Launcher {
	private static final Logger LOG = Logger.of(BedRemoveBed.class);

	@Parameter(names = { "-o", "--out" }, description = "Output is a setfile. "+ OPT_OUPUT_FILE_OR_STDOUT)
	private Path output = null;
	@Parameter(names = { "-exclude","-X","--exclude-bed" }, description = "Bed to be substracted from original input (gene)")
	private Path excludeBed = null;
	@Parameter(names = { "-input","-I","--include-bed" }, description = "Bed to be keep from original input (gene)")
	private Path includeBed = null;
	@Parameter(names = { "-R", "-r","--reference" }, description = DICTIONARY_SOURCE )
	private Path refFile = null;
	
	
	@Override
	public int doWork(final List<String> args) {
				
		try {
			final String input = oneFileOrNull(args);
			final SAMSequenceDictionary dict;
			if(this.refFile!=null)
				{
				dict = SAMSequenceDictionaryExtractor.extractDictionary(this.refFile);
				}
			else
				{
				dict = null;
				}
			final IntervalTreeMap<Locatable> excludeTree;
			final IntervalTreeMap<Locatable> includeTree;

			if(excludeBed!=null && includeBed==null) {
				excludeTree = new IntervalTreeMap<>();
				includeTree = null;
				try(BedLineReader blr=  new BedLineReader(excludeBed)) {
						if(dict!=null) {
							blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
							}
						while(blr.hasNext()) {
							final Interval r = new Interval(blr.next());
							excludeTree.put(r, r);
							}
						}
					}
			else if(includeBed!=null && excludeBed==null) {
					includeTree = new IntervalTreeMap<>();
					excludeTree = null;
					try(BedLineReader blr=  new BedLineReader(includeBed)) {
						if(dict!=null) {
							blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
							}
						while(blr.hasNext()) {
							final Interval r = new Interval(blr.next());
							includeTree.put(r, r);
							}
						}
				}
			else
				{
				LOG.error("both or none option --exclude and --include were specified.");
				return -1;
				}
			
			try(BedLineReader blr= StringUtils.isBlank(input)?
					new BedLineReader(stdin(), "STDIN"):
					new BedLineReader(Paths.get(input))) {
				if(dict!=null) {
					blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					}
				long file_id=0L;
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.output)) {
					while(blr.hasNext()) {
						final BedLine rec = blr.next();
						final BitSet bitSet = new BitSet(rec.getLengthOnReference());
						bitSet.set(0, bitSet.length(), excludeTree!=null);
						
						for(Locatable b:(excludeTree!=null?excludeTree:includeTree).getOverlapping(rec)) {
							final int x0 = Math.max(rec.getStart(), b.getStart());
							final int x1 = Math.min(rec.getEnd(), b.getEnd());
							bitSet.set(
								x0-rec.getStart(),
								x1-rec.getStart(),
								excludeTree==null
								);
							}
						int p = bitSet.nextSetBit(0);
						if(p==-1) continue;
						++file_id;
						pw.print(rec.getContig()+"_"+rec.getStart()+"_"+rec.getEnd()+"_"+String.format("%09d", file_id));
						int nfrags=0;
						while(p!=-1) {
							int p2 = bitSet.nextClearBit(p+1);
							if(p2==-1) p2=bitSet.length();
							pw.print(nfrags==0?"\t":",");
							pw.print(rec.getContig());
							pw.print("\t");
							pw.print(rec.getStart()+p);
							pw.print("\t");
							pw.print(rec.getStart()+p2-1);
							nfrags++;
							if(p2>=bitSet.length()) break;
							p = bitSet.nextSetBit(p2+1);
							}
						pw.println();
						} //end of p
					pw.flush();
					}// end of printwriter
				} // end of reader
		
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args) {
		new BedRemoveBed().instanceMainWithExit(args);
	}

}
