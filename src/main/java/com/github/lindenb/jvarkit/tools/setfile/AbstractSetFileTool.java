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
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.UnaryOperator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.setfile.SetFileReaderFactory;
import com.github.lindenb.jvarkit.setfile.SetFileRecord;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;


/**
base class for SetFile tools

**/
public abstract class AbstractSetFileTool extends Launcher {
	private static final Logger LOG = Logger.of(AbstractSetFileTool.class);
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxRef = null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT+ ". For action=cluster, output is: "+ArchiveFactory.OPT_DESC)
	protected Path outputFile=null;
	@Parameter(names={"-t","--trim-chr"},description="Remove chr prefix in chromosome names on output.")
	protected boolean trim_chr_prefix = false;

	
	/** global fai */
	private SAMSequenceDictionary theDict = null;
	/** global sorter */
	private Comparator<Locatable> theSorter = null;
	/** global dict converter */
	private UnaryOperator<String> theCtgConverter=null;
	
	protected SAMSequenceDictionary getSAMSequenceDictionary() {
		if(theDict==null) {
			if(this.faidxRef!=null) {
				this.theDict = SequenceDictionaryUtils.extractRequired(this.faidxRef);
				}
			}
		return theDict;
		}
	
	protected Logger getLogger() {
		return LOG;
		}
	protected Comparator<Locatable> getSorter() {
		if(this.theSorter==null) {
			final SAMSequenceDictionary dict= getSAMSequenceDictionary();
			if(dict==null) {
				this.theSorter=(A,B)->{
					int i=A.getContig().compareTo(B.getContig());
					if(i!=0) return i;
					i = Integer.compare(A.getStart(), B.getStart());
					if(i!=0) return i;
					return Integer.compare(A.getEnd(), B.getEnd());
					};
				}
			else
				{
				this.theSorter = new ContigDictComparator(dict).createLocatableComparator();
				}
			}
		return this.theSorter;
		}
	
	protected UnaryOperator<String> getContigConverter() {
		if(this.theCtgConverter==null) {
			final SAMSequenceDictionary dict= getSAMSequenceDictionary();
			if(dict==null) {
				this.theCtgConverter=A->A;
				}
			else
				{
				this.theCtgConverter= ContigNameConverter.fromOneDictionary(dict);
				}
			}
		return this.theCtgConverter;
		}
		
	protected String noChr(final String contig) {
		if(trim_chr_prefix && contig.toLowerCase().startsWith("chr")) {
			return contig.substring(3);
			}
		return contig;
		}
	
	protected List<Locatable> sortAndMerge(final List<? extends Locatable> L0) {
			final List<Locatable> L = new ArrayList<>(L0);
			// merge overlapping
			Collections.sort(L, getSorter());
			int i=0;
			while(i+1 < L.size()) {
				final Locatable xi = L.get(i  );
				final Locatable xj = L.get(i+1);
				if(xi.overlaps(xj)) {
					L.set(i, new SimpleInterval(
							xi.getContig(),
							Math.min(xi.getStart(),xj.getStart()),
							Math.max(xi.getEnd(),xj.getEnd())
							));
					L.remove(i+1);
				} else {
					i++;
				}
			}
		return L;
		}
	
	
	/** print SetFileRecord to pw */
	protected void print(PrintWriter pw,SetFileRecord setfile) {
		if(setfile.isEmpty()) return;
		pw.write(setfile.getName());
		for(int i=0;i< setfile.size();i++) {
			final Locatable rec = setfile.get(i);
			pw.write(i==0?"\t":",");
			pw.write(noChr(rec.getContig()));
			pw.write(":");
			pw.write(String.valueOf(rec.getStart()));
			pw.write("-");
			pw.write(String.valueOf(rec.getEnd()));
			}
		pw.write("\n");
		}
	
	
	protected CloseableIterator<SetFileRecord> openSetFileIterator(final List<String> args) throws IOException {
		CloseableIterator<SetFileRecord> iter  = null;
		final String input = oneFileOrNull(args);
		final SetFileReaderFactory srf  = new SetFileReaderFactory(getSAMSequenceDictionary());
		if(input==null) {
			iter = srf.open(IOUtils.openStdinForBufferedReader());
		} else {
			iter  =srf.open(IOUtils.openURIForBufferedReading(input));
		}
		return iter;
		}
	
	protected int beforeJob() {
		return 0;
		}
	
	protected void afterJob() {
		}
	}
