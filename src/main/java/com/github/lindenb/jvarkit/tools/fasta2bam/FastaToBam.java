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
package com.github.lindenb.jvarkit.tools.fasta2bam;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

public class FastaToBam extends Launcher {
	private static final Logger LOG = Logger.of(FastaToBam.class);

	@Parameter(names={"--instrument"},description="instrument")
	private String instrument="SEQ1";
	@Parameter(names={"--flowcell"},description="flowcell")
	private String flowcell="FC706VJ";
	@Parameter(names={"--run"},description="run id")
	private int run_id = 1;


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path fileout = null;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	
	public FastaToBam() {
		}
	
	private static class Pos2Base {
		int pos1;
		char base;
		boolean het=false;
		}
	
	private void flushReads(
				final SAMFileWriter w,
				final Integer pos1,
				final List<SAMRecord> buffer,
				final SAMRecordCoordinateComparator comparator
				) {
			List<SAMRecord> sort=null;
			if(pos1==null) {
				sort = new ArrayList<>(buffer);
				buffer.clear();
				}
			else
				{
				sort = null;
				int i=0;
				while(i< buffer.size()) {
					final SAMRecord rec = buffer.get(i);
					if(rec.getStart() < pos1   ) {
						if(sort==null) sort = new ArrayList<>();
						sort.add(rec);
						buffer.remove(i);
						}
					else
						{
						i++;
						}
					}
				if(sort==null) return;
				}
			Collections.sort(sort, comparator);
			for(SAMRecord rec: sort) {
				w.addAlignment(rec);
				}
			}
	
	private char baseAt(GenomicSequence ref, int pos1, List<Pos2Base> pos2base) {
		char base='x';
		while(pos2base.isEmpty() || pos2base.get(pos2base.size()-1).pos1 < pos1) {
			
			}
		return base;
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			final Random rnd = new Random();
			Path fasta = Paths.get( oneAndOnlyOneFile(args) ) ;
			try(ReferenceSequenceFile ref= ReferenceSequenceFileFactory.getReferenceSequenceFile(fasta)) {
				SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(ref);
				SAMFileHeader header = new SAMFileHeader(dict);
				final SAMProgramRecord prg = header.createProgramRecord();
				prg.setCommandLine(super.getProgramCommandLine());
				prg.setProgramName(getProgramName());
				prg.setProgramVersion(getGitHash());
				JVarkitVersion.getInstance().addMetaData(this, header);
				header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
				
				List<SAMReadGroupRecord> readGroups = new ArrayList<>();
				header.setReadGroups(readGroups);

				final DefaultSAMRecordFactory samRecordFactory = new DefaultSAMRecordFactory();
				
				final SAMRecordCoordinateComparator samRecordCoordinateComparator = new  SAMRecordCoordinateComparator();
				
				try(SAMFileWriter w= writingBamArgs.setReferencePath(fasta).openSamWriter(this.fileout, header, true)) {
					
					
					for(int tid=0; tid <dict.size();tid++) {
						final SAMSequenceRecord ssr = dict.getSequence(tid);
						final GenomicSequence genomicSequence = new GenomicSequence(ref, ssr.getContig());
						int read_length=150;
						double coverage = 10.0;
						double proba = 1.0/(read_length*coverage);
						final List<SAMRecord> buffer=new ArrayList<>();
						final List<Pos2Base> pos2base = new ArrayList<>();
						for(int x=1; x+read_length <= ssr.getSequenceLength();++x) {
							final int x1=x;
							final char ref_base = genomicSequence.charAt(x1-1);
							if(!AcidNucleics.isATGC(ref_base)) continue;
							
							pos2base.removeIf(P2B -> P2B.pos1 < x1);
							
							
							
							flushReads(w,x,buffer, samRecordCoordinateComparator);
							
							while(rnd.nextDouble() < proba) {
								final String readName = String.join(
										":",
										instrument,
										String.valueOf(run_id),
										String.valueOf(1+ rnd.nextInt(10)), //lane
										String.valueOf(1+ rnd.nextInt(1000)), //tile
										String.valueOf(1+ rnd.nextInt(1000)), //X
										String.valueOf(1+ rnd.nextInt(1000)) //Y
										);
								List<SAMRecord> pairs = new ArrayList<>();
								for(int side=0;side < 2;++side) {
									final SAMRecord rec = samRecordFactory.createSAMRecord(header);
									rec.setReadName(readName);
									StringBuilder seq  = new StringBuilder(read_length);
									StringBuilder qual = new StringBuilder(read_length);
									List<CigarOperator> ops= new ArrayList<>();
									
									for(int j=0;j< read_length;++j) {
										seq.append( baseAt(genomicSequence,x1,pos2base));
										qual.append('2');
										}
									rec.setReadString(seq.toString());
									rec.setBaseQualityString(qual.toString());
									
									
									rec.setAttribute("PG", prg.getId());
									rec.setCigar(Cigar.fromCigarOperators(ops));

									
									buffer.add(rec);
									}
								if(pairs.size()==1) {
									
									}
								else
									{
									pairs.get(0).setReadPairedFlag(true);
									pairs.get(1).setReadPairedFlag(true);
									
									
									Collections.shuffle(pairs, rnd);
									pairs.get(0).setFirstOfPairFlag(true);
									pairs.get(1).setSecondOfPairFlag(true);
									}
								buffer.addAll(pairs);
								}
							}
						flushReads(w,null,buffer, samRecordCoordinateComparator);
						}
					}
				}
			
			return 0;
			}
		catch(Throwable err) {
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new FastaToBam().instanceMainWithExit(args);
	}
}
