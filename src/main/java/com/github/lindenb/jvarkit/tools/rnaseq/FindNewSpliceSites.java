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
package com.github.lindenb.jvarkit.tools.rnaseq;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SortingCollection;

/**
BEGIN_DOC


## Example

```bash
$  java -jar dist/jvarkit.git findnewsplicesites  \
     --gtf /pah/to/file.gtf.gz \
      my.bam > out.sam
```

END_DOC
*/
@Program(name="findnewsplicesites",
	description="use the 'N' operator in the cigar string to find unknown splice sites",
	keywords={"bam","sam","rnaseq","splice","gtf"},
	creationDate="20140402",
	modificationDate="20230608",
	jvarkit_amalgamion = true,
	menu="BAM Manipulation"
	)
public class FindNewSpliceSites extends Launcher
	{
	private static final Logger LOG = Logger.of(FindNewSpliceSites.class);
	private static long GENERATE_ID=0L;
	
	@Parameter(names={"-o","-out","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-B","--bed"},description="Optional BED output")
	private Path bedOut = null;
	@Parameter(names={"-g","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfPath = null;
	@Parameter(names="-d",description="max distance between known splice site and cigar end")
	private int max_distance=0;
	@Parameter(names= {"-R","--reference"},description="For reading cram. "+ INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx;
	@Parameter(names= {"--debug"},description="debug",hidden = true)
	private boolean do_debug=false;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection= new WritingSortingCollection();
	
	
	private SAMFileWriter sfw=null;
	private SAMFileWriter weird=null;
	private final IntervalTreeMap<Interval> intronMap =new IntervalTreeMap<>();
	private Comparator<Junction> junctionComparator = null;
	
	private class Junction implements Locatable {
		long id=GENERATE_ID++;
		String contig;
		int start;
		int end;
		String name;
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getStart() {
			return start;
			}
		@Override
		public int getEnd() {
			return end;
			}
		public int compare1(final Junction other) {
			return junctionComparator.compare(this, other);
			}
		public int compare2(final Junction other) {
			final int i= compare1(other);
			if(i!=0) return i;
			return Long.compare(this.id,other.id);
			}
		}
	
	private class JunctionCodec extends AbstractDataCodec<Junction> {
		@Override
		public Junction decode(DataInputStream dis) throws IOException {
			final Junction j = new Junction();
			try {
				j.id=dis.readLong();
			} catch(EOFException err) {
				return null;
				}
			j.contig = dis.readUTF();
			j.start=dis.readInt();
			j.end=dis.readInt();
			j.name=dis.readUTF();
			return j;
			}
		@Override
		public void encode(DataOutputStream dos, Junction j) throws IOException {
			dos.writeLong(j.id);
			dos.writeUTF(j.contig);
			dos.writeInt(j.start);
			dos.writeInt(j.end);
			dos.writeUTF(j.name);

			}
		@Override
		public AbstractDataCodec<Junction> clone() {
			return new JunctionCodec();
			}
		
		}
		
	/** distance known splice site and cigar */
	private boolean is_close_to(int d1,int d2)
		{
		return Math.abs(d1-d2)<=this.max_distance;
		}
	


	private String findJunctionName(
			final Collection<Interval> introns,
			final int start1,
			final int end1
			)
		{
		return introns.
			stream().
			filter(
			intron-> 
				is_close_to(intron.getStart(),start1) &&
				is_close_to(intron.getEnd(),end1
				)).map(intron->intron.getName()).
			findFirst().
			orElse(null);
		}
	
	private void scanRead(
			final SAMRecord rec,
			final SAMSequenceDictionary dict,
			final SAMProgramRecord okPrg,
			final SortingCollection<Junction> sortingJunctions
			)
			{
		final Cigar cigar=rec.getCigar();
		
		final Collection<Interval> introns =this.intronMap.getOverlapping(rec);
		
		if(introns.isEmpty())
			{
			return;
			}

		boolean readSaved=false;
		int refPos1=rec.getAlignmentStart();
		
		
		for(int cIdx=0;cIdx< cigar.numCigarElements();++cIdx)
			{
			final CigarElement ce=cigar.getCigarElement(cIdx);
			switch(ce.getOperator())
				{
				case S: break;
				case I: break;
				case N:
					{
					if(cIdx+1<cigar.numCigarElements() &&
						cigar.getCigarElement(cIdx+1).getOperator().isAlignment())
						{
						final String  junctionName = findJunctionName(introns,refPos1,refPos1+ce.getLength()-1);
						
						
						if(junctionName==null && !readSaved) {
							rec.setAttribute("PG", okPrg.getId());
							this.sfw.addAlignment(rec);//unknown junction
							readSaved=true;
							}
						if(sortingJunctions==null)
							{
							return;
							}
						else
							{
							final Junction j = new Junction();
							j.contig = rec.getContig();
							j.start = refPos1;
							j.end = refPos1+ce.getLength()-1;
							j.name  = StringUtils.ifBlank(junctionName, ".");
							
							sortingJunctions.add(j);
							
							}
						}
					refPos1+=ce.getLength();	
					break;
					}
				case D:
				case X:
				case EQ:
				case M:
					{
					refPos1+=ce.getLength();
					break;
					}
				case H:case P: break;//ignore
				default:throw new IllegalStateException("operator not handled. ops.");
				}
			}
		}

	
	
	private static boolean isWeird(final SAMRecord rec,final SAMSequenceDictionary dict)
		{
		if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag() && 
				rec.getReferenceIndex().equals(rec.getMateReferenceIndex()) &&
				(
				rec.getReadNegativeStrandFlag()==rec.getMateNegativeStrandFlag() ||
				(rec.getReadNegativeStrandFlag()&& !rec.getMateNegativeStrandFlag() && rec.getAlignmentStart() < rec.getMateAlignmentStart()) ||
				(!rec.getReadNegativeStrandFlag() && rec.getMateNegativeStrandFlag() && rec.getAlignmentStart() > rec.getMateAlignmentStart())
				))
			{
			if(rec.getAlignmentStart() < rec.getMateAlignmentStart())
				{
				if(rec.getMateAlignmentStart() < rec.getAlignmentEnd()) return true;
				}
			if(rec.getAlignmentStart() > rec.getMateAlignmentStart())
				{
				if(rec.getMateAlignmentStart() + Math.abs(rec.getInferredInsertSize() ) > rec.getAlignmentStart()) return true;
				}
			return true;
			}
		return false;
		}
	
	private void scan(final SamReader in,
			final SAMProgramRecord okPrg,
			final SAMProgramRecord failPrg,
			final SortingCollection<Junction> junctionSorter
			) 
		{
		final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(in.getFileHeader());
		
		
		
		try(SAMRecordIterator iter=in.iterator() ) {
			final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
			
			while(iter.hasNext())
				{
				final SAMRecord rec=progress.apply(iter.next());
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				
				if(isWeird(rec,dict))
					{
					rec.setAttribute("PG", failPrg.getId());
					this.weird.addAlignment(rec);
					continue;
					}
				
				if(rec.getCigar().getCigarElements().
						stream().
						anyMatch(C->CigarOperator.N.equals(C.getOperator()))) {
					
					scanRead(rec,dict,okPrg,junctionSorter);
					}	
				}
			progress.close();
			}
		
		}
	@Override
	public int doWork(final List<String> args) {
		SortingCollection<Junction> junctionSorter = null;
		try
			{
			
			
			final SamReaderFactory srf = super.createSamReaderFactory();		
			if(this.faidx!=null) {
				srf.referenceSequence(this.faidx);
				}
			
			final String input = oneFileOrNull(args);
	
			final SAMFileHeader header0;
			
			try(SamReader sfr = input==null?srf.open(SamInputResource.of(stdin())):
					srf.open(SamInputResource.of(input))) {
				
				header0 =  sfr.getFileHeader();
				
				
				try(BufferedReader br = IOUtils.openPathForBufferedReading(this.gtfPath))
					{
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header0);
					final Comparator<Locatable> orderer = new ContigDictComparator(dict).createLocatableComparator();
					final GTFCodec codec = new GTFCodec();
					codec.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					final Map<String,List<Interval>> tx2exons = new HashMap<>(50_000);
					for(;;) {
						final String line = br.readLine();
						if(line==null) break;
						final GTFLine rec = codec.decode(line);
						if(rec==null) continue;
						if(!rec.isExon()) continue;
						if(!rec.hasAttribute("transcript_id")) {
							LOG.warn("no transcript_id in "+line);
							continue;
							}
						final String tr = rec.getAttribute("transcript_id");
						List<Interval> L = tx2exons.get(tr);
						if(L==null) {
							L=new ArrayList<>();
							tx2exons.put(tr, L);
							}
						L.add(new Interval(rec.getContig(),rec.getStart(),rec.getEnd(),rec.isNegativeStrand(),tr));
						}
					for(final String transcript_id: tx2exons.keySet() ) {
						final List<Interval> exons = tx2exons.get(transcript_id);
						if(exons.size()<2) continue;//no intron
						Collections.sort(exons,orderer);
						for(int i=0;i+1< exons.size();i++) {
							if(!exons.get(i).contigsMatch(exons.get(i+1))) {
								LOG.warn("multiple contigs in "+exons);
								continue;
								}
							if(exons.get(i).overlaps(exons.get(i+1))) {
								LOG.warn("overlapping exons: "+exons.get(i)+" and "+exons.get(i+1));
								continue;
								}
							if(!exons.get(i).getStrand().equals(exons.get(i+1).getStrand())) {
								LOG.warn("multiple strand exons: "+exons.get(i)+" and "+exons.get(i+1));
								continue;
								}
							
							final Interval intron = new Interval(
									exons.get(i).getContig(),
									exons.get(i  ).getEnd()+1,
									exons.get(i+1).getStart()-1,
									exons.get(0).isNegativeStrand(),
									exons.get(0).getName()+".Intron"+(
											exons.get(0).isPositiveStrand()?
											i+1:
											exons.size()-(i+1))
									);
							this.intronMap.put(intron,intron);
							}
						}
					if(do_debug) {
						for(Interval intron:this.intronMap.keySet()) {
							LOG.debug("found intron :"+intron);
							}
						}
					}
				
				final SAMFileHeader header1= header0.clone();
				final SAMProgramRecord p=header1.createProgramRecord();
				p.setCommandLine(getProgramCommandLine());
				p.setProgramVersion(getVersion());
				p.setProgramName(getProgramName());
				this.sfw=this.writingBamArgs.openSamWriter(outputFile, header1, true);
				
				final SAMFileHeader header2 = header0.clone();
				final SAMProgramRecord p2=header2.createProgramRecord();
				p2.setCommandLine(getProgramCommandLine());
				p2.setProgramVersion(getVersion());
				p2.setProgramName(getProgramName());
				this.weird=this.writingBamArgs.createSAMFileWriterFactory().makeSAMWriter(header2,true, new NullOuputStream());
				
				if(this.bedOut!=null) {
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(sfr.getFileHeader());
					this.junctionComparator = new ContigDictComparator(dict).createLocatableComparator();
					junctionSorter = SortingCollection.newInstance(
							Junction.class,
							new JunctionCodec(),
							(A,B)->A.compare2(B),
							this.writingSortingCollection.getMaxRecordsInRam(),
							this.writingSortingCollection.getTmpPaths()
							);
					}
				
				
				scan(sfr,p,p2,junctionSorter);
				}
			
			if(this.bedOut!=null) {
				junctionSorter.doneAdding();
				
				try(PrintWriter bedWriter =  super.openPathOrStdoutAsPrintWriter(this.bedOut)) {
		
					final String sample = StringUtils.ifBlank(
							header0.
							getReadGroups().
							stream().
							map(RG->RG.getSample()).
							filter(s->!StringUtils.isBlank(s)).
							collect(Collectors.toCollection(TreeSet::new)).
							stream().
							collect(Collectors.joining(";")),
							"."
							);
						
					
					try(CloseableIterator<Junction> iter=junctionSorter.iterator()) {
						final EqualRangeIterator<Junction> eq = new EqualRangeIterator<>(iter, (A,B)->A.compare1(B));
						while(eq.hasNext()) {
							final List<Junction> row = eq.next();
							final Junction first = row.get(0);
							bedWriter.print(first.getContig());
							bedWriter.print('\t');
							bedWriter.print(first.getStart()-1);
							bedWriter.print('\t');
							bedWriter.print(first.getEnd());
							bedWriter.print('\t');
							bedWriter.print(sample);
							bedWriter.print('\t');
							bedWriter.print(first.name);
							bedWriter.print('\t');
							bedWriter.print(row.size());
							bedWriter.println();
							}
						eq.close();
						}
					
					
					bedWriter.flush();
					}
				junctionSorter.cleanup();
				}
				
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(this.sfw!=null) try {this.sfw.close();} catch(Throwable e) {};
			if(this.weird!=null) try {this.weird.close();} catch(Throwable e) {};
			}
		}
	
	public static void main(String[] args) {
		new FindNewSpliceSites().instanceMainWithExit(args);

	}

}
