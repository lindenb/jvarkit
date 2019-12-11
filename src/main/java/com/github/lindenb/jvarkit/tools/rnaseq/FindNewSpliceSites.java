/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

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
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

/**
BEGIN_DOC


## Example

```bash
$  java -jar dist/findnewsplicesites.jar \
     --gtf http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.gtf.gz \
      hg19.bam > out.sam
```

END_DOC
*/
@Program(name="findnewsplicesites",
	description="use the 'N' operator in the cigar string to find unknown splice sites",
	keywords={"rnaseq","splice"},
	modificationDate="20190909"
	)
public class FindNewSpliceSites extends Launcher
	{
	private static final Logger LOG = Logger.build(FindNewSpliceSites.class).make();

	
	@Parameter(names={"-out","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-B","--bed"},description="Optional BED output")
	private Path bedOut = null;
	@Parameter(names={"-g","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfPath = null;
	@Parameter(names="-d",description="max distance between known splice site and cigar end")
	private int max_distance=10;
	@Parameter(names= {"-R","--reference"},description="For reading cram. "+ INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	
	private SAMFileWriter sfw=null;
	private SAMFileWriter weird=null;
	private final IntervalTreeMap<Interval> intronMap =new IntervalTreeMap<>();
	private PrintWriter bedWriter = null;
	
	private FindNewSpliceSites()
		{
		}
	
	/** distance known splice site and cigar */
	private boolean is_close_to(int d1,int d2)
		{
		return Math.abs(d1-d2)<=this.max_distance;
		}
	


	private boolean findJunction(
			final Collection<Interval> introns,
			final int start1,
			final int end1
			)
		{
		return introns.
			stream().
			anyMatch(
			intron-> 
				is_close_to(intron.getStart(),start1) &&
				is_close_to(intron.getEnd(),end1
				));
		}
	private void scanRead(
			final SAMRecord rec,
			final SAMSequenceDictionary dict,
			final SAMProgramRecord okPrg
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
						cigar.getCigarElement(cIdx+1).getOperator().isAlignment() &&	
						!findJunction(introns,refPos1,refPos1+ce.getLength()-1))
						{
						if(!readSaved) {
							rec.setAttribute("PG", okPrg.getId());
							this.sfw.addAlignment(rec);//unknown junction
							readSaved=true;
							}
						if(this.bedOut==null)
							{
							return;
							}
						else
							{
							this.bedWriter.print(rec.getContig());
							this.bedWriter.print('\t');
							this.bedWriter.print(refPos1-1);
							this.bedWriter.print('\t');
							this.bedWriter.print(refPos1+ce.getLength()-1);
							this.bedWriter.print('\t');
							this.bedWriter.print(rec.getReadName());
							this.bedWriter.print('\t');
							this.bedWriter.print(rec.getReadNegativeStrandFlag()?'-':'+');
							this.bedWriter.println();
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
			final SAMProgramRecord failPrg) 
		{
		final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(in.getFileHeader());
		if(dict==null) throw new RuntimeException("Sequence dictionary missing");
		try(SAMRecordIterator iter=in.iterator() ) {
			ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
			
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
					
					scanRead(rec,dict,okPrg);
					}	
				}
			progress.close();
			}
		
		}
	@Override
	public int doWork(final List<String> args) {
		SamReader sfr=null;
		try
			{
			try(GtfReader gftReader=new GtfReader(this.gtfPath))
				{
				gftReader.getAllGenes().stream().
					flatMap(G->G.getTranscripts().stream()).
					filter(T->T.getExonCount()>1).
					flatMap(T->T.getIntrons().stream()).
					map(T->T.toInterval()).
					forEach(T->
					{
					this.intronMap.put(T,T);
					});
				}
			
			this.bedWriter = this.bedOut==null?
					new PrintWriter(new NullOuputStream()):
					super.openPathOrStdoutAsPrintWriter(this.bedOut);
			
			final SamReaderFactory srf = super.createSamReaderFactory();		
			if(this.faidx!=null) {
				srf.referenceSequence(this.faidx);
				}
			
			final String input = oneFileOrNull(args);
	
			sfr = input==null?srf.open(SamInputResource.of(stdin())):
				srf.open(SamInputResource.of(input));
			
			SAMFileHeader header=sfr.getFileHeader().clone();
			final SAMProgramRecord p=header.createProgramRecord();
			p.setCommandLine(getProgramCommandLine());
			p.setProgramVersion(getVersion());
			p.setProgramName(getProgramName());
			this.sfw=this.writingBamArgs.openSamWriter(outputFile, header, true);
			
			header=sfr.getFileHeader().clone();
			final SAMProgramRecord p2=header.createProgramRecord();
			p2.setCommandLine(getProgramCommandLine());
			p2.setProgramVersion(getVersion());
			p2.setProgramName(getProgramName());
			this.weird=this.writingBamArgs.createSAMFileWriterFactory().makeSAMWriter(header,true, new NullOuputStream());
			
			scan(sfr,p,p2);
			sfr.close();
			
			this.bedWriter.flush();
			this.bedWriter.close();
			this.bedWriter=null;
			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
			CloserUtil.close(this.sfw);
			CloserUtil.close(this.weird);
			CloserUtil.close(this.bedWriter);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindNewSpliceSites().instanceMainWithExit(args);

	}

}
