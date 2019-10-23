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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;


/**
BEGIN_DOC

##Example

```

$ samtools view src/test/resources/toy.bam | tail -1
x6	0	ref2	14	30	23M	*	0	0	TAATTAAGTCTACAGAGCAACTA	???????????????????????	RG:Z:gid1

$ cat jeter.vcf
ref2	14	A	C

$ java -jar dist/biostar404363.jar -p jeter.vcf src/test/resources/toy.bam | tail -1
x6	0	ref2	14	30	23M	*	0	0	CAATTAAGTCTACAGAGCAACTA	???????????????????????	PG:Z:0	RG:Z:gid1	NM:i:1



```

END_DOC
*/
@Program(name="biostar404363",
	keywords={"sam","bam","variant"},
	description="introduce artificial mutation in bam",
	biostars=404363
	)
public class Biostar404363 extends Launcher {
	
	private static final Logger LOG = Logger.build(Biostar404363.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--position","--vcf"},description="File containing the positions to change: syntax (looks like a VCF line): 'CHROM\\tPOS\\t(ignored)\\tBASE",required=true)
	private Path vcfPath=null;

	@Parameter(names= {"-R","--reference"},description="For reading CRAM. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	private static class PosToChange extends SimplePosition {
		byte base = 'N';
		PosToChange(final String ctg,int pos) {
			super(ctg,pos);
		}
	}

	@Override
	public int doWork(final List<String> args) {
		SamReader in=null;
		SAMFileWriter out=null;
		final IntervalTreeMap<PosToChange> intervalTreeMap = new IntervalTreeMap<>();
		try {
			
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.faidx!=null) {
				srf.referenceSequence(this.faidx);
				writingBamArgs.setReferencePath(this.faidx);
				}
			final String input = oneFileOrNull(args);
			if(input==null) {
				in = srf.open(SamInputResource.of(stdin()));
				}
			else
				{
				in = srf.open(SamInputResource.of(input));
				}
			final SAMFileHeader header = in.getFileHeader();
			final SAMProgramRecord prg =  header.createProgramRecord();
			prg.setCommandLine(this.getProgramCommandLine());
			prg.setProgramName(this.getProgramName());
			prg.setProgramVersion(this.getGitHash());
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			
			BufferedReader br=IOUtils.openPathForBufferedReading(this.vcfPath);
			String line;
			while((line=br.readLine())!=null) {
				if(line.startsWith("#")) continue;
				final String tokens[] = CharSplitter.TAB.split(line);
				if(tokens.length<4 ||
						!StringUtils.isInteger(tokens[1]) || 
						tokens[3].length()!=1 || 
						!AcidNucleics.isATGCN(tokens[3].charAt(0))) {
						LOG.error("Syntax error "+line);
						return -1;
						}
				final String theContig = tokens[0];
				if(dict.getSequence(theContig)==null) throw new JvarkitException.ContigNotFoundInDictionary(theContig, dict);

				final int thePos = Integer.parseInt(tokens[1]);
				final PosToChange pos2change = new PosToChange(theContig, thePos);
				pos2change.base=(byte)Character.toUpperCase(tokens[3].charAt(0));
				
				intervalTreeMap.put(new Interval(pos2change),pos2change);
				}
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			out = this.writingBamArgs.openSamWriter(this.outputFile, header, true);
			final SAMRecordIterator iter=in.iterator();
			while(iter.hasNext()) {
				final SAMRecord rec= iter.next();
				final byte bases[]=rec.getReadBases();

				if(rec.getReadUnmappedFlag() ||
					bases==SAMRecord.NULL_SEQUENCE) {
					out.addAlignment(rec);
					continue;
					}
				
				final List<PosToChange> changes = intervalTreeMap.getOverlapping(new SimpleInterval(rec.getContig(), rec.getUnclippedStart(), rec.getUnclippedEnd())).
						stream().
						collect(Collectors.toList());
				if(changes.isEmpty()) {
					out.addAlignment(rec);
					continue;
					}
				
				
				int NM=0;
				if(rec.hasAttribute("NM")) NM=rec.getIntegerAttribute("NM");
				int readpos=0;
				int ref=rec.getUnclippedStart();
				boolean changed=false;
				final Cigar cigar = rec.getCigar();
				for(final CigarElement ce:cigar) {
					final CigarOperator op =ce.getOperator();
					
					if(op.equals(CigarOperator.S) || (op.consumesReadBases() &&  op.consumesReferenceBases() )) {
						{
						for(int i=0;i< ce.getLength();i++) {
							final int the_pos = ref+i;
							final byte the_base= bases[readpos+i]; 
							final PosToChange pos2change = changes.stream().
									filter(P->P.getPosition()==the_pos).
									filter(P->P.base!=the_base).
									findFirst().
									orElse(null);
							if(pos2change==null) continue;
							
							bases[readpos+i]=pos2change.base;
							NM++;
							changed=true;
							}
						}
					if(op.equals(CigarOperator.S) || op.consumesReadBases())
						readpos+=ce.getLength();
						}
					if(op.isClipping() || op.consumesReferenceBases()) {
						ref+=ce.getLength();
						}
					
					}
				
				if(changed) {
					rec.setReadBases(bases);
					rec.setAttribute("NM", NM);
					rec.setAttribute("PG", prg.getId());
					}
				out.addAlignment(rec);
				}
			in.close();
			in=null;
			out.close();
			out=null;
			
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(in);
			CloserUtil.close(out);
		}
	}
	
	public static void main(final String[] args) {
		new Biostar404363().instanceMainWithExit(args);
	}

}
