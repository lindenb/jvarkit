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
package com.github.lindenb.jvarkit.tools.textbam;


import java.nio.file.Path;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
/**
BEGIN_DOC

```
$ ./gradlew textbam && \
	java -jar dist/textbam.jar -R src/test/resources/rotavirus_rf.fa -p 'RF01:100' "Hello world" | samtools view -O BAM -o jeter.bam &&  \
	samtools index jeter.bam && \
	samtools tview ~/jeter.bam src/test/resources/rotavirus_rf.fa


samtools tview ~/jeter.bam src/test/resources/rotavirus_rf.fa -d T -p RF01:100
 101       111       121       131       141       151       161                
TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATC
................................................................................
................................................................................
 ...............................................................................
  ..............*....*......******......*...........*............*****..........
   .............*....*......*...........*...........*............*....*.........
    ............*....*......*...........*...........*...........*.....*.........
     ...........******......****........*...........*...........*.....*.........
      ..........*....*......*...........*...........*...........*.....*.........
       .........*....*......*...........*...........*...........*....*..........
        ........*....*......*...........*...........*............**..*..........
         .......*....*......******......******......******.........**...........
          ......................................................................
           .....................................................................
```

END_DOC
**/


@Program(name="texbam",
description="Write text in a bam. Mostly for fun...",
keywords={"fun","bam","sam","txt"},
creationDate="20220708",
modificationDate="20220708",
jvarkit_amalgamion = true,
menu="BAM Manipulation"
)
public class TextBam extends Launcher {
	private static final Logger LOG = Logger.build( TextBam.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path refPath = null;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--snp"},description="use SNV instead of deletion")
	private boolean use_snv_instead_of_del = false;
	@Parameter(names={"-p","--pos"},description="use this position instead of a random one. Syntax: CHROM:POS")
	private String userPos="";


	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	@Override
	public int doWork(final List<String> args) {
		try {
			final Hershey hershey = new Hershey();
			final String input = String.join(" ", args);
			if(StringUtils.isBlank(input)) {
				LOG.error("no text or text is blank");
				return -1;
				}
			int size=12;
			final List<BitSet> matrix = hershey.asciiArt(input, input.length()*size, size);
			for(BitSet row:matrix) {
				for(int x=0;x< row.length();x++) {
					System.err.print(row.get(x)?'#':' ');
					}
				System.err.println();
				}
			if(matrix.isEmpty()) {
				LOG.info("empty matrix");
				return -1;
				}
			final int matrixWidth = matrix.get(0).size();
			System.err.println("Matrix width: "+matrixWidth);
			final Random rnd = new Random();
			try(ReferenceSequenceFile fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.refPath)) {
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(fasta);
				final SAMSequenceRecord ssr;
				final int startpos;
				if(!StringUtils.isBlank(userPos)) {
					final SimplePosition sp = new SimplePosition(userPos);
					ssr = dict.getSequence(sp.getContig());
					if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(sp.getContig(), dict);
					startpos = sp.getPosition() ;
					if(startpos-matrix.size() <0 || startpos+matrixWidth>=ssr.getLengthOnReference()) {
						LOG.error("Cannot use "+sp+" string would be out of bounds");
						return -1;
						}
					}
				else
					{
					final List<SAMSequenceRecord>ssrList = dict.getSequences().
							stream().
							filter(SSR->SSR.getLengthOnReference() > matrixWidth+matrix.size()).
							collect(Collectors.toList());
					if(ssrList.isEmpty()) {
						LOG.info("no suitable chromosome");
						return -1;
						}
					ssr = ssrList.get(rnd.nextInt(ssrList.size()));
					startpos = rnd.nextInt(ssr.getLengthOnReference()-matrixWidth);
					}
				final SAMRecordFactory srf = new DefaultSAMRecordFactory();
				final GenomicSequence genomic = new GenomicSequence(fasta, ssr.getContig());
				final SAMFileHeader header = new SAMFileHeader(dict);
				final SAMReadGroupRecord rg = new SAMReadGroupRecord("S1");
				rg.setSample("sample1");
				rg.setLibrary("sample1");
				header.addReadGroup(rg);
				header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
				JVarkitVersion.getInstance().addMetaData(this, header);

				try(SAMFileWriter w = this.writingBamArgs.openSamWriter(this.outputFile, header, true)) {
					for(int row_index=0;row_index < matrix.size();row_index++) {
						final BitSet row = matrix.get(row_index);
						int nm=0;
						final int read_start = startpos + row_index;
						final SAMRecord rec = srf.createSAMRecord(header);
						rec.setReadName("HWI-ST975:118:C12H5ACXX:1:"+String.format("%03d", row_index)+":"+rnd.nextInt(1000)+":"+rnd.nextInt(1000));
						rec.setReferenceIndex(ssr.getSequenceIndex());
						rec.setAlignmentStart(Math.max(1, read_start));
						
						final StringBuilder basesBuilder = new StringBuilder(matrixWidth);
						final StringBuilder qualBuilder = new StringBuilder(matrixWidth);
						final List<CigarOperator> operators = new ArrayList<>();
						for(int i=-1 ; i< matrixWidth+1+row.size() ; i++) {
							final int refpos = i + read_start;
							final int textpos = refpos - startpos - matrix.size();
							final boolean in_ref = refpos>=0 && refpos< genomic.length();
							char c = Character.toUpperCase(in_ref?genomic.charAt(refpos):'N');
							char qual = (char)('I'+rnd.nextInt(5));
							final CigarOperator op = in_ref?CigarOperator.M:CigarOperator.S;
							
							if(i>=0 && i < matrixWidth && textpos>=0 && textpos<row.size() && row.get(textpos)) {
								nm++;
								if(use_snv_instead_of_del && AcidNucleics.isATGC(c)) {
									char compl;
									switch(c) {
										case 'A': compl = 'T' ; break;
										case 'T': compl = 'A' ; break;
										case 'G': compl = 'C' ; break;
										case 'C': compl = 'A' ; break;
										case 'N': compl = 'A' ; break;
										default: throw new IllegalStateException();
										}
									basesBuilder.append(compl);
									qualBuilder.append(qual);
									operators.add(CigarOperator.X);
									}
								else
									{
									operators.add(CigarOperator.D);
									}
								}
							else
								{
								basesBuilder.append(c);
								qualBuilder.append(qual);
								operators.add(op);
								}
							}
						rec.setMappingQuality(60);
						rec.setReadString(basesBuilder.toString());
						rec.setBaseQualityString(qualBuilder.toString());
						rec.setCigar(Cigar.fromCigarOperators(operators));
						rec.setAttribute(SAMTag.RG, rg.getId());
						rec.setAttribute(SAMTag.NM, nm);
						w.addAlignment(rec);
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
	
	public static void main(final String[] args) {
		new TextBam().instanceMainWithExit(args);
	}
}
