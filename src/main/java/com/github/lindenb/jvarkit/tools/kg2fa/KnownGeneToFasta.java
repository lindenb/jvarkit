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
package com.github.lindenb.jvarkit.tools.kg2fa;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
/**
 
## Example

```
$ java -jar dist/kg2fa.jar -R human_g1k_v37.fasta -D  --case 3 -L 0 | cut -c 1-200 | head

>ENST00000456328.2 585|ENST00000456328.2|chr1|+|11868|14409|11868|11868|3|11868,12612,13220,|12227,12721,14409,|0|DDX11L1|none|none|-1,-1,-1,
GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTCCCGTGTCCTTTCCACCGGGCCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCAGGTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGGAGTGG
>ENST00000607096.1 585|ENST00000607096.1|chr1|+|30365|30503|30365|30365|1|30365,|30503,|0|MIR1302-11|none|none|-1,
GGATGCCCAGCTAGTTTGAATTTTAGATAAACAACGAATAATTTCGTAGCATAAATATGTCCCAAGCTTAGTTTGGGACATACTTATGCTAAAAAACATTATTGGTTGTTTATCTGAGATTCAGAATTAAGCATTTTA
>ENST00000417324.1 585|ENST00000417324.1|chr1|-|34553|36081|34553|34553|3|34553,35276,35720,|35174,35481,36081,|0|FAM138A|none|none|-1,-1,-1,
CACACAACGGGGTTTCGGGGCTGTGGACCCTGTGCCAGGAAAGGAAGGGCGCAGCTCCTGCAATGCGGAGCAGCCAGGGCAGTGGGCACCAGGCTTTAGCCTCCCTTTCTCACCCTACAGAGGGCAGGCCCTTCAGCTCCATTCTCCTCCAAGGCTGCAGAGGGGGCAGGAATTGGGGGTGACAGGAGAGCTGTAAGGTC
>ENST00000335137.3 585|ENST00000335137.3|chr1|+|69090|70008|69090|70008|1|69090,|70008,|0|OR4F5|cmpl|cmpl|0,
ATGGTGACTGAATTCATTTTTCTGGGTCTCTCTGATTCTCAGGAACTCCAGACCTTCCTATTTATGTTGTTTTTTGTATTCTATGGAGGAATCGTGTTTGGAAACCTTCTTATTGTCATAACAGTGGTATCTGACTCCCACCTTCACTCTCCCATGTACTTCCTGCTAGCCAACCTCTCACTCATTGATCTGTCTCTGTC
>ENST00000466430.1 585|ENST00000466430.1|chr1|-|89294|120932|89294|89294|4|89294,92090,112699,120774,|91629,92240,112804,120932,|0|RP11-34P13.7|none|none|-1,-1,-1,-1,
CTGATCCATATGAATTCCTCTTATTAAGAAAAATAAAGCATCCAGGATTCAATGAAGAACTGACTATCACCTTGTTAATCATTCAGAAACATGTTGCAGGCTTAAGCCATTTTTGATATAGATACTGAAACAATTACTTGCTAAGAGCAAACTTGAAGgtatggataaggccctgagtcatcttcctgagctgaatgata
(...)

```
 
 
*/

@Program(name="kg2fa",
description="convert knownGenes to fasta",
keywords={"kg","knownGene","fasta"},
creationDate="2019-02-13"
)
public class KnownGeneToFasta
extends Launcher
{
private static final Logger LOG = Logger.build(KnownGeneToFasta.class).make();

@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
private File faidx = null;
@Parameter(names={"-D","--default"},description="Use default Known Gene source from UCSC.")
private boolean use_default_uri=false;
@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private File outputFile = null;
@Parameter(names={"--coding"},description="ignore non-coding transcripts.")
private boolean onlyCodingTranscript=false;
@Parameter(names={"-L"},description="fasta line length.")
private int fastaLineLen=50;
@Parameter(names={"--introns","--intron"},description="Remove introns")
private boolean remove_introns=false;
@Parameter(names={"--utrs","--utr"},description="Remove UTRs")
private boolean remove_utrs=false;
@Parameter(names={"--case","--style"},description="style: (0) do nothing (1): all UPPERCASE (2): all lowercase (3): exon UPPERCASE + intron LOWERCASE . Otherwise do nothing")
private int case_style = 0;
@Parameter(names={"--empty"},description="Discard empty files.")
private boolean prevent_empty = false;;
@Parameter(names={"--dict"},description="Write optional dict file")
private Path outputDict=null;



	@Override
	public int doWork(final List<String> args) {
		BufferedReader br=null;
		IndexedFastaSequenceFile indexedFastaSequenceFile = null;
		GenomicSequence genomicSequence = null;
		PrintWriter pw = null;
		BufferedWriter dictWriter =null;
		try {
			indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
			final ContigNameConverter refCtgNameConverter= ContigNameConverter.fromOneDictionary(refDict);

			
			
			final String input = oneFileOrNull(args);
			if(StringUtils.isBlank(input) && this.use_default_uri) {
				final String uri = KnownGene.getDefaultUri();
				LOG.info("default uri is "+uri);
				br= IOUtils.openURIForBufferedReading(uri);
				}
			else if(StringUtils.isBlank(input)) {
				br = IOUtils.openStreamForBufferedReader(stdin());
				}
			else if(this.use_default_uri)
				{
				LOG.error("input set and use_default_uri both set.");
				return -1;
				}
			else
				{
				br = IOUtils.openURIForBufferedReading(input);
				}
			
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			
			final SAMSequenceDictionaryCodec dictCodec;
			if(this.outputDict==null)
				{
				dictWriter=null;
				dictCodec=null;
				}
			else
				{
				dictWriter = Files.newBufferedWriter(this.outputDict);
				dictCodec = new  SAMSequenceDictionaryCodec(dictWriter);
				dictCodec.encodeHeaderLine(false);
				}
			
			final CharSplitter tab=CharSplitter.TAB;
			String line;
			final Set<String> contig_not_found = new HashSet<>();
			while((line=br.readLine())!=null) {
				if(pw.checkError()) break;
				if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
				final String tokens[]=tab.split(line);
				final KnownGene kg=new KnownGene(tokens);

				if((this.onlyCodingTranscript) && kg.getCdsStart()==kg.getCdsEnd())continue; 
				final String ctg = refCtgNameConverter.apply(kg.getContig());
				if(StringUtils.isBlank(ctg)) {
					if(contig_not_found.add(kg.getContig()))
						{
						LOG.warn(JvarkitException.ContigNotFoundInDictionary.getMessage(kg.getContig(), refDict));
						}
					continue;
					}
				kg.setChrom(ctg);
				
				final SAMSequenceRecord ssr = refDict.getSequence(ctg);
				if(kg.getTxEnd()>ssr.getSequenceLength()) {
					LOG.warn("knowngene "+kg+" ends "+kg.getTxEnd()+" beyond chromosome "+ctg+" length:"+ssr.getSequenceLength()+". Wrong REF ?");
					continue;
					}
				if( genomicSequence==null || !genomicSequence.getChrom().equals(ctg)) {
					/* now, we can change genomicSequence */
					genomicSequence = new GenomicSequence(indexedFastaSequenceFile, ctg);
					}
				int nBases=0;
				
				final Function<KnownGene,String> fastaHeader = K ->{
					final StringBuilder sb=new StringBuilder(">");
					sb.append(kg.getName());
					sb.append(" ");
					sb.append(String.join("|", tokens));
					return sb.toString();
					};
				
				
				if(!prevent_empty) {
					pw.print(fastaHeader.apply(kg));
				}
				if(this.fastaLineLen<=0) pw.println();
				
				
				for(int segIdx=0;segIdx< kg.getExonCount();++segIdx)
					{
					final int exon_index = kg.isPositiveStrand()?segIdx:(kg.getExonCount()-1)-segIdx;
					
					//exon and then intron
					for(int side=0;side<2;++side) {
						if(side==1 && this.remove_introns) continue;
						
						final KnownGene.Segment segment;
						//exon
						if(side==0) {
							segment = kg.getExon(exon_index);
							}
						else // intron
							{
							if(kg.isPositiveStrand()) {
								if(exon_index>=kg.getIntronCount()) continue;
								segment = kg.getIntron(exon_index);
								}
							else
								{
								if(exon_index==0) continue;
								segment = kg.getIntron(exon_index-1);
								}
							}
						
						
						final int exonLength = segment.getEnd()-segment.getStart();
						for(int x=0;x< exonLength;++x)
							{
							final int gpos0;
							if(kg.isPositiveStrand()) {
								gpos0 = segment.getStart()+x;
								}
							else
								{
								gpos0 = (segment.getEnd()-1)-x;
								}
							if(this.remove_utrs) {
								if(gpos0 < kg.getCdsStart()) continue;
								if(gpos0 >= kg.getCdsEnd()) continue;
							}
							
							char base= genomicSequence.charAt(gpos0);
							
							if(kg.isNegativeStrand()) {
								base=AcidNucleics.complement(base);
								}
							
							switch(this.case_style)
								{
								case 0: break;
								case 1: base= Character.toUpperCase(base);break;
								case 2: base= Character.toLowerCase(base);break;
								case 3: base= side==0?Character.toUpperCase(base):Character.toLowerCase(base);break;
								default:break;
								}
							
							if(base==0 && this.prevent_empty) {
								pw.print(fastaHeader.apply(kg));
								}
							if(this.fastaLineLen>0 && nBases%this.fastaLineLen==0) pw.println();
							pw.print(base);
							nBases++;
							}//end for x
						}//end for side
					}//end for exonIdx
				pw.println();
				if(nBases>0 || !this.prevent_empty) 
					{
					final SAMSequenceRecord ssrOut =new SAMSequenceRecord(kg.getName(), nBases);
					dictCodec.encodeSequenceRecord(ssrOut);
					}
				}
			pw.flush();
			pw.close();
			pw=null;
			if(dictWriter!=null) {
				dictWriter.flush();
				dictWriter.close();
				dictWriter=null;
				}

			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		}
		finally {
			CloserUtil.close(pw);
			CloserUtil.close(dictWriter);
			CloserUtil.close(br);
			CloserUtil.close(indexedFastaSequenceFile);
		}
	}

public static void main(final String[] args) {
	new KnownGeneToFasta().instanceMainWithExit(args);
	}
}
