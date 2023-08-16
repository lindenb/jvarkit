/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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

import java.io.File;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptReader;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
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
description="convert ucsc genpred to fasta",
keywords={"kg","knownGene","fasta","genpred"},
creationDate="20190213",
modificationDate="20230815"
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
@Parameter(names={"--empty"},description="Discard empty sequences.")
private boolean prevent_empty = false;
@Parameter(names={"-sql","--sql"},description= "SQL Schema URI. "+UcscTranscriptReader.SQL_DESC)
private String sqlUri="";


	private void writeFasta(PrintWriter pw, final Map<String,String> props, CharSequence seq) {
		if(seq.length()==0 && prevent_empty) return;
		props.put("length",String.valueOf(seq.length()));
		boolean first = true;
		for(String k: props.keySet()) {
			pw.print(first?">":"|");
			pw.print(k);
			pw.print(":");
			pw.print(props.get(k));
			first=false;
			}
		
		for(int i=0;i< seq.length();i++) {
			if(i%this.fastaLineLen==0) pw.println();
			pw.print(seq.charAt(i));
			}
		pw.println();
		}


	@Override
	public int doWork(final List<String> args) {
		try(ReferenceSequenceFile indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx)) {
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
			GenomicSequence genomicSequence = null;
			final String input = super.oneFileOrNull(args);
			try(PrintWriter pw = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
			final Set<String> contig_not_found = new HashSet<>();
			try(CloseableIterator<UcscTranscript> iter = StringUtils.isBlank(input)?
					UcscTranscriptCodec.makeIterator(stdin(),this.sqlUri):
					UcscTranscriptCodec.makeIterator(input,this.sqlUri)) {
				while(iter.hasNext()) {
				final UcscTranscript kg= iter.next();
				
				if((this.onlyCodingTranscript) && !kg.isProteinCoding()) continue; 
				final SAMSequenceRecord ssr = refDict.getSequence(kg.getContig());
				if(ssr==null) {
					if(!contig_not_found.contains(kg.getContig())) LOG.info("contig not found in dict:" + kg.getContig());
					contig_not_found.add(kg.getContig());
					continue;
					}
				if(kg.getTxEnd() > ssr.getSequenceLength()) {
					LOG.warn("knowngene "+kg+" ends "+kg.getTxEnd()+" beyond chromosome "+kg.getContig()+" length:"+ssr.getSequenceLength()+". Wrong REF ?");
					continue;
					}
				if( genomicSequence==null || !genomicSequence.hasName(kg.getContig())) {
					/* now, we can change genomicSequence */
					genomicSequence = new GenomicSequence(indexedFastaSequenceFile, kg.getContig());
					}
				Map<String,String> props = new LinkedHashMap<>();
				props.put("id",kg.getTranscriptId());
				if(!StringUtils.isBlank(kg.getName2())) props.put("gene",kg.getName2());
				props.put("strand",kg.getStrand().name());
				
				final UcscTranscript.MessengerRNA mRNA = kg.getMessengerRNA(genomicSequence);
				props.put("type","mRNA");
				
				writeFasta(pw,props,mRNA);
				if(mRNA.hasCodingRNA()) {
					final UcscTranscript.CodingRNA cDNA = mRNA.getCodingRNA();
					props.put("type","cDNA");
					writeFasta(pw,props,cDNA);
					final UcscTranscript.Peptide pep = cDNA.getPeptide();
					props.put("type","peptide");
					writeFasta(pw,props,pep);
					if(mRNA.hasUpstreamUntranslatedRNA()) {
						final UcscTranscript.UntranslatedRNA utr5 = mRNA.getUpstreamUntranslatedRNA();
						props.put("type","UTR5");
						writeFasta(pw,props,utr5);
						}
					if(mRNA.hasDownstreamUntranslatedRNA()) {
						final UcscTranscript.UntranslatedRNA utr3 = mRNA.getDownstreamUntranslatedRNA();
						props.put("type","UTR3");
						writeFasta(pw,props,utr3);
						}
					}
				} //end while(iter)
				} //end iter
			pw.flush();
			} // end pw
			return 0;
		} // end ref
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
	}

public static void main(final String[] args) {
	new KnownGeneToFasta().instanceMainWithExit(args);
	}
}
