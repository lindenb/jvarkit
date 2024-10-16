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
package com.github.lindenb.jvarkit.tools.gff2fa;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.tools.kg2fa.KnownGeneToFasta;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Constants;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

/**

BEGIN_DOC

Experimental.
Translates the output of bcftools consensus.
Output is a BED file containing the cDNA, the translated peptide and the mean codon usage .

## Example
```
samtools faidx ref.fa 8:11870-11890 |\
	bcftools consensus in.vcf.gz  -s SAMPLE -H A --include 'TYPE="snp"' |\
	java -jar dist/jvarkit.jar translategff3 --gff3 annot.gff3
```
END_DOC

*/
@Program(name="gff2fasta",
description="extract fasta from gtf",
keywords={"gff","gff3","fasta"},
creationDate="20241016",
modificationDate="20241016",
jvarkit_amalgamion = true
)
public class Gff3ToFasta extends Launcher {
	private static final Logger LOG = Logger.build(KnownGeneToFasta.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referencePath = null;

	
	private void echo(PrintWriter w, boolean reverse,FunctionalMap<String, String> props, String title,CharSequence seq) {
		w.print(">");
		w.print(props.entrySet().stream().map(KV->KV.getKey()+":"+KV.getValue()).collect(Collectors.joining("|")));
		if(reverse) seq=AcidNucleics.reverseComplement(seq);
		for(int i=0;i< seq.length();i++) {
			if(i%60==0) w.println();
			w.append(seq.charAt(i));
			}
		w.println();
		}
	

	
	private void recursive_transcript(Gff3Feature gene,Gff3Feature transcript,final PrintWriter pw,final ReferenceSequenceFile ref) {
		FunctionalMap<String, String> props  =FunctionalMap.make();
		
		boolean reverse= transcript.getStrand().equals(Strand.NEGATIVE);
		List<Gff3Feature> exons = new ArrayList<>();
		List<Gff3Feature> cdss = new ArrayList<>();
		for(Gff3Feature c:transcript.getChildren()) {
			if(c.getType().endsWith("exon")) exons.add(c);
			else if(c.getType().endsWith("CDS")) cdss.add(c);
			}
		Collections.sort(exons,LocatableUtils::compareTo);
		Collections.sort(cdss,LocatableUtils::compareTo);
		
		StringBuilder mRNA=new StringBuilder(exons.stream().mapToInt(E->E.getLengthOnReference()).sum());
		StringBuilder cdna=new StringBuilder(cdss.stream().mapToInt(E->E.getLengthOnReference()).sum());
		for(final Gff3Feature exon: exons) {
			mRNA.append(ref.getSubsequenceAt(exon.getContig(), exon.getStart(), exon.getEnd()));
			}
		echo(pw,reverse,props,"",mRNA);
		
		
		for(final Gff3Feature cds: cdss) {
			cdna.append(ref.getSubsequenceAt(cds.getContig(), cds.getStart(), cds.getEnd()));
			}
		echo(pw,reverse,props,"",cdna);
		
		if(!cdss.isEmpty()) {
			
			// [EX   ]===[    EXON  ]======
			// =============[ CDS   ]======
			// [UTR  ]=================
			StringBuilder leftutr=new StringBuilder();
			final Locatable firstCDS = cdss.get(0);
			for(Locatable exon:exons) {
				if(exon.getStart() > firstCDS.getEnd()) break;
				leftutr.append(ref.getSubsequenceAt(transcript.getContig(),
						exon.getStart(),
						Math.min(exon.getEnd(),firstCDS.getStart()-1)
						));
				}
			echo(pw,reverse,props,"",leftutr);

			// ======[    EXON  ]== [EXON]=
			// =========[ CDS ]============
			//======================[ UTR]===

			StringBuilder rightUTR=new StringBuilder();
			Locatable lastCDS = cdss.get(cdss.size()-1);
			for(Locatable exon:exons) {
				if(exon.getEnd() <= lastCDS.getEnd()) continue;
				rightUTR.append(ref.getSubsequenceAt(transcript.getContig(),
						Math.max(exon.getStart(),lastCDS.getEnd()+1),
						exon.getEnd()
						));
				}

			}
		}
	
	private void recursive(Gff3Feature feat,final PrintWriter pw,final ReferenceSequenceFile ref) {
		if(feat.getType().equals("gene")) {
			for(Gff3Feature c:feat.getChildren()) {
				recursive_transcript(feat,c,pw,ref);
				}
			}
		else
			{
			for(Gff3Feature c:feat.getChildren()) {
				recursive(c,pw,ref);
				}
			}
		}

	
	@Override
	public int doWork(List<String> args) {
		
		try {
		final String input = oneFileOrNull(args);
		final Set<String> keepAnnots  = new HashSet<>(Arrays.asList(
				Gff3Constants.ID_ATTRIBUTE_KEY,
				Gff3Constants.PARENT_ATTRIBUTE_KEY,
				Gff3Constants.NAME_ATTRIBUTE_KEY,
				"gene_name","gene_id","transcript_id","gene_biotype"
				));
		
		
		try(ReferenceSequenceFile fasta =ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referencePath)) {
			final Gff3Codec codec = new Gff3Codec(DecodeDepth.DEEP,ATT->!keepAnnots.contains(ATT));
		    LineIterator lr = IOUtils.openURIForLineIterator(input);
		    codec.readHeader(lr);
	    	try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
		        while(!codec.isDone(lr)) {
		                final Gff3Feature feat = codec.decode(lr);
		                if(feat==null) continue;
		                recursive(feat, pw,fasta);
		                }
		        pw.flush();
		    	}
	        codec.close(lr);
			}
	    
	    return 0;
		} catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new Gff3ToFasta().instanceMainWithExit(args);
	}

}
