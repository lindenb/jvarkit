/*

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
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.upstreamorf;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.function.UnaryOperator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.KozakSequence;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.gff.Gff3Constants;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import htsjdk.tribble.gff.Gff3Writer;

/**
BEGIN_DOC

## inspiration

Part of this code was inspired from: https://github.com/ImperialCardioGenetics/uORFs/blob/master/5primeUTRannotator/five_prime_UTR_annotator.pm

Wikipedia:

> An Upstream Open Reading Frame (uORF) is an open reading frame (ORF) within the 5' untranslated region (5'UTR) of an mRNA. uORFs can regulate eukaryotic gene expression.
> Translation of the uORF typically inhibits downstream expression of the primary ORF. In bacteria, uORFs are called leader peptides, and were originally discovered on the basis of their impact on the regulation of genes involved in the synthesis or transport of amino acids. 

## Examples

### Example 1

```
 java -jar dist/gff3upstreamorfasta  -R GRCh38.fa Homo_sapiens.GRCh38.107.chr.gff3.gz > uorf.gff3
```

note to self: test ENSG00000141736 https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003529

END_DOC

*/
@Program(name="gff3upstreamorf",
description="Takes a ucsc genpred file, scan the 5' UTRs and generate a GFF3 containing upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs ",
keywords={"gff","gff3","uorf","uorf"},
creationDate="20220724",
modificationDate="20230820",
jvarkit_amalgamion = true
)
public class Gff3UpstreamOrf extends Launcher
	{
	private static final Logger LOG = Logger.build(Gff3UpstreamOrf.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"--strength"},description="only accept events that are greater or equal to this Kozak strength.")
	private KozakSequence.Strength user_kozak_strength = KozakSequence.Strength.nil;
	@Parameter(names={"--break-original-orf"},description="if ATG(uORF) is in frame with original ORF , do not calculate the peptide beyond the original ATG.")
	private boolean break_in_original_orf = false;

	//private static int ID_GENERATOR=0;
	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;

	private static class Tupple {
		String contig;
		String geneName;
		int id;
		UcscTranscript transcript;
		UcscTranscript.CodingRNA microORF;
		UcscTranscript.Peptide uPep;
		KozakSequence kozak;
		UcscTranscript uTranscript;
		}
	
	
	private double kozakStrengthToScore(final KozakSequence.Strength f) {
		switch(f)
			{
			case nil: return 0;
			case Weak: return 10;
			case Moderate:  return 100;
			case Strong:return 1000;
			default: throw new IllegalStateException();
			}
		}
	
	
	
		

	
	/** mutated version of UpstreamORF */
	
	private boolean acceptKozak(final KozakSequence k) {
		return k.getStrength().compareTo(this.user_kozak_strength)<=0;
		}
	
	static private boolean containsOpenReadingFrame(
			final UcscTranscript.CodingRNA cDNA,
			final UcscTranscript.CodingRNA uorf
			) {
		if(cDNA.getStart()==uorf.getStart()) return true;
		// other must have atg in frame and ATG before the observed one
		final int my_atg0 = cDNA.convertCoding0ToMessenger0(0);
		final int orf_atg0 = uorf.convertCoding0ToMessenger0(0);
		if(orf_atg0%3 != my_atg0%3) return false;
		if(orf_atg0 < my_atg0) return false;
		int atg1genome = uorf.convertToGenomic0Coordinate0(0)+1;
		
		return atg1genome>= cDNA.getTranscript().getCdsStart() && atg1genome<= cDNA.getTranscript().getCdsEnd();
		}
	
	/*
	private GFF3Transcript withTranscript(Gff3Feature transcript,final String fixedContig) {
		if(transcript.getStrand()==Strand.NONE) return null;
		final GFF3Transcript tr = new GFF3Transcript(transcript,fixedContig);
		for(final Locatable cds:tr.CDS) {
			if(!tr.exons.stream().anyMatch(E->E.contains(cds))) {
				LOG.error("in "+tr.getId()+" cds "+cds+" in not covered by some exon.");
				return null;
				}
			}
		
		return (tr.CDS.isEmpty() || tr.exons.isEmpty()) ? null: tr;
		}*/

	private void withMicroORF(final Tupple tupple ,Gff3Writer w) throws IOException {
		
		final String transcriptId= tupple.transcript.getTranscriptId()+"."+ tupple.id;
		final Map<String,List<String>> tr_attributes = new HashMap<>();
		tr_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("transcript:"+ transcriptId));
		tr_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(transcriptId));
		tr_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, Collections.singletonList("gene:"+tupple.geneName));

		tr_attributes.put("biotype", Collections.singletonList("protein_coding"));
		tr_attributes.put("kozak", Collections.singletonList(String.valueOf(tupple.kozak)));
		tr_attributes.put("kozak-strength", Collections.singletonList(String.valueOf(tupple.kozak.getStrength())));
		tr_attributes.put("peptide", Collections.singletonList(tupple.uPep.toString()));
		tr_attributes.put("peptide-length", Collections.singletonList(String.valueOf(tupple.uPep.length())));
		tr_attributes.put("transcript_id", Collections.singletonList(String.valueOf(tupple.transcript.getTranscriptId())));
		//tr_attributes.put("original-atg-pos", Collections.singletonList(String.valueOf(transcript.getATGPos1InGenome())));
		//tr_attributes.put("atg-pos", Collections.singletonList(String.valueOf(getATGPos1InGenome())));
		//tr_attributes.put("in-frame-with-transcript-atg", Collections.singletonList(String.valueOf(this.getATGPos0InRNA()%3==this.transcript.getATGPos0InRNA()%3)));

		final Gff3FeatureImpl trFeat = new Gff3FeatureImpl(
			tupple.contig,
			Gff3UpstreamOrf.class.getSimpleName(),
			"transcript",
			tupple.transcript.getStart(),
			tupple.transcript.getEnd(),
			kozakStrengthToScore(tupple.kozak.getStrength()),
			tupple.transcript.getStrand(),
			0,
			tr_attributes
			);
		w.addFeature(trFeat);

		for(UcscTranscript.Exon exon: tupple.uTranscript.getExons()) {
			//if( tupple.transcript.isPositiveStrand() && exon.getStart() > tupple.microORF.getEnd()) continue;
			//if( tupple.transcript.isNegativeStrand() && tupple.microORF.getStart() > exon.getEnd() ) continue;
			final Map<String,List<String>> ex_attributes = new HashMap<>();
			ex_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, Collections.singletonList(trFeat.getID()));
			ex_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("exon:"+exon.getName()));
			ex_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(exon.getName()));
			final Gff3FeatureImpl exFeat = new Gff3FeatureImpl(
					tupple.contig,
					Gff3UpstreamOrf.class.getSimpleName(),
					"exon",
					exon.getStart(),
					exon.getEnd(),
					0.0,
					exon.getStrand(),
					0,
					ex_attributes
					);
			w.addFeature(exFeat);
			}

		final List<UcscTranscript.CDS> cdss = tupple.uTranscript.getCDS().stream().sorted((A,B)->{
                        if(A.isPositiveStrand()) {
                                return Integer.compare(A.getStart(),B.getStart());
                                }
                        else
                                {
                                return Integer.compare(B.getStart(),A.getStart());
                                }
                        }).collect(Collectors.toList());

		for(int i=0;i+1< cdss.size();i++) {
			for(int j=i+1;j< cdss.size();j++) {
			if(cdss.get(i).overlaps(cdss.get(j))) throw new IllegalStateException("overlappng CDS:" +cdss.get(i)+" "+cdss.get(j));
			}
		}

		int posInCDNA =0;
		for(UcscTranscript.CDS cds: cdss) {
			if(!cds.overlaps(tupple.microORF)) throw new IllegalStateException("boum vs cds");
			//if( tupple.transcript.isPositiveStrand() && cds.getStart() > tupple.microORF.getEnd()) continue;
			//if( tupple.transcript.isNegativeStrand() && tupple.microORF.getStart() > cds.getEnd() ) continue;
			//if(cds.isNegativeStrand()) posInCDNA += (cds.getLengthOnReference()-1);
			final Map<String, List<String>> cds_attributes = new HashMap<>();
			cds_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, Collections.singletonList(trFeat.getID()));
			cds_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("cds:"+cds.getName()));
			cds_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(cds.getName()));
			cds_attributes.put("pos-in-cdna", Collections.singletonList(""+posInCDNA));
			cds_attributes.put("pos-in-cdna-mod3", Collections.singletonList(""+(posInCDNA%3)));

			int phase = posInCDNA%3;
			/* https://genomic.social/@scottcain/110923216702965764 */
			if(phase!=0) phase=3-phase;

			final Gff3FeatureImpl cdsfeat = new Gff3FeatureImpl(
				tupple.contig,
				Gff3UpstreamOrf.class.getSimpleName(),
				"CDS",
				cds.getStart() ,
				cds.getEnd() ,
				kozakStrengthToScore(tupple.kozak.getStrength()) ,
				tupple.uTranscript.getStrand() ,
				phase ,
				cds_attributes
				);
			posInCDNA += cds.getLengthOnReference();
			w.addFeature(cdsfeat);
			}
		}

	private void withGene(List<UcscTranscript> transcripts ,UnaryOperator<String> ctgConverter, final Gff3Writer w)  throws IOException {
		if(transcripts.isEmpty()) return;
		final String ctg = ctgConverter.apply(transcripts.get(0).getContig());
		if(StringUtils.isBlank(ctg)) return ;
		final String gene = transcripts.get(0).getName2();
		if(genomicSequence==null || !genomicSequence.getContig().equals(ctg)) {
			genomicSequence= new GenomicSequence(this.indexedFastaSequenceFile, ctg);
			}
		final List<Tupple> microORFSs = new ArrayList<>();
		for(UcscTranscript tr: transcripts) {
			final UcscTranscript.MessengerRNA mRNA = tr.getMessengerRNA(this.genomicSequence);
			if(!tr.hasUTR5()) continue;
			final UcscTranscript.UntranslatedRNA upstream= mRNA.getUpstreamUntranslatedRNA();
			for(UcscTranscript.CodingRNA microORF : upstream.getORFs()) {
				boolean flag=true;
				for(UcscTranscript tr2: transcripts) {
					if(tr2==tr) continue;
					if(!tr2.isCoding()) continue;
					final UcscTranscript.MessengerRNA mRNA2 = tr2.getMessengerRNA(this.genomicSequence);
					if(containsOpenReadingFrame(mRNA2.getCodingRNA(),microORF)) {
						flag=false;
						break;
						}
					}
				if(!flag) continue;
				final Tupple t  = new Tupple();
				t.kozak = microORF.getKozakSequence();
				if(!acceptKozak(t.kozak)) continue;
				t.geneName = gene;
				t.contig = ctg;
				t.transcript = tr;

				t.microORF = microORF;
				t.uTranscript = microORF.getTranscript();
				if(t.uTranscript == tr) throw new IllegalStateException("boum");
				t.uPep = microORF.getPeptide();
				t.id = microORFSs.size()+1;
				microORFSs.add(t);
				}
			}
		if(microORFSs.isEmpty()) return;
		
		final Map<String,List<String>> g_atts = new HashMap<>();
		g_atts.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("gene:"+gene));
		g_atts.put("gene_name", Collections.singletonList(gene));
			
	
		
		final Gff3FeatureImpl gFeat = new Gff3FeatureImpl(
				ctg,
				Gff3UpstreamOrf.class.getSimpleName(),
				"gene",
				microORFSs.stream().mapToInt(T->T.microORF.getStart()).min().getAsInt(),
				microORFSs.stream().mapToInt(T->T.microORF.getEnd()).max().getAsInt(),
				0.0,
				transcripts.get(0).getStrand(),
				0,
				g_atts
				);
		w.addFeature(gFeat);
		
		for(Tupple t :microORFSs) {
			withMicroORF(t, w);
		}
	}
	
	@Override
	public int doWork(final List<String> args) {		
		try {
			this.indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			final ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(refDict);
			final String input = oneFileOrNull(args);
			
			final Map<String,List<UcscTranscript>> gene2transcripts = new HashMap<>(20_000);
			try(CloseableIterator<UcscTranscript> iter = StringUtils.isBlank(input)?UcscTranscriptCodec.makeIterator(stdin()):UcscTranscriptCodec.makeIterator(input)) {
				while(iter.hasNext()) {
					final UcscTranscript tr = iter.next();
					if(!tr.hasUTR5()) continue;
					if(StringUtils.isBlank(tr.getName2())) continue;
					final String key = tr.getContig()+"~"+tr.getName2();
					List<UcscTranscript> L = gene2transcripts.get(key);
					if(L==null) {
						L=new ArrayList<>();
						gene2transcripts.put(key,L);
						}
					L.add(tr);
					}
				}
			if(gene2transcripts.isEmpty()) {
				LOG.error("No transcript was found");
				return -1;
				}
			try(OutputStream ps = super.openPathOrStdoutAsStream(this.outputFile)) {
				try(Gff3Writer w = new Gff3Writer(ps) {
					protected String escapeString(String s) {
						return s;
					}
				}) {
					for(final KozakSequence.Strength f:KozakSequence.Strength.values()) {
						w.addComment("kozak."+f.name()+"="+kozakStrengthToScore(f));
						}
					final String gtfSource = getProgramName().toLowerCase();
					for(final SAMSequenceRecord ssr: refDict.getSequences()) {
						w.addComment("contig "+ssr.getSequenceName()+": length:"+ssr.getSequenceLength());
					}
					w.addComment(gtfSource+":"+JVarkitVersion.getInstance().toString());
					if(!StringUtils.isBlank(input)) {
						w.addComment("genpred:"+input);
					}
					for(final String geneName: gene2transcripts.keySet()) {
						withGene(gene2transcripts.get(geneName),ctgConverter,w);
						}
					}
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
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		
		}
	
	public static void main(final String[] args) {
		new Gff3UpstreamOrf().instanceMainWithExit(args);
	}
}
