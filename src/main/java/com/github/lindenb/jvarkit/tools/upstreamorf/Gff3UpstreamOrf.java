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


Input is a UCSC "genpred/knowngene" file, but if you only have a gff/gff3 file, you can use `gff2kg` to create one.

## Examples


```
wget -O - "https://hgdownload.soe.ucsc.edu/goldenPath/hg38//database/wgEncodeGencodeBasicV47.txt.gz" |\
	gunzip -c |\
	java -jar dist/jvarkit.jar gff3upstreamorf  -R GRCh38.fa  > uorf.gff3

bcftools csq --ncsq 1000 -l -f GRCh38.fa -g uorf.gff3 in.bcf |\
	bcftools annotate --rename-annots <(echo -e 'INFO/BCSQ\tUTR_BCSQ')

	
```

note to self: test ENSG00000141736 https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003529

END_DOC

*/
@Program(name="gff3upstreamorf",
description="Takes a ucsc genpred file, scan the 5' UTRs and generate a GFF3 containing upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs ",
keywords={"gff","gff3","uorf","uorf"},
creationDate="20220724",
modificationDate="20250131",
jvarkit_amalgamion = true
)
public class Gff3UpstreamOrf extends Launcher
	{
	private static final Logger LOG = Logger.of(Gff3UpstreamOrf.class);
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"--strength"},description="only accept events that are greater or equal to this Kozak strength.")
	private KozakSequence.Strength user_kozak_strength = KozakSequence.Strength.nil;
	//@Parameter(names={"--break-original-orf"},description="if ATG(uORF) is in frame with original ORF , do not calculate the peptide beyond the original ATG.")
	//private boolean break_in_original_orf = false;

	private static int ID_GENERATOR=0;
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
	
	
	static Map<String,List<String>> gff3Map(final Map<String,Object> h) {
		return h.entrySet().stream().collect(
				Collectors.toMap(K->K.getKey(),
				K->Collections.singletonList(String.valueOf(K.getValue())))
				);
		}
	
	/** return true if standard cDNA contains the micro ORF */
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
		
		return atg1genome>= cDNA.getTranscript().getCdsStart() && 
			   atg1genome<= cDNA.getTranscript().getCdsEnd();
		}
	

	private void withMicroORF(final Tupple tupple ,Gff3Writer w) throws IOException {
		
		final String transcriptId= tupple.transcript.getTranscriptId()+".uORF."+tupple.kozak.getStrength().name()+"."+ tupple.id;
		final Map<String,Object> tr_attributes = new HashMap<>();
		tr_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, "transcript:"+ transcriptId);
		tr_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, transcriptId);
		tr_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY,"gene:"+tupple.geneName);

		tr_attributes.put("biotype", "protein_coding");
		tr_attributes.put("kozak", tupple.kozak);
		tr_attributes.put("kozak-strength", tupple.kozak.getStrength());
		tr_attributes.put("peptide", tupple.uPep.toString());
		tr_attributes.put("peptide-length", tupple.uPep.length());
		tr_attributes.put("transcript_id", tupple.transcript.getTranscriptId());
		
		final Gff3FeatureImpl trFeat = new Gff3FeatureImpl(
			tupple.contig,
			Gff3UpstreamOrf.class.getSimpleName(),
			"transcript",
			tupple.transcript.getStart(),
			tupple.transcript.getEnd(),
			kozakStrengthToScore(tupple.kozak.getStrength()),
			tupple.transcript.getStrand(),
			-1,
			gff3Map(tr_attributes)
			);
		w.addFeature(trFeat);

		for(UcscTranscript.Exon exon: tupple.uTranscript.getExons()) {
			//if( tupple.transcript.isPositiveStrand() && exon.getStart() > tupple.microORF.getEnd()) continue;
			//if( tupple.transcript.isNegativeStrand() && tupple.microORF.getStart() > exon.getEnd() ) continue;
			final Map<String,Object> ex_attributes = new HashMap<>();
			ex_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, tr_attributes.get(Gff3Constants.ID_ATTRIBUTE_KEY));
			ex_attributes.put("transcript_id", ex_attributes.get(Gff3Constants.PARENT_ATTRIBUTE_KEY));

			ex_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, "EXON"+(++ID_GENERATOR));
			ex_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, ex_attributes.get(Gff3Constants.ID_ATTRIBUTE_KEY));
			final Gff3FeatureImpl exFeat = new Gff3FeatureImpl(
					tupple.contig,
					Gff3UpstreamOrf.class.getSimpleName(),
					"exon",
					exon.getStart(),
					exon.getEnd(),
					0.0,
					exon.getStrand(),
					-1,
					gff3Map(ex_attributes)
					);
			w.addFeature(exFeat);
			}

		final List<UcscTranscript.CDS> cdsList = tupple.uTranscript.getCDS().stream().
				sorted((A,B)->Integer.compare(A.getStart(),B.getStart())).
				collect(Collectors.toList());


		int posInCDNA =0;
		for(int x=0; x<cdsList.size();++x) {
			final UcscTranscript.CDS cds = cdsList.get(x);
			if(!cds.overlaps(tupple.microORF)) throw new IllegalStateException("boum vs cds");
			final Map<String,Object> cds_attributes = new HashMap<>();
			cds_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, tr_attributes.get(Gff3Constants.ID_ATTRIBUTE_KEY));
			cds_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, "CDS"+(++ID_GENERATOR));
			cds_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY,cds_attributes.get(Gff3Constants.ID_ATTRIBUTE_KEY));
			cds_attributes.put("pos-in-cdna", posInCDNA);
			cds_attributes.put("pos-in-cdna-mod3", posInCDNA%3);
			cds_attributes.put("transcript_id", cds_attributes.get(Gff3Constants.PARENT_ATTRIBUTE_KEY));

			
			
			int n_codon=0;
			
			if(tupple.transcript.isPositiveStrand()) {
				for(int y=0;y< x;y++) {
					n_codon+= cdsList.get(y).getLengthOnReference();
					}
				}
			else
				{
				for(int y=cdsList.size()-1;y>x ;y--) {
					n_codon+= cdsList.get(y).getLengthOnReference();
					}
				}
			
			
			/**
			 * https://github.com/The-Sequence-Ontology/Specifications/issues/20
			 * 
			 * For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
			 * The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the 
			 * beginning of this feature to reach the first base of the next codon. 
			 * In other words, a phase of "0" indicates that the next codon begins at the first base of the region described 
			 * by the current line, a phase of "1" indicates that the next codon begins at the second base of this region, 
			 * and a phase of "2" indicates that the codon begins at the third base of this region. 
			 * This is NOT to be confused with the frame, which is simply start modulo 3.
			 */
			int phase=0;
			while((n_codon+phase)%3!=0) {
				phase++;
				}

			final Gff3FeatureImpl cdsfeat = new Gff3FeatureImpl(
				tupple.contig,
				Gff3UpstreamOrf.class.getSimpleName(),
				"CDS",
				cds.getStart() ,
				cds.getEnd() ,
				kozakStrengthToScore(tupple.kozak.getStrength()) ,
				tupple.uTranscript.getStrand() ,
				phase ,
				gff3Map(cds_attributes)
				);
			posInCDNA += cds.getLengthOnReference();
			w.addFeature(cdsfeat);
			}
		}

	/** process each transcripts in a gene */
	private void withGene(
			final ReferenceSequenceFile indexedFastaSequenceFile,
			final List<UcscTranscript> transcripts ,
			final UnaryOperator<String> ctgConverter,
			final Gff3Writer w
			)  throws IOException {
		if(transcripts.isEmpty()) return;
		final String ctg = ctgConverter.apply(transcripts.get(0).getContig());
		if(StringUtils.isBlank(ctg)) return ;
		final String gene = transcripts.get(0).getName2();
		if(genomicSequence==null || !genomicSequence.getContig().equals(ctg)) {
			genomicSequence= new GenomicSequence(indexedFastaSequenceFile, ctg);
			}
		final List<Tupple> microORFSs = new ArrayList<>();
		
		for(UcscTranscript tr: transcripts) {
			if(!tr.hasUTR5()) continue;
			final UcscTranscript.MessengerRNA mRNA = tr.getMessengerRNA(this.genomicSequence);
			final UcscTranscript.UntranslatedRNA upstream= mRNA.getUpstreamUntranslatedRNA();
			for(UcscTranscript.CodingRNA microORF : upstream.getORFs()) {
				boolean flag=true;
				// check this micro ORF is not the frame of another standard transcript 
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
		
		final Map<String,Object> g_atts = new HashMap<>();
		g_atts.put(Gff3Constants.ID_ATTRIBUTE_KEY, "gene:"+gene);
		g_atts.put("gene_name", gene);
		g_atts.put("biotype", "protein_coding");
		g_atts.put("gene_biotype", "protein_coding");
			
	
		
		final Gff3FeatureImpl gFeat = new Gff3FeatureImpl(
				ctg,
				Gff3UpstreamOrf.class.getSimpleName(),
				"gene",
				microORFSs.stream().mapToInt(T->T.microORF.getStart()).min().getAsInt(),
				microORFSs.stream().mapToInt(T->T.microORF.getEnd()).max().getAsInt(),
				0.0,
				transcripts.get(0).getStrand(),
				-1,
				gff3Map(g_atts)
				);
		w.addFeature(gFeat);
		
		for(Tupple t :microORFSs) {
			withMicroORF(t, w);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {		
		try {
			try(ReferenceSequenceFile indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx) ) {
				final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
				final ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(refDict);
				final String input = oneFileOrNull(args);
				
				final Map<String,List<UcscTranscript>> gene2transcripts = new HashMap<>(20_000);
				
				
				String prev_chr = "";
				try(CloseableIterator<UcscTranscript> iter = StringUtils.isBlank(input)?UcscTranscriptCodec.makeIterator(stdin()):UcscTranscriptCodec.makeIterator(input)) {
					while(iter.hasNext()) {
						final UcscTranscript tr = iter.next();
						if(!tr.isCoding()) continue;
						if(!tr.hasUTR5()) continue;
						if(StringUtils.isBlank(tr.getName2())) continue;
						final String key = tr.getContig()+"~"+tr.getName2();
						List<UcscTranscript> L = gene2transcripts.get(key);
						if(L==null) {
							L=new ArrayList<>();
							gene2transcripts.put(key,L);
							}
						L.add(tr);
						
						// keep memory low if possible
						if(!tr.getContig().equals(prev_chr)) {
							System.gc();
							prev_chr = tr.getContig();
							}
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
							withGene(
								indexedFastaSequenceFile,
								gene2transcripts.get(geneName),
								ctgConverter,
								w
								);
							}
						}
					}
				return 0;
				}
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new Gff3UpstreamOrf().instanceMainWithExit(args);
	}
}
