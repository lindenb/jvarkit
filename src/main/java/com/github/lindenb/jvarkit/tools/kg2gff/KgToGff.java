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
package com.github.lindenb.jvarkit.tools.kg2gff;

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Constants;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import htsjdk.tribble.gff.Gff3Writer;
/**

BEGIN_DOC

# Warning

ouput is escaped for UTF8, some characters like ':' might be converted to a hexadecimal encoding.

# Example

```
headjeter.kg -n1 | java -jar dist/jvarkit.jar kg2gff 
##gff-version 3.1.25
chr1	ucsc	gene	55039456	55064852	.	+	.	ID=gene%3APCSK9;Name=PCSK9;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;gene_type=protein_coding
chr1	ucsc	mRNA	55039456	55064852	.	+	.	ID=ENST00000710286.1;Parent=gene%3APCSK9;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	exon	55039456	55040044	.	+	.	ID=ENST00000710286.1%3AE0;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE0
chr1	ucsc	exon	55043843	55044034	.	+	.	ID=ENST00000710286.1%3AE1;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE1
chr1	ucsc	exon	55046523	55046646	.	+	.	ID=ENST00000710286.1%3AE2;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE2
chr1	ucsc	exon	55052278	55052411	.	+	.	ID=ENST00000710286.1%3AE3;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE3
chr1	ucsc	exon	55052650	55052791	.	+	.	ID=ENST00000710286.1%3AE4;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE4
chr1	ucsc	exon	55055993	55056189	.	+	.	ID=ENST00000710286.1%3AE5;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE5
chr1	ucsc	exon	55057331	55057514	.	+	.	ID=ENST00000710286.1%3AE6;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE6
chr1	ucsc	exon	55058036	55058209	.	+	.	ID=ENST00000710286.1%3AE7;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE7
chr1	ucsc	exon	55058499	55058647	.	+	.	ID=ENST00000710286.1%3AE8;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE8
chr1	ucsc	exon	55059486	55059663	.	+	.	ID=ENST00000710286.1%3AE9;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE9
chr1	ucsc	exon	55061375	55061556	.	+	.	ID=ENST00000710286.1%3AE10;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE10
chr1	ucsc	exon	55063369	55064852	.	+	.	ID=ENST00000710286.1%3AE11;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE11
chr1	ucsc	CDS	55039481	55040044	.	+	0	ID=ENST00000710286.1%3AC1;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55043843	55044034	.	+	0	ID=ENST00000710286.1%3AC2;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55046523	55046646	.	+	0	ID=ENST00000710286.1%3AC3;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55052278	55052411	.	+	2	ID=ENST00000710286.1%3AC4;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55052650	55052791	.	+	0	ID=ENST00000710286.1%3AC5;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55055993	55056189	.	+	2	ID=ENST00000710286.1%3AC6;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55057331	55057514	.	+	0	ID=ENST00000710286.1%3AC7;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55058036	55058209	.	+	2	ID=ENST00000710286.1%3AC8;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55058499	55058647	.	+	2	ID=ENST00000710286.1%3AC9;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55059486	55059663	.	+	0	ID=ENST00000710286.1%3AC10;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55061375	55061556	.	+	2	ID=ENST00000710286.1%3AC11;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55063369	55063584	.	+	0	ID=ENST00000710286.1%3AC12;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	5_prime_UTR	55039456	55039480	.	+	.	ID=ENST00000710286.1%3AU55039456;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	3_prime_UTR	55063585	55064852	.	+	.	ID=ENST00000710286.1%3AU55063585;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
```

END_DOC
**/
@Program(
		name="kg2gff",
		description="Convert UCSC genpred/knowngene file to gff3 or gtf",
		creationDate="20210106",
		modificationDate="20250324",
		keywords= {"gff","gff3","ucsc","genpred"},
		biostars = {9610527},
		jvarkit_amalgamion = true
		)
public class KgToGff extends Launcher {
	private static final Logger LOG = Logger.of(KgToGff.class);
	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--source"},description= "label for column 'source' in the output")
	private String source = "ucsc";
	@Parameter(names={"--score"},description= "default score (-1: undefined)")
	private int score = -1;
	@Parameter(names={"--coding"},description= "select coding transcript only.")
	private boolean coding_only = false;
	@Parameter(names={"--gtf"},description= "use the GTF writer instead of GFF3")
	private boolean use_gtf_writer = false;
	@Parameter(names={"--modulo3"},description= "discard transcripts where CDS length isn't a modulo 3 (eg. remove xeno-transcript mapped on another build) ")
	private boolean check_cds_modulo_3 = false;

	
	private static int ID_GENERATOR=0;
	
	private static interface GXXWriter extends Closeable {
		public void addFeature(Gff3Feature feature)  throws IOException ;
		}
	
	
	private class GFFWriter implements GXXWriter {
		private final Gff3Writer delegate;
		GFFWriter(OutputStream os) {
			delegate = new Gff3Writer(os);
			}
		@Override
		public void addFeature(Gff3Feature feature) throws IOException {
			delegate.addFeature(feature);
			}
		@Override
		public void close() throws IOException {
			delegate.close();
			}
		}
	
	private class GTFWriter implements GXXWriter {
	    private final OutputStream out;
		GTFWriter(OutputStream os) {
			out=os;
			}
		@Override
		public void addFeature(Gff3Feature feature) throws IOException {
			write(escapeString(feature.getContig()));
			out.write('\t');
			write(escapeString(feature.getSource()));
			out.write('\t');
			write(escapeString(feature.getType()));
			out.write('\t');
			write(Integer.toString(feature.getStart()));
			out.write('\t');
			write(Integer.toString(feature.getEnd()));
			out.write('\t');
			write(feature.getScore() < 0 ? Gff3Constants.UNDEFINED_FIELD_VALUE : Double.toString(feature.getScore()));
			out.write('\t');
			out.write(feature.getStrand().toString().getBytes());
			out.write('\t');
			write(feature.getPhase() < 0 ? Gff3Constants.UNDEFINED_FIELD_VALUE : Integer.toString(feature.getPhase()));
			out.write('\t');
			if (feature.getAttributes().isEmpty()) {
	            out.write(Gff3Constants.UNDEFINED_FIELD_VALUE.getBytes());
				}
			else
				{
				boolean first=true;
				for(Map.Entry<String, List<String>> kv:feature.getAttributes().entrySet()) {
					if(!first)  out.write(' ');
					first=false;
					write(escapeString(kv.getKey()));
					for(String v:kv.getValue()) {
						out.write(' ');
						out.write(StringUtils.doubleQuote(escapeString(v)).getBytes());
						}
					
					out.write(';');
					}
				}
			out.write('\n');
			}
		
		private void write(final String s) throws IOException {
			out.write(escapeString(s).getBytes());
			}
		@Override
		public void close() throws IOException {
			out.flush();
			out.close();
			}
	    private String escapeString(final String s) {
	        try {
	            return URLEncoder.encode(s, "UTF-8").replace("+", " ");
	        	} catch (final UnsupportedEncodingException ex) {
	            throw new TribbleException("Encoding failure", ex);
	        	}
	    	}
		}

	private static Map<String,List<String>> convertMap(final Map<String,String> map) {
		final Map<String,List<String>> atts = new LinkedHashMap<>(map.size());
		for(final String key: map.keySet()) {
			atts.put(key,Collections.singletonList(map.get(key)));
			}
		return atts;
		}
	
	private void process(final GXXWriter out,final List<UcscTranscript> kgs) {
		try {
		final char delim = ':';
		final double score= this.score;
		final int UNDEFINED_PHASE=-1;
		final Map<String,String> atts = new LinkedHashMap<>();
		final UcscTranscript first = kgs.get(0);
		
		final String bioType = !first.isProteinCoding()?"misc_RNA":"protein_coding";
		final String gene_name = StringUtils.ifBlank(first.getName2(),"G"+(++ID_GENERATOR));
		final String gene_id = "GENE"+(++ID_GENERATOR);
		final Strand strand = Strand.decode(first.getStrand().encodeAsChar());
		atts.put("ID", gene_id);
		atts.put("Name",gene_name);
		atts.put("biotype",bioType);
		atts.put("gene_id",gene_id);
		atts.put("gene_name",gene_name);
		atts.put("gene_type",bioType);
		
		out.addFeature(new Gff3FeatureImpl(
				first.getContig(),
				this.source,
				"gene",
				kgs.stream().mapToInt(K->K.getTxStart()).min().getAsInt(),
				kgs.stream().mapToInt(K->K.getTxEnd()).max().getAsInt(),
				score,
				strand,
				UNDEFINED_PHASE,
				convertMap(atts)
				));

		for(UcscTranscript kg:kgs) {
			final String transcriptId = kg.getTranscriptId()+"."+(++ID_GENERATOR);
			atts.clear();
			atts.put("ID", transcriptId);
			atts.put("Parent", gene_id);
			atts.put("Name",transcriptId);
			atts.put("biotype",bioType);
			atts.put("gene_id",gene_id);
			atts.put("gene_name",gene_name);
			atts.put("transcript_id", transcriptId);
			atts.put("transcript_name", kg.getTranscriptId());
			out.addFeature(new Gff3FeatureImpl(
					kg.getContig(),
					this.source,
					(this.use_gtf_writer?"transcript": "mRNA"),
					kg.getTxStart(),
					kg.getTxEnd(),
					score,
					strand,
					UNDEFINED_PHASE,
					convertMap(atts)
					));

			for(UcscTranscript.Exon exon: kg.getExons()) {
				atts.clear();
				atts.put("ID",kg.getTranscriptId()+delim+"E"+exon.getIndex());
				atts.put("Parent",  transcriptId);
				atts.put("Name",kg.getTranscriptId());
				atts.put("biotype",bioType);
				atts.put("gene_id",gene_id);
				atts.put("gene_id",gene_id);
				atts.put("gene_name",gene_name);
				atts.put("transcript_id",transcriptId);
				atts.put("exon_id",kg.getTranscriptId()+delim+"E"+exon.getIndex());
				out.addFeature(new Gff3FeatureImpl(
						kg.getContig(),
						this.source,
						"exon",
						exon.getStart(),
						exon.getEnd(),
						score,
						strand,
						UNDEFINED_PHASE,
						convertMap(atts)
						));
				
				atts.remove("exon_id");
				atts.remove("Name");
				}
			if(kg.isProteinCoding()) {
				final List<UcscTranscript.CDS> cdsList = kg.getCDSs();
				if(check_cds_modulo_3) {
					if(cdsList.stream().mapToInt(CDS->CDS.getLengthOnReference()).sum()%3!=0) continue;
					}
				
				
				
				for(int x=0;x< cdsList.size();x++) {
						final UcscTranscript.CDS cds = cdsList.get(x);
						
						int n_codon=0;
						
						if(strand.equals(Strand.FORWARD)) {
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

						atts.put("ID","CDS"+(++ID_GENERATOR));
						out.addFeature(new Gff3FeatureImpl(
								kg.getContig(),
								this.source,
								"CDS",
								cds.getStart(),
								cds.getEnd(),
								score,
								strand,
								phase, 
								convertMap(atts)
								));
						}
				
					for(UcscTranscript.UTR utr : kg.getUTRs())  {
						atts.put("ID","UTR"+(++ID_GENERATOR));
						out.addFeature(new Gff3FeatureImpl(
								kg.getContig(),
								this.source,
								/** bcftools csq wants the 3|5_prime notation */
								utr.isUTR3()?(this.use_gtf_writer?"three_prime_utr":"3_prime_UTR"):(this.use_gtf_writer?"five_prime_utr":"5_prime_UTR"),
								utr.getStart(),
								utr.getEnd(),
								score,
								strand,
								UNDEFINED_PHASE,
								convertMap(atts)
								));
						}
					for(UcscTranscript.Codon codon : kg.getCodons())  {
						atts.put("ID","codon"+(++ID_GENERATOR));
						out.addFeature(new Gff3FeatureImpl(
								kg.getContig(),
								this.source,
								/** bcftools csq wants the 3|5_prime notation */
								(codon.isStartCodon()?"start_codon":"stop_codon"),
								codon.getStart(),
								codon.getEnd(),
								score,
								strand,
								UNDEFINED_PHASE,
								convertMap(atts)
								));
						}
					
				}
			}
		} catch(final Throwable err) {
			throw new RuntimeIOException(err);
			}
		}
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = super.oneFileOrNull(args);
			try(CloseableIterator<UcscTranscript> iter= (input!=null?UcscTranscriptCodec.makeIterator(input):UcscTranscriptCodec.makeIterator(stdin()))) {
				final Map<String,List<UcscTranscript>> name2kgs = new HashMap<>(60_000);
				while(iter.hasNext()) {
					final UcscTranscript kg=iter.next();
					if(coding_only && !kg.isProteinCoding()) continue;
					final String gene = kg.getContig()+"~"+kg.getName2();
					List<UcscTranscript> L = name2kgs.get(gene);
					if(L==null) {
						L = new ArrayList<>();
						name2kgs.put(gene,L);
						}
					L.add(kg);
					}
				
				try(OutputStream os=super.openPathOrStdoutAsPrintStream(this.outputFile)) {
					try( GXXWriter gffw = use_gtf_writer?new GTFWriter(os):new GFFWriter(os)) {
						for(List<UcscTranscript> L:name2kgs.values()) {
							process(gffw,L);
							}
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
		new KgToGff().instanceMainWithExit(args);

	}

}
