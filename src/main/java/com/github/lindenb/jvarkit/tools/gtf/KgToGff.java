/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gtf;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene.CodingRNA;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import htsjdk.tribble.gff.Gff3Writer;
/**

BEGIN_DOC

```

```

END_DOC
**/
@Program(
		name="kg2gff",
		description="Convert UCSC knowGene file to gff",
		creationDate="20210106",
		modificationDate="20210106",
		generate_doc=false,
		keywords= {"gff","ucsc"}
		)
public class KgToGff extends Launcher {
	private static final Logger LOG = Logger.build(KgToGff.class).make();

	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--source"},description= "gff source field")
	private String source = "ucsc";
	@Parameter(names={"--coding"},description= "select coding transcript only.")
	private boolean coding_only = false;

	private static int ID_GENERATOR=0;
	
	private KgToGff() {}

	private static Map<String,List<String>> convertMap(final Map<String,String> map) {
		final Map<String,List<String>> atts = new LinkedHashMap<>(map.size());
		for(final String key: map.keySet()) {
			atts.put(key,Collections.singletonList(map.get(key)));
			}
		return atts;
		}
	
	private void process(final Gff3Writer out,final String line) {
		try {
		final char delim = ':';
		final double UNDEFINED_SCORE=-1;
		final int UNDEFINED_PHASE=-1;
		final String tokens[] = CharSplitter.TAB.split(line);
		final KnownGene kg = new KnownGene(tokens);
		if(coding_only && kg.isNonCoding()) return;
		final Map<String,String> atts = new LinkedHashMap<>();
		final int lclid = (++ID_GENERATOR);
		final String bioType = kg.isNonCoding()?"misc_RNA":"protein_coding";
		final String geneId = kg.getName()+".g"+ lclid;
		final Strand strand = Strand.decode(kg.getStrand().encodeAsChar());
		atts.put("ID", "gene"+delim+geneId);
		atts.put("Name",kg.getName());
		atts.put("biotype",bioType);
		atts.put("gene_id",geneId);
		
		out.addFeature(new Gff3FeatureImpl(
				kg.getContig(),
				this.source,
				"gene",
				kg.getTxStart()+1,
				kg.getTxEnd(),
				UNDEFINED_SCORE,
				strand,
				UNDEFINED_PHASE,
				convertMap(atts)
				));
		
		//final String transcriptId =;

		atts.clear();
		atts.put("ID", "transcript"+delim+geneId);
		atts.put("Parent", "gene"+delim+geneId);
		atts.put("Name",kg.getName());
		atts.put("biotype",bioType);
		atts.put("transcript_id", kg.getName()+".t"+lclid);
		out.addFeature(new Gff3FeatureImpl(
				kg.getContig(),
				this.source,
				"mRNA",
				kg.getTxStart()+1,
				kg.getTxEnd(),
				UNDEFINED_SCORE,
				strand,
				UNDEFINED_PHASE,
				convertMap(atts)
				));
		
		final CodingRNA cDNA = (kg.isNonCoding()?null:kg.getCodingRNA());
		
		for(KnownGene.Exon exon: kg.getExons()) {
			atts.clear();
			atts.put("ID",kg.getName()+delim+"E"+exon.getIndex());
			atts.put("Parent",  "transcript"+delim+geneId);
			atts.put("Name",kg.getName());
			atts.put("biotype",bioType);
			atts.put("exon_id",kg.getName()+delim+"E"+exon.getIndex());
			out.addFeature(new Gff3FeatureImpl(
					kg.getContig(),
					this.source,
					"exon",
					exon.getStart()+1,
					exon.getEnd(),
					UNDEFINED_SCORE,
					strand,
					UNDEFINED_PHASE,
					convertMap(atts)
					));
			if(cDNA!=null && !(exon.getEnd() <= kg.getCdsStart() || kg.getCdsEnd()<=exon.getStart())) {
				final int cdsStart = Math.max(exon.getStart(),kg.getCdsStart());
				final int cdsEnd = Math.min(exon.getEnd(),kg.getCdsEnd());
				if(cdsStart >= cdsEnd ) continue;
				atts.clear();
				atts.put("Parent", "transcript"+delim+geneId);
				atts.put("ID", "CDS"+delim+kg.getName()+delim+"CDS"+exon.getIndex());
				atts.put("protein_id",kg.getName());
				int i= 0;
				int firstExonPos=-1;
				int phase = -1;
				while(i< cDNA.length()) {
					int pos = cDNA.convertToGenomicCoordinate(i);
					if( cdsStart<=pos && pos < cdsEnd && firstExonPos==-1) {
						firstExonPos=pos;
						if( firstExonPos == -1) {
							LOG.error("Cannot extract first base of exon!!!");
							}
						 while(i< cDNA.length()) {
							pos = cDNA.convertToGenomicCoordinate(i);
							if(i%3==0) { phase = Math.abs( pos - firstExonPos ); break;}
							i++;
							}
						break;
						}
					i++;
					}
				if(phase<0  || phase>2) {
					LOG.warn("cannot get phase for "+kg.getName()+" "+exon.getName()+
						" kg="+ kg.getCdsStart()+"-"+kg.getCdsEnd()+
						" ex="+exon.getStart()+"-"+exon.getEnd()+
						" rn:"+cdsStart+"-"+cdsEnd +
						" first:"+firstExonPos +
						" phase="+phase+" leng="+cDNA.length()+" "+line);
					phase = UNDEFINED_PHASE;
					}
				out.addFeature(new Gff3FeatureImpl(
						kg.getContig(),
						this.source,
						"CDS",
						cdsStart+1,
						cdsEnd,
						UNDEFINED_SCORE,
						strand,
						phase,
						convertMap(atts)
						));
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
			try(BufferedReader br = super.openBufferedReader(input)) {
				try( Gff3Writer gffw = new Gff3Writer(super.openPathOrStdoutAsPrintStream(this.outputFile))) {
					br.lines().forEach(L->process(gffw,L));
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
