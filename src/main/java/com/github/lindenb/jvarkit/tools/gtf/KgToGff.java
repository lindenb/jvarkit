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
package com.github.lindenb.jvarkit.tools.gtf;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import htsjdk.tribble.gff.Gff3Writer;
/**

BEGIN_DOC

# Warning

ouput is escaped for UTF8, some characters like ':' might be converted to a hexadecimal encoding.
Frame is missing/wrong for CDS. I need to fix this.

# Example

```
$ head -n1  ~/jeter.kg | java -jar dist/kg2gff.jar   | sed 's/%3A/:/g'  

##gff-version 3.1.25
chr22	ucsc	gene	18317101	18507325	.	-	.	ID=gene:uc002znh.2.g1;Name=uc002znh.2.1;biotype=protein_coding;gene_id=uc002znh.2.g1
chr22	ucsc	mRNA	18317101	18507325	.	-	.	ID=transcript:uc002znh.2.g1;Parent=gene:uc002znh.2.g1;Name=uc002znh.2.1;biotype=protein_coding;transcript_id=uc002znh.2.t1
chr22	ucsc	exon	18317101	18317315	.	-	.	ID=uc002znh.2:E0;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E0
chr22	ucsc	CDS	18317216	18317315	.	-	1	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS0;protein_id=uc002znh.2
chr22	ucsc	exon	18324588	18324783	.	-	.	ID=uc002znh.2:E1;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E1
chr22	ucsc	CDS	18324588	18324783	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS1;protein_id=uc002znh.2
chr22	ucsc	exon	18347665	18347752	.	-	.	ID=uc002znh.2:E2;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E2
chr22	ucsc	CDS	18347665	18347752	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS2;protein_id=uc002znh.2
chr22	ucsc	exon	18348690	18348778	.	-	.	ID=uc002znh.2:E3;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E3
chr22	ucsc	CDS	18348690	18348778	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS3;protein_id=uc002znh.2
chr22	ucsc	exon	18354603	18354789	.	-	.	ID=uc002znh.2:E4;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E4
chr22	ucsc	CDS	18354603	18354789	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS4;protein_id=uc002znh.2
chr22	ucsc	exon	18368644	18368817	.	-	.	ID=uc002znh.2:E5;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E5
chr22	ucsc	CDS	18368644	18368817	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS5;protein_id=uc002znh.2
chr22	ucsc	exon	18369936	18369998	.	-	.	ID=uc002znh.2:E6;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E6
chr22	ucsc	CDS	18369936	18369998	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS6;protein_id=uc002znh.2
chr22	ucsc	exon	18370089	18370201	.	-	.	ID=uc002znh.2:E7;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E7
chr22	ucsc	CDS	18370089	18370201	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS7;protein_id=uc002znh.2
chr22	ucsc	exon	18371800	18371996	.	-	.	ID=uc002znh.2:E8;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E8
chr22	ucsc	CDS	18371800	18371996	.	-	1	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS8;protein_id=uc002znh.2
chr22	ucsc	exon	18374251	18374398	.	-	.	ID=uc002znh.2:E9;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E9
chr22	ucsc	CDS	18374251	18374398	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS9;protein_id=uc002znh.2
chr22	ucsc	exon	18376574	18376670	.	-	.	ID=uc002znh.2:E10;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E10
chr22	ucsc	CDS	18376574	18376670	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS10;protein_id=uc002znh.2
chr22	ucsc	exon	18378050	18378176	.	-	.	ID=uc002znh.2:E11;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E11
chr22	ucsc	CDS	18378050	18378176	.	-	1	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS11;protein_id=uc002znh.2
chr22	ucsc	exon	18379012	18379127	.	-	.	ID=uc002znh.2:E12;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E12
chr22	ucsc	CDS	18379012	18379127	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS12;protein_id=uc002znh.2
chr22	ucsc	exon	18379490	18379747	.	-	.	ID=uc002znh.2:E13;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E13
chr22	ucsc	CDS	18379490	18379747	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS13;protein_id=uc002znh.2
chr22	ucsc	exon	18382214	18382314	.	-	.	ID=uc002znh.2:E14;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E14
chr22	ucsc	CDS	18382214	18382314	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS14;protein_id=uc002znh.2
chr22	ucsc	exon	18383608	18383763	.	-	.	ID=uc002znh.2:E15;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E15
chr22	ucsc	CDS	18383608	18383763	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS15;protein_id=uc002znh.2
chr22	ucsc	exon	18384644	18384745	.	-	.	ID=uc002znh.2:E16;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E16
chr22	ucsc	CDS	18384644	18384745	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS16;protein_id=uc002znh.2
chr22	ucsc	exon	18385397	18385513	.	-	.	ID=uc002znh.2:E17;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E17
chr22	ucsc	CDS	18385397	18385513	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS17;protein_id=uc002znh.2
chr22	ucsc	exon	18387398	18387605	.	-	.	ID=uc002znh.2:E18;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E18
chr22	ucsc	CDS	18387398	18387605	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS18;protein_id=uc002znh.2
chr22	ucsc	exon	18389315	18389652	.	-	.	ID=uc002znh.2:E19;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E19
chr22	ucsc	CDS	18389315	18389578	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS19;protein_id=uc002znh.2
chr22	ucsc	exon	18507047	18507325	.	-	.	ID=uc002znh.2:E20;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E20
```

END_DOC
**/
@Program(
		name="kg2gff",
		description="Convert UCSC genpred file to gff3",
		creationDate="20210106",
		modificationDate="20230817",
		keywords= {"gff","gff3","ucsc","genpred"},
		jvarkit_amalgamion = true
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
	

	private static Map<String,List<String>> convertMap(final Map<String,String> map) {
		final Map<String,List<String>> atts = new LinkedHashMap<>(map.size());
		for(final String key: map.keySet()) {
			atts.put(key,Collections.singletonList(map.get(key)));
			}
		return atts;
		}
	
	private void process(final Gff3Writer out,final List<UcscTranscript> kgs) {
		try {
		final char delim = ':';
		final double UNDEFINED_SCORE=-1;
		final int UNDEFINED_PHASE=-1;
		final Map<String,String> atts = new LinkedHashMap<>();
		final int lclid = (++ID_GENERATOR);
		final UcscTranscript first = kgs.get(0);
		
		final String bioType = !first.isProteinCoding()?"misc_RNA":"protein_coding";
		final String geneId = first.getName2();
		final Strand strand = Strand.decode(first.getStrand().encodeAsChar());
		final String protName = (StringUtils.isBlank(first.getName2())?first.getTranscriptId():first.getName2());
		atts.put("ID", "gene"+delim+geneId);
		atts.put("Name",protName+"."+lclid);
		atts.put("biotype",bioType);
		atts.put("gene_id",geneId);
		
		out.addFeature(new Gff3FeatureImpl(
				first.getContig(),
				this.source,
				"gene",
				kgs.stream().mapToInt(K->K.getTxStart()).min().getAsInt(),
				kgs.stream().mapToInt(K->K.getTxEnd()).max().getAsInt(),
				UNDEFINED_SCORE,
				strand,
				UNDEFINED_PHASE,
				convertMap(atts)
				));

		for(UcscTranscript kg:kgs) {
			atts.clear();
			atts.put("ID", "transcript"+delim+geneId);
			atts.put("Parent", "gene"+delim+geneId);
			atts.put("Name",protName+"."+lclid);
			atts.put("biotype",bioType);
			atts.put("transcript_id", kg.getTranscriptId()+".t"+lclid);
			out.addFeature(new Gff3FeatureImpl(
					kg.getContig(),
					this.source,
					"mRNA",
					kg.getTxStart(),
					kg.getTxEnd(),
					UNDEFINED_SCORE,
					strand,
					UNDEFINED_PHASE,
					convertMap(atts)
					));

			for(UcscTranscript.Exon exon: kg.getExons()) {
				atts.clear();
				atts.put("ID",kg.getTranscriptId()+delim+"E"+exon.getIndex());
				atts.put("Parent",  "transcript"+delim+geneId);
				atts.put("Name",kg.getTranscriptId());
				atts.put("biotype",bioType);
				atts.put("exon_id",kg.getTranscriptId()+delim+"E"+exon.getIndex());
				out.addFeature(new Gff3FeatureImpl(
						kg.getContig(),
						this.source,
						"exon",
						exon.getStart(),
						exon.getEnd(),
						UNDEFINED_SCORE,
						strand,
						UNDEFINED_PHASE,
						convertMap(atts)
						));
				}
			if(kg.isProteinCoding()) {
				for(UcscTranscript.CDS cds: kg.getCDSs()) {
						out.addFeature(new Gff3FeatureImpl(
								kg.getContig(),
								this.source,
								"CDS",
								cds.getStart(),
								cds.getEnd(),
								UNDEFINED_SCORE,
								strand,
								-1, //TODO
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
				
				try( Gff3Writer gffw = new Gff3Writer(super.openPathOrStdoutAsPrintStream(this.outputFile))) {
					for(List<UcscTranscript> L:name2kgs.values()) {
						process(gffw,L);
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
