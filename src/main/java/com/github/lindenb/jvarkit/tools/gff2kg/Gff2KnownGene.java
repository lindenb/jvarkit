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
package com.github.lindenb.jvarkit.tools.gff2kg;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

/**
BEGIN_DOC

## Example

```
$  curl -s "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz" |\
	gunzip -c |\
	java -jar dist/jvarkit.jar gff2kg
(...)
1826	ENST00000367917.3	chr1	+	162760522	162782607	162760590	162782210	8	162760522,162762448,162766374,162767591,162769532,162774056,162775183,162782087	162760625,162762652,162766467,162767706,162769727,162774113,162775282,162782607	gene_id=ENSG00000132196.9;transcript_id=ENST00000367917.3;gene_type=protein_coding;gene_status=KNOWN;gene_name=HSD17B7;transcript_type=protein_coding;transcript_name=HSD17B7-201;protein_id=ENSP00000356894.3;havana_gene=OTTHUMG00000034420.6;	ENST00000367917.3
(...)
```

In the UCSC (not the structure of konwGene, but we can validate intervals):

```
$ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -D hg19 -e 'select * from wgEncodeGencodeBasicV19 where name="ENST00000367917.3"' | cat
bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
1826	ENST00000367917.3	chr1	+	162760522	162782607	162760590	162782210	8	162760522,162762448,162766374,162767591,162769532,162774056,162775183,162782087,	162760625,162762652,162766467,162767706,162769727,162774113,162775282,162782607,	0	HSD17B7	cmpl	cmpl	0,2,2,2,0,0,0,0,
```

### From ensembl 

```
$	wget -O -  "ftp://ftp.ensembl.org/pub/grch37/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz" |\
	gunzip -c |\
	java -jar dist/jvarkit.jar gff2kg
```

## see also

  * Ensembl vs UCSC  [https://twitter.com/yokofakun/status/743751004785545218](https://twitter.com/yokofakun/status/743751004785545218)



END_DOC
 */
@Program(name="gff2kg",
		description="Convert GFF3 format to UCSC knownGene format.",
		keywords={"gff",",gtf","knownGene","ucsc","convert"},
		creationDate="20160404",
		modificationDate="20220725",
		biostars=276099,
		jvarkit_amalgamion = true
		)
public class Gff2KnownGene extends Launcher {
	private static final Logger LOG = Logger.build(Gff2KnownGene.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-bin","--bin"},description="Insert  UCSC 'bin' column as the first column.")
	private boolean writeBin = false;
	@Parameter(names={"-bed12","--bed12"},description="Ouput bed.")
	private boolean writeBed12 = false;
	
	private static final Set<String> gene_keys = new HashSet<>(Arrays.asList(
		"gene_id","transcript_type","gene_name","gene_status",
		"gene_type","transcript_id","havana_gene","havana_transcript",
		"transcript_name","protein_id","ccdsid"
	   ));
	
    private void dump(final Gff3Feature feat,final PrintWriter pw) {
    	
    	if(feat.getType().equals("gene")) {
    		for(Gff3Feature transcript :feat.getChildren()) {
    			String firstTranscriptName = transcript.getAttribute("transcript_id").stream().findFirst().orElse(transcript.getID());
    			if(firstTranscriptName.startsWith("transcript:")) {
    				firstTranscriptName = firstTranscriptName.substring(11);
    				}
    			final List<Gff3Feature> exons = new ArrayList<>();
    			final List<Gff3Feature> cds = new ArrayList<>();
    			for(Gff3Feature f2 :transcript.getChildren()) {
    				if(f2.getType().equals("exon")) {
    					exons.add(f2);
    					}
    				else if(f2.getType().equals("CDS")) {
    					cds.add(f2);
    					}
        			}
    			if(exons.isEmpty()) continue;
    			exons.sort((o1,o2)->o1.getStart()-o2.getStart());
    			cds.sort((o1,o2)->o1.getStart()-o2.getStart());
    			if(this.writeBed12) {
    				final int transcript_start0 = transcript.getStart()-1;
    				pw.print(transcript.getContig());
					pw.print("\t");
					pw.print(transcript_start0);
					pw.print("\t");
					pw.print(transcript.getEnd());
					pw.print("\t");
					pw.print(firstTranscriptName);
					pw.print("\t");
					pw.print(1);//score
					pw.print("\t");
					switch(transcript.getStrand()) {
						case NEGATIVE: pw.print("-"); break;
						case POSITIVE: pw.print("+"); break;
						default: throw new IllegalArgumentException(transcript.toString());
						}
					pw.print("\t");
					pw.print(cds.stream().mapToInt(F->F.getStart()-1).min().orElse(transcript_start0));//thickStart
					pw.print("\t");
					pw.print(cds.stream().mapToInt(F->F.getEnd()).max().orElse(transcript_start0));//thickEnd
					pw.print("\t");
					pw.print("0,0,255");
					pw.print("\t");
					pw.print(exons.size());
					pw.print("\t");
					pw.print(exons.stream().map(E->String.valueOf(E.getLengthOnReference())).collect(Collectors.joining(",")));
					pw.print("\t");
					pw.print(exons.stream().map(E->String.valueOf((E.getStart()-1)-transcript_start0)).collect(Collectors.joining(",")));
					pw.println();
    			} else {
					if(this.writeBin) {
						pw.print(GenomicIndexUtil.regionToBin(transcript.getStart()-1, transcript.getEnd()));
						pw.print("\t");
						}
					pw.print(firstTranscriptName);
					pw.print("\t");
					pw.print(transcript.getContig());
					pw.print("\t");
					switch(transcript.getStrand()) {
						case NEGATIVE: pw.print("-"); break;
						case POSITIVE: pw.print("+"); break;
						default: throw new IllegalArgumentException(transcript.toString());
						}
					pw.print("\t");
					pw.print(transcript.getStart()-1);
					pw.print("\t");
					pw.print(transcript.getEnd());
					pw.print("\t");
					
					pw.print(cds.stream().mapToInt(F->F.getStart()-1).min().orElse(transcript.getStart()-1));
					pw.print("\t");
					pw.print(cds.stream().mapToInt(F->F.getEnd()).max().orElse(transcript.getStart()-1));
					pw.print("\t");
					
					pw.print(exons.size());
					pw.print("\t");
					pw.print(exons.stream().map(F->String.valueOf(F.getStart()-1)).collect(Collectors.joining(",")));
					pw.print("\t");
					pw.print(exons.stream().map(F->String.valueOf(F.getEnd())).collect(Collectors.joining(",")));
					pw.print("\t");
					
					
					for(String key : gene_keys)  {
						final List<String> values = transcript.getAttribute(key);
						if(values==null || values.size()!=1) continue;	
						pw.print(values.get(0));
						pw.print(";");
						}
					
					pw.print("\t");
					pw.print(firstTranscriptName);
					pw.println();
					} /* end not bed12 */
    			}
    		}
    	else {
    		for(Gff3Feature c:feat.getChildren()) {
    			dump(c,pw);
    			}
    	}
    }

	@Override
	public int doWork(final List<String> args) {		
		LineIterator in =null;
		try {
			final Gff3Codec gff3codec = new Gff3Codec(DecodeDepth.DEEP);
			final String input = oneFileOrNull(args);
			in = (input==null?
					IOUtils.openStdinForLineIterator():
					IOUtils.openURIForLineIterator(input)
					);
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				gff3codec.readHeader(in);
				while(!gff3codec.isDone(in)) {
					final Gff3Feature feat = gff3codec.decode(in);
					if(feat==null) continue;
					dump(feat,pw);
					}
				pw.flush();
				gff3codec.close(in);
				}
			return 0;
		} catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			}
		}
		
	public static void main(final String[] args) throws IOException
		{
		new Gff2KnownGene().instanceMainWithExit(args);
		}

	}
