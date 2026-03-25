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
package com.github.lindenb.jvarkit.tools.gwascat2bed;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.util.Interval;
/**
 
BEGIN_DOC

## Example

```
 $ unzip -p  "gwas-catalog.zip "gwas-catalog-download-associations-alt-full.tsv" |\
 	java -jar dist/jvarkit.jar gwascatalog2bed --gtf  hs37d5.gtf.gz --header 2> /dev/null  |\
 	verticalize 

>>> 2
$1                      #contig : chr18
$2                        start : 11981023
$3                          end : 12030876
$4        DATE ADDED TO CATALOG : 2008-06-16
$5                     PUBMEDID : 17434096
$6                 FIRST AUTHOR : Matarin M
$7                         DATE : 2007-05-06
$8                      JOURNAL : Lancet Neurol
$9                         LINK : www.ncbi.nlm.nih.gov/pubmed/17434096
$10                       STUDY : A genome-wide genotyping study in patients with ischaemic stroke: initial analysis and data release.
$11               DISEASE/TRAIT : Stroke
$12         INITIAL SAMPLE SIZE : 249 European ancestry cases, 268 European ancestry controls
$13     REPLICATION SAMPLE SIZE : NA
$14                      REGION : 18p11.21
$15                      CHR_ID : 18
$16                     CHR_POS : 11987273
$17            REPORTED GENE(S) : IMPA2
$18                 MAPPED_GENE : IMPA2
$19            UPSTREAM_GENE_ID : 
$20          DOWNSTREAM_GENE_ID : 
$21                SNP_GENE_IDS : ENSG00000141401
$22      UPSTREAM_GENE_DISTANCE : 
$23    DOWNSTREAM_GENE_DISTANCE : 
$24   STRONGEST SNP-RISK ALLELE : rs7506045-?
$25                        SNPS : rs7506045
$26                      MERGED : 0
$27              SNP_ID_CURRENT : 7506045
$28                     CONTEXT : intron_variant
$29                  INTERGENIC : 0
$30       RISK ALLELE FREQUENCY : 0.10
$31                     P-VALUE : 7E-7
$32                 PVALUE_MLOG : 6.154901959985743
$33              P-VALUE (TEXT) : 
$34                  OR or BETA : 5.39
$35               95% CI (TEXT) : [2.77-10.5]
$36  PLATFORM [SNPS PASSING QC] : Illumina [408803]
$37                         CNV : N
$38                MAPPED_TRAIT : stroke
$39            MAPPED_TRAIT_URI : http://www.ebi.ac.uk/efo/EFO_0000712
$40             STUDY ACCESSION : GCST000032
$41       GENOTYPING TECHNOLOGY : Genome-wide genotyping array
<<< 2
```

END_DOC

*/
@Program(name="gwascatalog2bed",
	description="convert gwas catalog TSV association file to BED using a gtf file.",
	keywords={"bed","gwas","gwascatalog"},
	modificationDate="20260325",
	creationDate = "20260325"
	)
public class GwasCatalogToBed extends Launcher {
	private static final Logger LOG = Logger.of(GwasCatalogToBed.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output=null;
	@Parameter(names={"-G","--gtf"},description="GTF file",required = true)
	private Path gtf=null;
	@Parameter(names={"--header"},description="print header")
	private boolean with_header = false;
	@Parameter(names={"--pos","--position"},description="use CHR_ID/CHR_POS to get a position. Beware, check the build version.")
	private boolean enable_pos = false;

	public GwasCatalogToBed() {
	}

	private String removeVersion(final String s) {
		return s.trim().replaceAll("\\.[0-9]+$", "");
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Map<String,Interval> name2gene=new HashMap<>(50_000);
			boolean warning_position_printed=false;
			try(BufferedReader br = IOUtils.openPathForBufferedReading(gtf)) {
				final GTFCodec codec= new GTFCodec();
				for(;;) {
					String line=br.readLine();
					if(line==null) break;
					final GTFLine rec = codec.decode(line);
					if(rec==null) continue;
					if(!(rec.isTranscript() || rec.isGene())) continue;
					final String id ;
					if(rec.isTranscript()) {
						id = removeVersion(rec.getTranscriptId());
						}
					else
						{
						id = removeVersion(rec.getGeneId());
						}
					if(StringUtils.isBlank(id)) continue;
					name2gene.put(id, new Interval(rec));
					}
				}
			if(name2gene.isEmpty()) {
				LOG.error("No gene was found in "+this.gtf);
				return -1;
				}
			final ContigNameConverter ctgConverter = ContigNameConverter.fromContigSet(name2gene.values().stream().map(G->G.getContig()).collect(Collectors.toSet()));
			try(BufferedReader br = super.openBufferedReader(super.oneFileOrNull(args))) {
				String line = br.readLine();
				if(StringUtils.isBlank(line)) throw new IOException("Cannot read first line of gwas cat");
				final FileHeader header = new FileHeader(line, CharSplitter.TAB);
				header.assertColumnExists("UPSTREAM_GENE_ID");
				header.assertColumnExists("DOWNSTREAM_GENE_ID");
				header.assertColumnExists("SNP_GENE_IDS");
				header.assertColumnExists("UPSTREAM_GENE_DISTANCE");
				header.assertColumnExists("DOWNSTREAM_GENE_DISTANCE");
				header.assertColumnExists("CHR_ID");
				header.assertColumnExists("CHR_POS");
				try(PrintWriter w = super.openPathOrStdoutAsPrintWriter(this.output)) {
					if(with_header) {
						w.println("#contig\tstart\tend\t"+line);
						}
					while((line=br.readLine())!=null) {
						final FileHeader.RowMap row= header.toMap(line);
						final String srcContig = ctgConverter.apply(row.getOrDefault("CHR_ID", ""));
						
						final String contig;
						int chromStart;
						int chromEnd;
						final String upstream_gene_id = removeVersion(row.getOrDefault("UPSTREAM_GENE_ID",""));
						final String downtream_gene_id = removeVersion(row.getOrDefault("DOWNSTREAM_GENE_ID",""));
						final List<Interval> snp_genes = Arrays.stream(CharSplitter.COMMA.split(row.getOrDefault("SNP_GENE_IDS",""))).
								map(S->removeVersion(S))
								.filter(S->!StringUtils.isBlank(S))
								.map(S->name2gene.get(S))
								.filter(R->R!=null)
								.collect(Collectors.toList());
						if(!StringUtils.isBlank(upstream_gene_id) && !StringUtils.isBlank(downtream_gene_id) && snp_genes.isEmpty()) {
							final Interval rL = name2gene.get(upstream_gene_id);
							if(rL==null) {
								LOG.warn("cannot find upstream gene for "+upstream_gene_id+" in "+row);
								continue;
								}
							final int xL =Integer.parseInt(StringUtils.ifBlank(row.getOrDefault("UPSTREAM_GENE_DISTANCE",""),"0"));
							
							
							final Interval rR = name2gene.get(downtream_gene_id);
							if(rR==null) {
								LOG.warn("cannot find downstream gene for "+downtream_gene_id+" in "+row);
								continue;
								}
							if(!rL.contigsMatch(rR)) {
								LOG.warn("contig mismatch in "+rL+" vs "+rR+" in "+row);
								continue;
								}
							final int xR =Integer.parseInt(StringUtils.ifBlank(row.getOrDefault("DOWNSTREAM_GENE_DISTANCE",""),"0"));

							contig = rL.getContig();
							chromStart = rL.getStart() + xL;
							chromEnd = rR.getEnd() - xR;
							if(chromStart> chromEnd) {
								LOG.warn("chromStart > chromEnd in "+row);
								continue;
								}
							}
						else if(!snp_genes.isEmpty()) {
							if(!StringUtils.isBlank(srcContig) &&  snp_genes.stream().map(R->R.getContig()).collect(Collectors.toSet()).size()!=1) {
								snp_genes.stream()
									.filter(R->!R.getContig().equals(srcContig))
									.forEach(V->LOG.info("ignoring "+V+" because contig does not match "+srcContig+" for "+row));
								snp_genes.removeIf(R->!R.getContig().equals(srcContig));
								}
							
							if(snp_genes.stream().map(R->R.getContig()).collect(Collectors.toSet()).size()!=1) {
								LOG.warn("multiple contigs found for "+row +" "+snp_genes);
								continue;
								}
							if(snp_genes.isEmpty()) {
								LOG.warn("no gene on "+srcContig+" was found for "+row);
								continue;
								}
							contig = snp_genes.get(0).getContig();
							chromStart = snp_genes.stream().mapToInt(R->R.getStart()).min().getAsInt();
							chromEnd = snp_genes.stream().mapToInt(R->R.getEnd()).max().getAsInt();
							}
						else if(!StringUtils.isBlank(srcContig) &&
								!StringUtils.isBlank(row.getOrDefault("CHR_POS", ""))) {
							final List<Integer> positions = Arrays.stream(CharSplitter.SEMICOLON.split(row.getOrDefault("CHR_POS", "")))
									.map(S->S.trim())
									.filter(S->!StringUtils.isBlank(S))
									.map(S->Integer.parseInt(S))
									.collect(Collectors.toList());
							if(positions.isEmpty()) {
								LOG.warn("cannot find positions in "+row);
								continue;
								}
							contig = srcContig;
							chromStart = positions.stream().mapToInt(R->R.intValue()).min().getAsInt();
							chromEnd = positions.stream().mapToInt(R->R.intValue()).max().getAsInt();
							if(!this.enable_pos) {
								LOG.info("got position but --position is not set");
								}
							else
								{
								if(!warning_position_printed) {
									LOG.warn("Printing from position. Check the build ! "+row);
									warning_position_printed= true;
									}
								}
							}
						else
							{
							LOG.warn("no gene was found for "+row.toString());
							continue;
							}
						
						w.print(contig);
						w.print("\t");
						w.print(chromStart-1/*bed*/);
						w.print("\t");
						w.print(chromEnd);
						w.print("\t");
						w.print(line);
						w.println();
						}
					w.flush();
					}
				}
			
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
		}
		}
	
	public static void main(String[] args) {
		new GwasCatalogToBed().instanceMainWithExit(args);

	}

}
