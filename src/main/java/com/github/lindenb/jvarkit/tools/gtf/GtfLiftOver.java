/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;

/**

BEGIN_DOC

```
$ gunzip -c ~/data/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz | grep ^X | head 
X	Regulatory_Build	CTCF_binding_site	100072201	100072600	.	.	.	ID=CTCF_binding_site:ENSR00000911312;bound_end=100072600;bound_start=100072201;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100076001	100076400	.	.	.	ID=CTCF_binding_site:ENSR00000911313;bound_end=100076400;bound_start=100076001;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100080601	100080800	.	.	.	ID=CTCF_binding_site:ENSR00000911314;bound_end=100080800;bound_start=100080601;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100099201	100099400	.	.	.	ID=CTCF_binding_site:ENSR00000911321;bound_end=100099400;bound_start=100099201;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	10015601	10016000	.	.	.	ID=CTCF_binding_site:ENSR00000901123;bound_end=10016000;bound_start=10015601;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	10018001	10018400	.	.	.	ID=CTCF_binding_site:ENSR00000901124;bound_end=10018400;bound_start=10018001;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100185401	100185800	.	.	.	ID=CTCF_binding_site:ENSR00000911326;bound_end=100185800;bound_start=100185401;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	10021601	10021800	.	.	.	ID=CTCF_binding_site:ENSR00000901125;bound_end=10021800;bound_start=10021601;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100351201	100351400	.	.	.	ID=CTCF_binding_site:ENSR00001161915;bound_end=100351400;bound_start=100351201;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100368001	100368800	.	.	.	ID=CTCF_binding_site:ENSR00000342428;bound_end=100368800;bound_start=100368001;description=CTCF binding site;feature_type=CTCF Binding Site

$ java -jar dist/gtfliftover.jar -x jeter.boum -f ../htsjdk/src/test/resources/htsjdk/samtools/liftover/hg18ToHg19.over.chain ~/data/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz  2> /dev/null | head
chrX	Regulatory_Build	CTCF_binding_site	100185545	100185944	.	.	.	ID=CTCF_binding_site:ENSR00000911312;bound_end=100072600;bound_start=100072201;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100189345	100189744	.	.	.	ID=CTCF_binding_site:ENSR00000911313;bound_end=100076400;bound_start=100076001;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100193945	100194144	.	.	.	ID=CTCF_binding_site:ENSR00000911314;bound_end=100080800;bound_start=100080601;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100212545	100212744	.	.	.	ID=CTCF_binding_site:ENSR00000911321;bound_end=100099400;bound_start=100099201;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	10055601	10056000	.	.	.	ID=CTCF_binding_site:ENSR00000901123;bound_end=10016000;bound_start=10015601;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	10058001	10058400	.	.	.	ID=CTCF_binding_site:ENSR00000901124;bound_end=10018400;bound_start=10018001;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100298745	100299144	.	.	.	ID=CTCF_binding_site:ENSR00000911326;bound_end=100185800;bound_start=100185401;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	10061601	10061800	.	.	.	ID=CTCF_binding_site:ENSR00000901125;bound_end=10021800;bound_start=10021601;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100464545	100464744	.	.	.	ID=CTCF_binding_site:ENSR00001161915;bound_end=100351400;bound_start=100351201;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100481345	100482144	.	.	.	ID=CTCF_binding_site:ENSR00000342428;bound_end=100368800;bound_start=100368001;description=CTCF binding site;feature_type=CTCF Binding Site
```

when the GTF contains a gene, all sub-section MUST be re-mapped.

```
$ gunzip -c ~/src/jvarkit-git/src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | head
22	ensembl_havana	exon	41697526	41697776	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; exon_id "ENSE00001307363"; exon_version "4"; tag "basic";
22	ensembl_havana	five_prime_utr	41697526	41697776	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
22	ensembl_havana	gene	41697526	41756151	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
22	ensembl_havana	transcript	41697526	41756151	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
22	havana	exon	41697719	41697776	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000486331"; transcript_version "1"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-002"; transcript_source "havana"; transcript_biotype "retained_intron"; havana_transcript "OTTHUMT00000320697"; havana_transcript_version "1"; exon_id "ENSE00001942555"; exon_version "1";
22	havana	transcript	41697719	41732847	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000486331"; transcript_version "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-002"; transcript_source "havana"; transcript_biotype "retained_intron"; havana_transcript "OTTHUMT00000320697"; havana_transcript_version "1";
22	ensembl_havana	exon	41716659	41716717	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; exon_number "2"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; exon_id "ENSE00003628979"; exon_version "1"; tag "basic";
22	ensembl_havana	five_prime_utr	41716659	41716664	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
22	havana	exon	41716659	41716717	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000486331"; transcript_version "1"; exon_number "2"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-002"; transcript_source "havana"; transcript_biotype "retained_intron"; havana_transcript "OTTHUMT00000320697"; havana_transcript_version "1"; exon_id "ENSE00003530265"; exon_version "1";
22	ensembl	CDS	41716665	41716717	.	+	0	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000351589"; transcript_version "4"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; protein_id "ENSP00000263243"; protein_version "5"; tag "basic";

$ java -jar dist/gtfliftover.jar -f ../htsjdk/src/test/resources/htsjdk/samtools/liftover/hg18ToHg19.over.chain ~/src/jvarkit-git/src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz  --header 2> /dev/null  | head
chr22	ensembl_havana	exon	43367582	43367832	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; exon_id "ENSE00001307363"; exon_version "4"; tag "basic";
chr22	ensembl_havana	five_prime_utr	43367582	43367832	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
chr22	ensembl_havana	gene	43367582	43426207	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
chr22	ensembl_havana	transcript	43367582	43426207	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
chr22	havana	exon	43367775	43367832	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000486331"; transcript_version "1"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-002"; transcript_source "havana"; transcript_biotype "retained_intron"; havana_transcript "OTTHUMT00000320697"; havana_transcript_version "1"; exon_id "ENSE00001942555"; exon_version "1";
```

END_DOC
*/
@Program(
		name="gtfliftover",
		description="LiftOver GTF file.",
		creationDate="20190823",
		modificationDate="20190823",
		keywords= {"gtf","liftover"}
		)
public class GtfLiftOver
	extends Launcher
	{
	private static final Logger LOG = Logger.build(GtfLiftOver.class).make();
	private final GTFCodec gtfCodec = new GTFCodec();
	private LiftOver liftOver = null;
	private ContigNameConverter contigNameConverter = null;
	
	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-f","--chain"},description="LiftOver file.",required=true)
	private File liftOverFile = null;
	@Parameter(names={"-x","--failed"},description="write  failing the liftOver here. Optional.")
	private Path failedFile = null;
	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;
	@Parameter(names={"-H","--header"},description= "include gtf header")
	private boolean include_header = false;

	public GtfLiftOver()
		{
		}

	private Optional<Interval> doliftOver(final Locatable gtfLine) {
		final String ctg = this.contigNameConverter.apply(gtfLine.getContig());
		if(StringUtils.isBlank(ctg)) return Optional.empty();
		return Optional.ofNullable(this.liftOver.liftOver(new Interval(ctg,gtfLine.getStart(),gtfLine.getEnd()),this.userMinMatch));
		}
	
	private String applyLift(final Locatable loc,final String line) {
		final String tokens[] = CharSplitter.TAB.split(line,6);
		tokens[0]=loc.getContig();
		tokens[3]=String.valueOf(loc.getStart());
		tokens[4]=String.valueOf(loc.getEnd());
		return String.join("\t", tokens);
		}
	
	@Override
	public int doWork(final List<String> args)  {
		BufferedReader in = null;
		boolean got_gtf_record = false;
		PrintWriter p = null;
		PrintWriter fail = null;
		try {
			fail = this.failedFile==null?
					new PrintWriter(new NullOuputStream()):
					super.openPathOrStdoutAsPrintWriter(this.failedFile)
					;
			p = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			
			this.liftOver = new LiftOver(this.liftOverFile);
			this.contigNameConverter =  ContigNameConverter.fromContigSet(this.liftOver.getContigMap().keySet());
			
			in = openBufferedReader(oneFileOrNull(args));
			final Map<String,List<String>> gene2lines = new HashMap<>(50_000);
			for(;;)
				{
				final String line = in.readLine();
				if(line==null) break;
				if(StringUtils.isBlank(line)) continue;
				if(line.startsWith("#")) {
					if(!got_gtf_record && include_header) {
						fail.println(line);
						p.println(line);
						}
					continue;
					}
				final GTFLine gtfLine = this.gtfCodec.decode(line);
				if(gtfLine==null) {
					fail.println(line);
					continue;
					}
				got_gtf_record = true;
				
				final String gene_key = gtfLine.getAttribute("gene_id");
				if(StringUtils.isBlank(gene_key))
					{
					Optional<Interval> lifted = doliftOver(gtfLine);
					if(lifted.isPresent()) {
						p.println(applyLift(lifted.get(),line));
						}
					else
						{
						fail.println(line);
						}
					}
				else
					{
					List<String> lines = gene2lines.get(gene_key);
					if(lines==null) {
						lines = new ArrayList<>();
						gene2lines.put(gene_key,lines);
						}
					lines.add(line);
					}
				}
			
			for(final List<String> lines:gene2lines.values()) {
			final List<Optional<Interval>> intervals = lines.stream().map(L->doliftOver(gtfCodec.decode(L))).filter(R->R.isPresent()).collect(Collectors.toList());
			/* all have been mapped */
			if(intervals.size()!=lines.size() || 
				/* all mapped to same contig */
				intervals.stream().map(R->R.get().getContig()).collect(Collectors.toSet()).size()!=1) {
				for(final String L: lines) {
					fail.println(L);
					}
				}
			else
				{
				for(int i=0;i< intervals.size();i++) {
					p.println(applyLift(intervals.get(i).get(),lines.get(i)));
					}
				}
			}
			
			fail.flush();
			fail.close();
			fail = null;
			p.flush();
			p.close();
			p = null;
			return RETURN_OK;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(fail);
			CloserUtil.close(p);
			}
		}
	
	public static void main(final String[] args)
		{
		new GtfLiftOver().instanceMainWithExit(args);
		}
	}
