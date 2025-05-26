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
package com.github.lindenb.jvarkit.tools.kg2bed;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;
import java.util.function.Predicate;

import org.apache.commons.jexl2.JexlContext;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.jexl.JexlPredicate;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptReader;
import com.beust.jcommander.Parameter;


/**

BEGIN_DOC

##  See also

**ucsc** tools :  genePredToBed

### Example

```
$ curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" |\
  gunzip -c |\
  java -jar dist/jvarkit.jar kg2bed
chr1	11873	14409	+	uc001aaa.3	TRANSCRIPT	uc001aaa.3
chr1	11873	12227	+	uc001aaa.3	EXON	Exon 1
chr1	12227	12612	+	uc001aaa.3	INTRON	Intron 1
chr1	11873	12227	+	uc001aaa.3	UTR	UTR3
chr1	12612	12721	+	uc001aaa.3	EXON	Exon 2
chr1	12721	13220	+	uc001aaa.3	INTRON	Intron 2
chr1	12612	12721	+	uc001aaa.3	UTR	UTR3
chr1	13220	14409	+	uc001aaa.3	EXON	Exon 3
chr1	13220	14409	+	uc001aaa.3	UTR	UTR3
chr1	11873	14409	+	uc010nxr.1	TRANSCRIPT	uc010nxr.1
chr1	11873	12227	+	uc010nxr.1	EXON	Exon 1
chr1	12227	12645	+	uc010nxr.1	INTRON	Intron 1
chr1	11873	12227	+	uc010nxr.1	UTR	UTR3
chr1	12645	12697	+	uc010nxr.1	EXON	Exon 2
chr1	12697	13220	+	uc010nxr.1	INTRON	Intron 2
```

END_DOC
*/


@Program(name="kg2bed",
	description="converts UCSC knownGenes file to BED.",
	keywords={"ucsc","bed","knownGenes"},
	biostars= {151628,9557497},
	creationDate = "20140311",
	modificationDate="20230815",
	jvarkit_amalgamion = true
	)
public class KnownGenesToBed extends Launcher
	{
	private static final Logger LOG = Logger.of(KnownGenesToBed.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--exclude","--hide"},description="don't show the following items (comma separated, one of 'INTRON,UTR,CDS,EXON,TRANSCRIPT,NON_CODING,CODING'). Empty don't hide anything")
	private String hideStr="";
	@Parameter(names={"-s","--select"},description="JEXL select expression. Object 'kg' is an instance of KnownGene (https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/ucsc/KnownGene.java)." +JexlPredicate.OPT_WHAT_IS_JEXL)
	private String selectExpr="";
	@Parameter(names={"-sql","--sql"},description= "SQL Schema URI. "+UcscTranscriptReader.SQL_DESC)
	private String sqlUri="";

	
	
	private void print(	final PrintWriter out, final UcscTranscript kg,final int start,final int end,final String type,final String name)
		{
		if(start>=end) return;
		out.print(kg.getContig());
		out.print('\t');
		out.print(start);
		out.print('\t');
		out.print(end);
		out.print('\t');
		out.print(kg.isPositiveStrand()?'+':'-');
		out.print('\t');
		out.print(kg.getTranscriptId());
		out.print('\t');
		out.print(type);
		out.print('\t');
		out.print(name);
		out.println();
		}
	// 
	private void scan(final PrintWriter out,final UcscTranscriptCodec codec,final BufferedReader r) throws IOException
		{
		boolean hide_introns = false;
		boolean hide_utrs = false;
		boolean hide_cds = false;
		boolean hide_exons = false;
		boolean hide_transcripts = false;
		boolean hide_non_coding = false;
		boolean hide_coding = false;
		final Predicate<JexlContext> predicate = StringUtils.isBlank(selectExpr)?
				KG->true:
				new JexlPredicate(this.selectExpr)
				; 

		for(String str : CharSplitter.COMMA.split(this.hideStr)) {
			if(StringUtils.isBlank(str)) continue;
			str = str.trim().toUpperCase();
			if(str.equals("INSTRON") || str.equals("INSTRONS")) hide_introns=true;
			if(str.equals("UTR") || str.equals("UTRs")) hide_utrs=true;
			if(str.equals("CDS")) hide_cds=true;
			if(str.equals("EXON") || str.equals("EXONS")) hide_exons=true;
			if(str.equals("TRANSCRIPT") || str.equals("TRANSCRIPTS")) hide_transcripts=true;
			if(str.equals("NON_CODING")) hide_non_coding=true;
			if(str.equals("CODING")) hide_coding=true;
			}

		String line;
		while((line=r.readLine())!=null)
			{
			final UcscTranscript kg= codec.decode(line);
			if(kg==null) continue;

			if(hide_coding && kg.isProteinCoding()) continue;
			if(hide_non_coding && !kg.isProteinCoding()) continue;

			if(!predicate.test(kg.asJexlContext())) continue;

			if(!hide_transcripts) {
				print(out, kg,kg.getBedStart(),kg.getBedEnd(),"TRANSCRIPT",kg.getTranscriptId());
				}
			if(!hide_exons) {
				for(final UcscTranscript.Exon exon: kg.getExons())
					{
					print(out, kg,exon.getBedStart(),exon.getBedEnd(),"EXON",exon.getName());
					}
				}
			if(!hide_introns) {
				for(final UcscTranscript.Intron intron: kg.getIntrons())
					{
					print(out, kg,intron.getBedStart(),intron.getBedEnd(),"INTRON",intron.getName());
					}
				}

			if(!hide_cds && kg.isProteinCoding()) {
				for(final UcscTranscript.CDS cds: kg.getCDS())
					{
					print(out, kg,cds.getBedStart(),cds.getBedEnd(),"CDS",cds.getName());
					}
				}
			if(!hide_utrs && kg.isProteinCoding()) {
				for(final UcscTranscript.UTR utr: kg.getUTRs())
					{
					print(out, kg,utr.getBedStart(),utr.getBedEnd(),"UTR",utr.getName());
					}
				}
			}
		}

	@Override
	public int doWork(final List<String> args) {
		try
			{
			final UcscTranscriptCodec codec = StringUtils.isBlank(this.sqlUri)? new UcscTranscriptCodec():new UcscTranscriptCodec(sqlUri);

			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				if(args.isEmpty())
					{
					try(BufferedReader r =IOUtils.openStreamForBufferedReader(super.stdin())) {
						scan(out,codec,r);
						}
					}
				else
					{
					for(final String filename:args)
						{
						try(BufferedReader r=IOUtils.openURIForBufferedReading(filename)) {
							scan(out,codec,r);
							}
						}
					}
				out.flush();
				}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}

	
	public static void main(final String[] args) {
	new KnownGenesToBed().instanceMainWithExit(args);

	}

}
