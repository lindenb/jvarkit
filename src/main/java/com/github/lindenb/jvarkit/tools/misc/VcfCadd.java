/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter.OnNotFound;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
BEGIN_DOC
## Example

```bash
$ java -Dhttp.proxyHost=my.proxy.host.fr -Dhttp.proxyPort=1234 -jar dist/vcfcadd.jar \
	-u "http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz"  \
	src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz 2> /dev/null | ~/package/bcftools/bcftools annotate -x '^INFO/CADD_SCORE,INFO/CADD_PHRED'

##fileformat=VCFv4.2
(...)
##INFO=<ID=CADD_PHRED,Number=A,Type=Float,Description="PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc.  URI was http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz">
##INFO=<ID=CADD_SCORE,Number=A,Type=Float,Description="Score suggests that that variant is likely to be  observed (negative values) vs simulated(positive values).However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more likely to have deleterious effects. URI was http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz">
##VcfCaddCmdLine=-u http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	905606	rs540662886	G	C,A	41743.9	PASS	CADD_PHRED=3.426,.;CADD_SCORE=0.082875,.
(...)
1	905621	rs368876607	G	A	14291.5	PASS	CADD_PHRED=6.025;CADD_SCORE=0.334762
(...)
1	905669	rs111483874	C	G,T	86574.3	PASS	CADD_PHRED=12.77,.;CADD_SCORE=1.39614,.
(...)
1	905723	rs150703609	G	A	15622.1	PASS	CADD_PHRED=23.7;CADD_SCORE=4.05532
1	905726	rs751084833	C	T,A	8733.36	PASS	.
1	905727	rs761609807	G	A	12936.9	PASS	.
(..)
```

## History

  * 2018-04-25 : changing INFO -type to 'A', splitting into two CADD_score/phred and adding dict converter


END_DOC
*/
@Program(
	name="vcfcadd",
	description= "Annotate VCF with  Combined Annotation Dependent Depletion (CADD) (Kircher & al. "+
		"A general framework for estimating the relative pathogenicity of human genetic variants. "+
		"Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892." +
		"PubMed PMID: 24487276.",
	keywords={"vcf","prediction","cadd","annotation"}
	)

public class VcfCadd extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfCadd.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	public static final String DEFAULT_URI="http://krishna.gs.washington.edu/download/CADD/v1.2/whole_genome_SNVs.tsv.gz";
	private TabixFileReader tabix=null;
	@Parameter(names={"-u","--uri","--tabix"},description="Combined Annotation Dependent Depletion (CADD) Tabix file URI ")
	private String ccaduri=DEFAULT_URI;
	@Parameter(names={"-S","--score","--score-tag"},description="INFO tag for score")
	private String CADD_FLAG_SCORE = "CADD_SCORE";
	@Parameter(names={"-P","--phred","--phred-tag"},description="INFO tag for phred")
	private String CADD_FLAG_PHRED = "CADD_PHRED";
	@Parameter(names={"-d","--buffer-size"},description="Buffer size / processing window size")
	private int buffer_distance=1000;

	
	private final Pattern TAB=Pattern.compile("[\t]");
	private ContigNameConverter convertToCaddContigs = null;
	
	/**
	$ wget -q -O - "http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz" | gunzip  -c | head 

## CADD v1.3 (c) University of Washington and Hudson-Alpha Institute for Biotechnology 2013-2015. All rights reserved.
#Chrom	Pos	Ref	Alt	RawScore	PHRED
1	10177	A	AC	-0.402095	0.382
1	10235	T	TA	-0.196330	1.117
1	10352	T	TA	-0.240617	0.900
1	10505	A	T	0.618129	8.281
1	10506	C	G	0.789734	9.401
1	10511	G	A	0.772599	9.294
1	10539	C	A	1.054015	10.96
1	10542	C	T	1.113857	11.29

	
	*/
	private static class Record
		{
		final int pos;
		final Allele ref;
		final Allele alt;
		final float score;
		final float phred;
		Record(final String tokens[]) {
			if(tokens.length!=6) throw new JvarkitException.TokenErrors("Bad CADD line . Expected 6 fields",tokens);
			this.pos= Integer.parseInt(tokens[1]);
			this.ref = Allele.create(tokens[2],true);
			this.alt = Allele.create(tokens[3],false);
			this.score = Float.parseFloat(tokens[4]);
			this.phred = Float.parseFloat(tokens[5]);
			}
		}
	
	public VcfCadd()
		{
		}
	
	private final Set<String> contigsNotFounds = new HashSet<>();
	
	private void runTabix(final List<VariantContext> buffer)
			throws IOException
		{
		if(buffer.isEmpty()) return;
		final String contigVcf = buffer.get(0).getContig();
		final String chromCadd = this.convertToCaddContigs.apply(contigVcf);
		if(StringUtil.isBlank(chromCadd)) {
			if(!contigsNotFounds.contains(contigVcf) && this.contigsNotFounds.size()<10000 ) {
				contigsNotFounds.add(contigVcf);
				LOG.warning("Cannot find or convert VCF contig "+contigVcf+" to cadd contigs");
				}
			return;
			}
		
		
		int start=Integer.MAX_VALUE;
		int end=0;
		for(final VariantContext ctx:buffer)
			{
			start= Math.min(start, ctx.getStart());
			end= Math.max(end, ctx.getEnd());
			}
		LOG.info("Scanning "+contigVcf+":"+start+"-"+end);
		final Map<ContigPosRef,List<Record>> caddMap = new HashMap<>(buffer.size());
		
		
		Iterator<String> iter = this.tabix.iterator(
				chromCadd,
				(int)Math.max(1,start),
				(int)Math.min(Integer.MAX_VALUE,end)
				);
		while(iter.hasNext() )
			{
			final String line=iter.next();
			if(line.startsWith("#") || StringUtil.isBlank(line))continue;
			final String tokens[] = TAB.split(line);
			final Record rec=new Record(tokens);
			if(!tokens[0].equals(chromCadd)) throw new IllegalStateException("Expected CADD contig "+chromCadd+" in "+line);
			
			final ContigPosRef cpr = new ContigPosRef(contigVcf,rec.pos,rec.ref);
			List<Record> L = caddMap.get(cpr);
			if(L==null) {
				L = new ArrayList<>();
				caddMap.put(cpr,L);
				}		
			L.add(rec);
			}
		CloserUtil.close(iter);
		
		
		for(int i=0;i< buffer.size();++i)
			{
			final VariantContext ctx=buffer.get(i);
			final List<Record> cadd_rec_for_ctx = caddMap.get(new ContigPosRef(ctx));
			if(cadd_rec_for_ctx==null || cadd_rec_for_ctx.isEmpty()) continue;
			
			final List<Float> cadd_array_score=new ArrayList<>();
			final List<Float> cadd_array_phred=new ArrayList<>();
			boolean got_non_null = false;
			for(final Allele alt:ctx.getAlternateAlleles())
				{
				final Record rec = cadd_rec_for_ctx.
						stream().
						filter(REC->REC.alt.equals(alt)).
						findAny().
						orElse(null);
				if(rec==null) {
					cadd_array_score.add(null);
					cadd_array_phred.add(null);
					}
				else {
					got_non_null = true;
					cadd_array_score.add(rec.score);
					cadd_array_phred.add(rec.phred);
					}
				}
				
			
			if(!cadd_array_score.isEmpty() && got_non_null)
				{
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(this.CADD_FLAG_SCORE, cadd_array_score);
				vcb.attribute(this.CADD_FLAG_PHRED, cadd_array_phred);
				buffer.set(i, vcb.make());
				}
			}

		
		}
	@Override
	protected int doVcfToVcf(final String inputName, VcfIterator in, VariantContextWriter out) {
			
			try {
			final VCFHeader header=in.getHeader();
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header).logger(LOG);
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			if(header.getInfoHeaderLine(this.CADD_FLAG_PHRED)!=null) {
				throw new JvarkitException.DuplicateVcfHeaderInfo(header, this.CADD_FLAG_PHRED);
			}
			if(header.getInfoHeaderLine(this.CADD_FLAG_SCORE)!=null) {
				throw new JvarkitException.DuplicateVcfHeaderInfo(header, this.CADD_FLAG_SCORE);
			}

			
			header.addMetaDataLine(new VCFInfoHeaderLine(
					this.CADD_FLAG_SCORE,
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"Score suggests that that variant is likely to be  observed (negative values) vs simulated(positive values)."+
					"However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more likely to have deleterious effects." +
					" URI was " +this.ccaduri
					));
			header.addMetaDataLine(new VCFInfoHeaderLine(
					this.CADD_FLAG_PHRED,
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc. " +
					" URI was " +this.ccaduri
					));
			
			out.writeHeader(header);
			final List<VariantContext> buffer= new ArrayList<>();
			for(;;)
				{	
				VariantContext ctx=null;
				if(in.hasNext())
					{
					ctx=progress.watch(in.next());
					}
				if( ctx==null ||
					(!buffer.isEmpty() &&
					(buffer.get(0).getContig().equals(ctx.getContig())) && (ctx.getEnd()-buffer.get(0).getStart())>buffer_distance))
					{
					if(!buffer.isEmpty())
						{
						runTabix(buffer);
						for(final VariantContext c:buffer)
							{
							out.add(c);
							}
						}
					if(ctx==null) break;
					buffer.clear();
					}
				buffer.add(ctx);
				}
			progress.finish();
			return 0;
			} 
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.CADD_FLAG_PHRED.equals(CADD_FLAG_SCORE)) {
				LOG.error("tag phred same as tag score");
				return -1;
				}
			if(this.ccaduri==null || !this.ccaduri.endsWith(".gz"))
				{
				LOG.error("CCAD uri should end with gz. got "+this.ccaduri);
				return -1;
				}
			
			LOG.info("Loading index for "+this.ccaduri+". Please wait...");
			this.tabix = new TabixFileReader(this.ccaduri);
			this.convertToCaddContigs = ContigNameConverter.fromContigSet(this.tabix.getChromosomes());
			this.convertToCaddContigs.setOnNotFound(OnNotFound.SKIP);
			LOG.info("End loading index");
			
			return doVcfToVcf(args,outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(tabix);
			tabix=null;
			}
		}

	public static void main(final String[] args)
		{
		new VcfCadd().instanceMainWithExit(args);
		}
	}
