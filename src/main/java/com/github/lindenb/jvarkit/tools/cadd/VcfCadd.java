/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.cadd;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.variant.vcf.VCFIterator;

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
## Note to self

I got problem with the certificate. Fixed with `-Dcom.sun.security.enableAIAcaIssuers=true -Dcom.sun.net.ssl.checkRevocation=false `

END_DOC
*/
@Program(
	name="vcfcadd",
	description= "Annotate VCF with  Combined Annotation Dependent Depletion (CADD) (Kircher & al. "+
		"A general framework for estimating the relative pathogenicity of human genetic variants. "+
		"Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892." +
		"PubMed PMID: 24487276.",
	creationDate="20140218",
	modificationDate="20220119",
	keywords={"vcf","prediction","cadd","annotation"}
	)

public class VcfCadd extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VcfCadd.class).make();

	/** global can be used by vcf2r for Matilde */
	public static final String DEFAULT_CADD_FLAG_SCORE = "CADD_SCORE";
	/** global can be used by vcf2r for Matilde */
	public static final String DEFAULT_CADD_FLAG_PHRED = "CADD_PHRED";
	
	public static final String DEFAULT_URI="https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz";
	private TabixFileReader tabix=null;
	@Parameter(names={"-u","--uri","--tabix","--cadd"},description="Combined Annotation Dependent Depletion (CADD) Tabix file URI ")
	private String ccaduri=DEFAULT_URI;
	@Parameter(names={"-S","--score","--score-tag"},description="INFO tag for score")
	private String CADD_FLAG_SCORE = DEFAULT_CADD_FLAG_SCORE;
	@Parameter(names={"-P","--phred","--phred-tag"},description="INFO tag for phred")
	private String CADD_FLAG_PHRED = DEFAULT_CADD_FLAG_PHRED;
	@Parameter(names={"-d","--buffer-size"},description="Buffer size / processing window size. " + DistanceParser.OPT_DESCRIPTION , converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int buffer_distance= 1_000;
	@Parameter(names={"-f","--fields"},description=
			"Other Fields to be included. See the header of http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs_inclAnno.tsv.gz . Multiple separeted by space, semicolon or comma."
			+ " Warning: This tool currently uses the first CHROM/POS/REF/ALT values it finds while I saw some duplicated fields in 'whole_genome_SNVs_inclAnno.tsv.gz'.")
	private String otherFieldsStr = "";
	@Parameter(names={"-na"},description="value  used for 'allele-not-found'.")
	private float NA_value = -999f;

	
	private final CharSplitter TAB = CharSplitter.TAB;
	private ContigNameConverter convertToCaddContigs = null;
	private List<String> headerLabel = null;
	private Map<String,Integer> headerLabel2column = null;
	private int column_index_for_Alt = -1;
	private int column_index_for_RawScore = -1;
	private int column_index_for_PHRED = -1;
	private final Set<String> userFields = new HashSet<>();
	private final Map<String,VCFHeaderLineCount> annoField2linecount = new HashMap<>();
	
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
	private class Record implements Locatable {
		final String contig;
		final int pos;
		final Allele ref;
		final Allele alt;
		final float score;
		final float phred;
		final Map<String,String> otherKeyValues;
		Record(final String contig,final List<String> tokens) {
			if(tokens.size()<6) throw new JvarkitException.TokenErrors("Bad CADD line . Expected at least 6 fields",tokens);
			this.contig = contig;
			this.pos= Integer.parseInt(tokens.get(1));
			this.ref = Allele.create(tokens.get(2),true);
			this.alt = Allele.create(tokens.get(VcfCadd.this.column_index_for_Alt),false);
			this.score = Float.parseFloat(tokens.get(VcfCadd.this.column_index_for_RawScore));
			this.phred = Float.parseFloat(tokens.get(VcfCadd.this.column_index_for_PHRED));
			if(userFields.isEmpty())
				{
				this.otherKeyValues  = null;
				}
			else
				{
				otherKeyValues= new HashMap<>(userFields.size());
				for(final String uf: userFields) {
					Integer i = headerLabel2column.get(uf);
					if(i==null || i.intValue()>=tokens.size()) continue;
					final String v =  tokens.get(headerLabel2column.get(uf));
					this.otherKeyValues.put(uf, v);
					}
				}
			}
		@Override
		public String getContig()
			{
			return this.contig;
			}
		@Override
		public int getStart()
			{
			return this.pos;
			}
		@Override
		public int getEnd()
			{
			return this.pos;
			}
		}
	
	public VcfCadd()
		{
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private Locatable lastInterval = null;
	private final List<Record> buffer = new ArrayList<>();
	
	private List<Record> query(final VariantContext loc) {
		final String chromCadd = this.convertToCaddContigs.apply(loc.getContig());
		if(StringUtils.isBlank(chromCadd)) return Collections.emptyList();
		final Locatable query = new SimpleInterval(chromCadd,loc.getStart(),loc.getEnd());
		if(this.lastInterval==null || !this.lastInterval.contains(query)) {
			this.buffer.clear();
			this.lastInterval = new SimpleInterval(chromCadd, Math.max(0, loc.getStart()-1), Math.max(loc.getEnd(), loc.getStart()+this.buffer_distance));
			Iterator<String> iter = this.tabix.iterator(
					this.lastInterval.getContig(),
					this.lastInterval.getStart(),
					this.lastInterval.getEnd()
					);
			while(iter.hasNext()) {
				final String line=iter.next();
				if(line.startsWith("#") || StringUtil.isBlank(line))continue;
				final List<String> tokens = TAB.splitAsStringList(line);
			
				if(!AcidNucleics.isATGC(tokens.get(2)))
					{
					LOG.warn("REF allele not suitable in  line "+line+". skipping");
					continue;
					}
				if(!AcidNucleics.isATGC(tokens.get(this.column_index_for_Alt)))
					{
					LOG.warn("ALT allele not suitable in  line "+line+". skipping");
					continue;
					}
			
				if(!tokens.get(0).equals(chromCadd)) throw new IllegalStateException("Expected CADD contig "+chromCadd+" in "+line);

				final Record rec=new Record(loc.getContig(),tokens);

				this.buffer.add(rec);
				}
			}
		return this.buffer.stream().filter(R->
				R.getStart() == loc.getStart() &&
				R.ref.equals(loc.getReference()) &&
				loc.getAlternateAlleles().contains(R.alt)
				).collect(Collectors.toList());
		}

	private void  checkRange(float f) {
		if(f <= this.NA_value || f >= this.NA_value) {
			throw new IllegalArgumentException("Got score = "+f +" but the value for N/A is "+this.NA_value);
			}
		}

	private VariantContext runTabix(final VariantContext ctx) throws IOException {
			final List<Record> cadd_rec_for_ctx = query(ctx);
			if(cadd_rec_for_ctx.isEmpty()) return ctx;
			
			final List<Float> cadd_array_score=new ArrayList<>();
			final List<Float> cadd_array_phred=new ArrayList<>();
			boolean got_non_null = false;
			for(final Allele alt:ctx.getAlternateAlleles()) {
				final Record rec = cadd_rec_for_ctx.
						stream().
						filter(REC->REC.alt.equals(alt)).
						findAny().
						orElse(null);
				if(rec==null) {
					cadd_array_score.add(NA_value);
					cadd_array_phred.add(NA_value);
					}
				else {
					got_non_null = true;
					checkRange(rec.score);
					checkRange(rec.phred);
					cadd_array_score.add(rec.score);
					cadd_array_phred.add(rec.phred);
					}
				}

			final Map<String,List<String>> cadd_array_other = new HashMap<>();

			if(got_non_null && !this.userFields.isEmpty())
				{
				for(final String key: this.userFields) {
					
					if(getFieldHeaderLineCount(key).equals(VCFHeaderLineCount.A)) {
						final List<String> vals = new ArrayList<>();
						cadd_array_other.put(key, vals);
						for(final Allele alt:ctx.getAlternateAlleles())
							{
							final String rec = cadd_rec_for_ctx.
									stream().
									filter(REC->REC.alt.equals(alt)).
									filter(REC->REC.otherKeyValues!=null).
									filter(REC->REC.otherKeyValues.containsKey(key)).
									map(REC->REC.otherKeyValues.get(key)).
									findAny().
									orElse(null);
							if(StringUtil.isBlank(rec) || rec.equals("NA") ) {
								vals.add(null);
								}
							else
								{
								vals.add(rec);
								}
							}
						}
					else
						{
						cadd_array_other.put(key, new ArrayList<>(
							cadd_rec_for_ctx.
								stream().
								filter(R->R.otherKeyValues!=null).
								map(R->R.otherKeyValues.get(key)).
								filter(S->!(StringUtil.isBlank(S) || S.equals(".")|| S.equals("NA"))).
								collect(Collectors.toSet())
								));
						}
					
					}
					
				}
			if(!cadd_array_score.isEmpty() && got_non_null)
				{
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(this.CADD_FLAG_SCORE, cadd_array_score);
				vcb.attribute(this.CADD_FLAG_PHRED, cadd_array_phred);
				for(final String uf:this.userFields)
					{
					final List<String> cadd_other = cadd_array_other.get(uf);
					if(cadd_array_other==null || cadd_array_other.isEmpty()) continue;
					final List<String> array2 = new ArrayList<>(cadd_other.size());

					boolean all_na=true;
					boolean all_null=true;
					for(String v:cadd_other)
						{
						if(StringUtil.isBlank(v))
							{
							array2.add(".");
							}
						else
							{
							all_null = false;
							String s= v;
							if(StringUtil.isBlank(s)) s="NA";
							if(!s.equals("NA")) all_na=false;
							s=VCFUtils.escapeInfoField(s.replace(',','&'));
							array2.add(s);
							}
						}
					
					if(array2.isEmpty() || all_na || all_null) continue;
					vcb.attribute("CADD_"+uf, cadd_other);
					}
				
				return vcb.make();
				}
		return ctx;
		}
	
	private VCFHeaderLineCount getFieldHeaderLineCount(final String field) {
		return this.annoField2linecount.getOrDefault(field, VCFHeaderLineCount.UNBOUNDED);
	}
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator in, final VariantContextWriter out) {
			try {
			final VCFHeader header=in.getHeader();
			final OrderChecker<VariantContext> orderChecker = new OrderChecker<>(
					SequenceDictionaryUtils.extractRequired(header),
					false);
			
			
			JVarkitVersion.getInstance().addMetaData(this, header);
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
					" URI was " +this.ccaduri + 
					". We use "+NA_value+" for unknown value."
					));
			header.addMetaDataLine(new VCFInfoHeaderLine(
					this.CADD_FLAG_PHRED,
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc. " +
					" URI was " +this.ccaduri + 
					". We use "+NA_value+" for unknown value."
					));
			for(final String uf: this.userFields) {
				header.addMetaDataLine(new VCFInfoHeaderLine(
						"CADD_"+uf,
						this.getFieldHeaderLineCount(uf),
						VCFHeaderLineType.String,
						"User field extracted from " +this.ccaduri
						));
				}
			
			out.writeHeader(header);
			while(in.hasNext()) {
				out.add(runTabix(orderChecker.apply(in.next())));
				}
			return 0;
			} 
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected int beforeVcf()
		{
		try
			{
			if(this.CADD_FLAG_PHRED.equals(CADD_FLAG_SCORE)) {
				LOG.error("tag phred same as tag score");
				return -1;
				}
			if(StringUtils.isBlank(this.ccaduri) || !this.ccaduri.endsWith(".gz"))
				{
				LOG.error("CCAD uri should end with gz. got "+this.ccaduri);
				return -1;
				}
			
			this.tabix = new TabixFileReader(this.ccaduri);
			this.convertToCaddContigs = ContigNameConverter.fromContigSet(this.tabix.getChromosomes());
			
			for(;;)
				{
				final String head = this.tabix.readLine();
				if(!head.startsWith("#")) break;
				if(head.startsWith("#Chrom\tPos\tRef") ||
				   head.startsWith("#Chr\tPos\tRef"))
					{
					this.headerLabel = Arrays.asList(TAB.split(head));
					break;
					}
				}
			if(this.headerLabel==null || this.headerLabel.isEmpty())
				{
				LOG.error("Cannot read tabix file header");
				return -1;
				}
			this.headerLabel2column= new HashMap<>(this.headerLabel.size());
			for(int i=0;i< this.headerLabel.size();i++)
				{
				this.headerLabel2column.put(this.headerLabel.get(i), i);
				}
			Integer col = this.headerLabel2column.get("#Chrom");
			if(col==null) col = this.headerLabel2column.get("#Chr");
			if(col==null || col.intValue()!=0) {
				LOG.error("Illegal Header: Cannot get column index of #Chrom at 0");
				return -1;
				}
			col = this.headerLabel2column.get("Pos");
			if(col==null || col.intValue()!=1) {
				LOG.error("Illegal Header: Cannot get column index of 'Pos' at 1");
				return -1;
				}
			col = this.headerLabel2column.get("Ref");
			if(col==null || col.intValue()!=2) {
				LOG.error("Illegal Header: Cannot get column index of 'Ref' at 2");
				return -1;
				}
			col = this.headerLabel2column.get("Alt");
			if(col==null ) {
				LOG.error("Illegal Header: Cannot get column index of 'Alt'");
				return -1;
				}
			this.column_index_for_Alt = col.intValue();
			
			col = this.headerLabel2column.get("RawScore");
			if(col==null ) {
				LOG.error("Illegal Header: Cannot get column index of 'RawScore'");
				return -1;
				}
			this.column_index_for_RawScore = col.intValue();
			
			col = this.headerLabel2column.get("PHRED");
			if(col==null ) {
				LOG.error("Illegal Header: Cannot get column index of 'PHRED'");
				return -1;
				}
			this.column_index_for_PHRED = col.intValue();

			for(final String field:this.otherFieldsStr.split("[ \t,;]+"))
				{
				if(StringUtil.isBlank(field)) continue;
				col = this.headerLabel2column.get(field);
				if(col==null) {
					LOG.error("Illegal Header: Cannot find user's field "+field+ "in "+this.headerLabel);
					return -1;
					}
				this.userFields.add(field);
				}
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected void afterVcf() {
		if(this.tabix!=null) this.tabix.close();
		super.afterVcf();
		}

	public static void main(final String[] args)
		{
		new VcfCadd().instanceMainWithExit(args);
		}
	}
