/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
BEGIN_DOC
 
## Manifest
 
 the manifest is a tab delimited file containing 3 columns. It's used to map a contig to a URI
 
   * 1st column is a keyword 'exome' or 'genome'
   * 2d column is a contig name e.g: '1' .  Use '*' for 'any' chromosome
   * 3d column is a URL or file path where to find the data
 

```
$ cat ./gnomad.manifest 
exome   *       /commun/data/pubdb/broadinstitute.org/gnomad/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.gz
genome  *       /commun/data/pubdb/broadinstitute.org/gnomad/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.gz
```

## Example:
 
 ```
  curl -s "https://storage.googleapis.com/gnomad-public/release-170228/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz" |\
     gunzip -c | head -n 400 |\
     java  -jar ~/src/jvarkit-git/dist/vcfgnomad.jar -ac -gf IN_GNOMAD 

 (...)
 1	13595	.	AGT	A	379.68	AC0;IN_GNOMAD;RF	AB_HIST_ALL=0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_MEDIAN=1.44068e-01;AC=0;AC_AFR=0;AC_AMR=0;AC_ASJ=0;AC_EAS=0;AC_FIN=0;AC_Female=0;AC_Male=0;AC_NFE=0;AC_OTH=0;AC_POPMAX=.;AC_SAS=0;AC_raw=1;AF=0.00000e+00;AF_AFR=0.00000e+00;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=0.00000e+00;AF_FIN=0.00000e+00;AF_Female=0.00000e+00;AF_Male=0.00000e+00;AF_NFE=0.00000e+00;AF_OTH=0.00000e+00;AF_POPMAX=.;AF_SAS=0.00000e+00;AF_raw=9.99900e-06;AN=50778;AN_AFR=4986;AN_AMR=10892;AN_ASJ=1274;AN_EAS=7560;AN_FIN=694;AN_Female=24940;AN_Male=25838;AN_NFE=17556;AN_OTH=1486;AN_POPMAX=.;AN_SAS=6330;AN_raw=100010;AS_FilterStatus=RF|AC0;AS_RF=1.49748e-01;BaseQRankSum=-4.60000e-01;CSQ=-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034|YES|||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|6/6||ENST00000450305.2:n.561_562delTG||558-559||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.847_848delTG||844-845||||||1||1||deletion|1|HGNC|37102|YES||||||||||||3|||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene|||||||||||1|807|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|3/3||ENST00000515242.2:n.840_841delTG||837-838||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/4||ENST00000518655.2:n.678_679delTG||675-676||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene|||||||||||1|814|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site|||||||||||1||||deletion|1||||||||||||||||||||||||||||||||||||||||||||;ClippingRankSum=5.63000e-01;DP=2519792;DP_HIST_ALL=20921|3680|466|85|62|97|652|4365|4551|3656|2891|2039|1464|1114|954|811|688|497|352|310;DP_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;DP_MEDIAN=118;DREF_MEDIAN=3.98107e-38;FS=1.59250e+01;GC=25389,0,0;GC_AFR=2493,0,0;GC_AMR=5446,0,0;GC_ASJ=637,0,0;GC_EAS=3780,0,0;GC_FIN=347,0,0;GC_Female=12470,0,0;GC_Male=12919,0,0;GC_NFE=8778,0,0;GC_OTH=743,0,0;GC_SAS=3165,0,0;GC_raw=50004,1,0;GQ_HIST_ALL=11211|8535|2038|2055|803|203|195|95|28|49|65|37|115|64|88|117|164|34|237|23872;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;GQ_MEDIAN=99;Hom=0;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_Female=0;Hom_Male=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=0;Hom_raw=0;InbreedingCoeff=-4.37000e-02;MQ=3.15600e+01;MQRankSum=-8.97000e-01;POPMAX=.;QD=3.22000e+00;ReadPosRankSum=-1.23200e+00;SOR=1.09000e-01;VQSLOD=-1.83100e+00;VQSR_NEGATIVE_TRAIN_SITE;VQSR_culprit=QD;gnomad.exome.AC_AFR=0;gnomad.exome.AC_AMR=0;gnomad.exome.AC_ASJ=0;gnomad.exome.AC_EAS=0;gnomad.exome.AC_FIN=0;gnomad.exome.AC_Female=0;gnomad.exome.AC_Male=0;gnomad.exome.AC_NFE=0;gnomad.exome.AC_OTH=0;gnomad.exome.AC_raw=1;gnomad.exome.AN_AFR=4986;gnomad.exome.AN_AMR=10892;gnomad.exome.AN_ASJ=1274;gnomad.exome.AN_EAS=7560;gnomad.exome.AN_FIN=694;gnomad.exome.AN_Female=24940;gnomad.exome.AN_Male=25838;gnomad.exome.AN_NFE=17556;gnomad.exome.AN_OTH=1486;gnomad.exome.AN_raw=100010;gnomad.genome.AC_AFR=0;gnomad.genome.AC_AMR=0;gnomad.genome.AC_ASJ=0;gnomad.genome.AC_EAS=0;gnomad.genome.AC_FIN=0;gnomad.genome.AC_Female=0;gnomad.genome.AC_Male=0;gnomad.genome.AC_NFE=0;gnomad.genome.AC_OTH=0;gnomad.genome.AC_raw=1;gnomad.genome.AN_AFR=8680;gnomad.genome.AN_AMR=794;gnomad.genome.AN_ASJ=224;gnomad.genome.AN_EAS=1592;gnomad.genome.AN_FIN=3490;gnomad.genome.AN_Female=13274;gnomad.genome.AN_Male=16168;gnomad.genome.AN_NFE=13754;gnomad.genome.AN_OTH=908;gnomad.genome.AN_raw=30500

 
 ```

## Note to self: Another alternative with VariantAnnotator,
manifestFile
but I think it slower...

(javascript / Makefile generation)

```javascript
out.print(" ${java.exe} -jar ${gatk.jar} -R $(REF) -L $(addsuffix .tmp.vcf,$@) -T VariantAnnotator --variant $(addsuffix .tmp.vcf,$@) -o $(addsuffix .tmp2.vcf,$@) --resourceAlleleConcordance ");

out.print(" --resource:gnomad_exome /commun/data/pubdb/broadinstitute.org/gnomad/release-170228/vcf/exome/gnomad.exomes.r2.0.1.sites.vcf.gz ");
out.print("$(foreach A,${GFIELDS}, -E gnomad_exome.${A} ) ");

var genome="/commun/data/pubdb/broadinstitute.org/gnomad/release-170228/vcf/genome/gnomad.genomes.r2.0.1.sites."+chrom+".vcf.gz";

out.print("$(if $(realpath "+genome+"), --resource:gnomad_genome  "+genome+"  $(foreach A,${GFIELDS}, -E gnomad_genome.${A} ) )");
```

## History

  * 20181214 : keep gnomad FILTERs
  * 20181127 : rewritten for gnomad 2.1

END_DOC
 */
@Program(name="vcfgnomad",
	description="Peek annotations from gnomad",
	keywords={"vcf","annotation","gnomad"},
	modificationDate="20190308"
)
public class VcfGnomad extends Launcher{
	
	private static final Logger LOG = Logger.build(VcfGnomad.class).make();
	/** 'ome'-type section */
	/* private */ enum OmeType {exome,genome};
	
	private enum GnomadVersion {
		v2_0,
		v2_1
		}
	private final GnomadVersion DEFAULT_GNOMAD_VERSION = GnomadVersion.v2_0;
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-m","--manifest"},description="manifest file descibing how to map a contig to an URI . 3 columns: 1) exome|genome 2) contig 3) path or URL.")
	private File manifestFile=null;
	@Parameter(names={"--bufferSize"},description="When we're looking for variant in Gnomad, load the variants for 'N' bases instead of doing a random access for each variant. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int gnomadBufferSize= 10_000;
	@Parameter(names={"-filtered","--filtered"},description="Skip Filtered User Variants")
	private boolean skipFiltered=false;
	@Parameter(names={"-filteredGnomad","--filteredGnomad"},description="[20170706] Skip Filtered GNOMAD Variants")
	private boolean filteredGnomad=false;
	@Parameter(names={"-noMultiAltGnomad","--noMultiAltGnomad"},description="[20170706] Skip Multi Allelic GNOMAD Variants")
	private boolean noMultiAltGnomad=false;
	@Parameter(names={"--noUpdateId"},description="do Not Update ID if it is missing in user's variant")
	private boolean doNotUpdateId=false;
	@Parameter(names={"-gf","--gnomadFilter"},description="if defined, add this FILTER when any variant [CHROM:POS:REF] is found in nomad")
	private String inGnomadFilterName=null;
	@Parameter(names={"-of","--overlapFilter"},description="if defined, add this FILTER when any variant overlapping [CHROM:POS] is found in nomad")
	private String overlapGnomadFilterName=null;
	@Parameter(names={"-gnomad-filter-prefix","--gnomad-filter-prefix"},description="[20181214] if not empty, include the Gnomad FILTERs using this prefix.")
	private String filteredInGnomadFilterPrefix="GNOMAD";
	@Parameter(names={"--genome"},description="[20180327] For @MKarakachoff : genome only, don't use 'exome data'")
	private boolean useGenomeOnly = false;
	@Parameter(names={"--exclude"},description="[20180327] exclude gnomad INFO field matching this regular expression. Empty: accept all")
	private String excludePatternStr = "controls|non_cancer|non_neuro|non_topmed";
	@Parameter(names={"--ani"},description="[20190311] for allele numbers 'AN' to be variant-count-type=Integer (not 'A' as declared in gnomad)")
	private boolean alleleNumber_is_integer = false;
	@Parameter(names={"--ignore-error0"},description="[20190429] ignore error when gnomad/INFO is found twice for the same position. I found the error after a liftover to hg38. see https://twitter.com/yokofakun/status/1122814203381858305")
	private boolean ignore_info_found_twice = false;
	
	/** entries mapping chromosome/type->vcf.gz */
	private List<ManifestEntry> manifestEntries=new ArrayList<>();
	private GnomadVersion gnomadVersion = DEFAULT_GNOMAD_VERSION;


	
	private class ManifestEntry
	implements Closeable
		{
		OmeType omeType;
		String contig;
		String uri;
		/** when using Tabix reader */
		private TabixVcfFileReader gnomad_tabix=null;
		
		private Interval lastInterval = null;
		final List<VariantContext> buffer = new ArrayList<>();
		/** convert to gnomad notation. Since I lift overred the VCF to hg38 */
		private ContigNameConverter ctgNameConverter = ContigNameConverter.getIdentity();
		
		@Override
		public void close() {
			CloserUtil.close(this.gnomad_tabix);
			this.buffer.clear();
			this.lastInterval =  null;
			this.gnomad_tabix = null;
			}
		
		public void open()
			{
			try {
				this.gnomad_tabix=new TabixVcfFileReader(this.uri);
				this.ctgNameConverter = ContigNameConverter.fromContigSet(this.gnomad_tabix.getChromosomes());
				}
			catch(final IOException err)
				{
				throw new RuntimeIOException("Cannot open "+this.uri,err);
				}
			}
		
		public VCFHeader getHeader() {
			return Objects.requireNonNull(this.gnomad_tabix).getHeader();
			}	
		
		private String normalizeContig(final String s) {
			return this.ctgNameConverter.apply(s);
			}
		
		boolean acceptContig(final String userCtg) {
			final String normContig = this.normalizeContig(userCtg);
			if(StringUtil.isBlank(normContig)) return false;
			if(this.contig.equals("*")) return true;
			return this.contig.equals(userCtg);
			}
		
		/** find matching variant in tabix file, use a buffer to avoid multiple random accesses */
		final List<VariantContext> findOverlapping(final VariantContext userVariantCtx)
			{
			if(!acceptContig(userVariantCtx.getContig())) return Collections.emptyList();
			final String normContig = this.normalizeContig(userVariantCtx.getContig());
			
			if(this.lastInterval==null ||
				!this.lastInterval.getContig().equals(normContig) ||
				!CoordMath.encloses(lastInterval.getStart(), lastInterval.getEnd(), userVariantCtx.getStart(), userVariantCtx.getEnd())
				)
				{
				this.buffer.clear();
				this.lastInterval = new Interval(
						normContig,
						Math.max(0, userVariantCtx.getStart()-10),
						userVariantCtx.getEnd()+ VcfGnomad.this.gnomadBufferSize
					);
				final Iterator<VariantContext> iter= this.gnomad_tabix.iterator(
						this.lastInterval.getContig(),
						this.lastInterval.getStart(),
						this.lastInterval.getEnd()
						);
				while(iter.hasNext())
					{
					final VariantContext ctx = iter.next();
					if( VcfGnomad.this.filteredGnomad && ctx.isFiltered()) continue;
					if( VcfGnomad.this.noMultiAltGnomad && ctx.getAlternateAlleles().size()>1) continue;
					this.buffer.add(ctx);
					}
				CloserUtil.close(iter);
				}
			
			return this.buffer.stream().
					filter(V->V.getContig().equals(normContig) &&  CoordMath.overlaps(V.getStart(), V.getEnd(), userVariantCtx.getStart(), userVariantCtx.getEnd())	).
					collect(Collectors.toList());
			}
		
		}
	
	
	private class InfoField
		{
		final OmeType ome;
		final VCFInfoHeaderLine original;
		private VCFInfoHeaderLine _outputheader = null;
		InfoField(final VCFInfoHeaderLine original, final OmeType ome) {
			this.original=original;
			this.ome=ome;
			}
		public String getOutputTag() {
			return "gnomad_"+ this.ome.name()+"_"+this.original.getID().toUpperCase();
			}
		
		public VCFInfoHeaderLine geOutputHeaderLine() {
			if(this._outputheader==null) {
					final String desc = "["+ome.name()+"]"+this.original.getDescription();
					if((this.original.isFixedCount() && this.original.getCount()==1 && gnomadVersion.equals(GnomadVersion.v2_0)) ||
						this.getOutputLineCount()==VCFHeaderLineCount.INTEGER && gnomadVersion.equals(GnomadVersion.v2_1))
						{
						this._outputheader = new VCFInfoHeaderLine(
								this.getOutputTag(),
								1,
								this.original.getType(),
								desc);
						}
					else
						{
						this._outputheader = new VCFInfoHeaderLine(
								this.getOutputTag(),
								this.getOutputLineCount(),
								this.original.getType(),
								desc);
						}
				}
			return this._outputheader;
			}
		/*
		 bcftools view /commun/data/pubdb/broadinstitute.org/gnomad/release-181127/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.gz |  grep -i AN_nfe,
		##INFO=<ID=AN_nfe,Number=A,Type=Integer,Description="Total number of alleles in samples of non-Finnish European ancestry">
		##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description="Alternate allele count for samples of non-Finnish European ancestry">

		 
		 
		 */
		public VCFHeaderLineCount getOutputLineCount() {
			if(this.original.isFixedCount() && 
				this.original.getCount()==1 && 
				gnomadVersion.equals(GnomadVersion.v2_0))
				{
				return VCFHeaderLineCount.INTEGER;
				}
			if(	VcfGnomad.this.alleleNumber_is_integer &&
					this.original.getCountType().equals(VCFHeaderLineCount.A) &&
					(this.original.getID().startsWith("AN_") || this.original.getID().contains("_AN_"))&&
					gnomadVersion.equals(GnomadVersion.v2_1))
					{
					return VCFHeaderLineCount.INTEGER;
					}
			return VCFHeaderLineCount.A;
			}
		
		Object getDefault() {
			switch(original.getType())
				{
				case Float: return 0.0f;
				case Integer : return 0;
				case String: return ".";
				default: throw new IllegalStateException(original.toString());
				}
			}
		
		Object parse(Object o) {
			if(o==null) return o;
			String str = String.valueOf(o);
			if(str.isEmpty() || str.equals(".")) return str;
			switch(original.getType())
				{
				case Float: return Float.parseFloat(str);
				case Integer : return Integer.parseInt(str);
				case String: return str;
				default: throw new IllegalStateException(original.toString());
				}
			}
		
		@Override
		public String toString() {
			return "\nINPUT : " +original.toString()+"\nOUTPUT: "+geOutputHeaderLine();
			}
		}
		
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter out
			)
		{
		final VCFHeader h0 = iter.getHeader();
		if(!SequenceDictionaryUtils.isGRCh37(h0)) {
			LOG.warn("Input is NOT GRCh37 ?");
			// can be lift over
			}
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(h0).logger(LOG).build();
		
		
		final VCFHeader h2 = new VCFHeader(h0);
		final Set<String> all_filters_from_gnomad = new HashSet<>();
		final List<InfoField> infoFields = new ArrayList<>();
		for(final OmeType ome: OmeType.values())
			{
			
			final ManifestEntry entry = this.manifestEntries.stream().
					filter(M->M.omeType.equals(ome)).
					findFirst().
					orElse(null);
			if(entry==null) continue;
			entry.open();
			final VCFHeader header= entry.getHeader();
			
			this.gnomadVersion = header.getInfoHeaderLines().stream().anyMatch(V->V.getID().equals("non_neuro_AC_nfe_male"))?
					GnomadVersion.v2_1:
					GnomadVersion.v2_0
					;
			LOG.debug("identified as gnomad version "+this.gnomadVersion);
			
			if(!StringUtil.isBlank(this.filteredInGnomadFilterPrefix))
				{
				for(final VCFFilterHeaderLine fh: header.getFilterLines())
					{
					if(fh.getID().equals(VCFConstants.PASSES_FILTERS_v4)) continue;
					final String fid = this.filteredInGnomadFilterPrefix+"_"+ome.name().toUpperCase()+"_"+ fh.getID();
					final VCFFilterHeaderLine fh2 = new VCFFilterHeaderLine(
							fid,
							"[gnomad-"+ome.name()+"]" + fh.getDescription()
							);
					h2.addMetaDataLine(fh2);
					all_filters_from_gnomad.add(fid);
					}
				}
			
			
			final Predicate<VCFInfoHeaderLine> acceptInfoTag;
			if(StringUtil.isBlank(this.excludePatternStr))
				{
				acceptInfoTag = T->true;
				}
			else
				{
				final Pattern pat = Pattern.compile(this.excludePatternStr);
				acceptInfoTag = T->!pat.matcher(T.getID()).find();
				}
			
			switch(this.gnomadVersion) {
				case v2_0:
					header.getInfoHeaderLines().
						stream().
						filter(acceptInfoTag).
						filter(FH->
							FH.getID().equals("AC") ||
							FH.getID().equals("AN")  ||
							FH.getID().equals("AF") ||
							FH.getID().startsWith("AC_") ||
							FH.getID().startsWith("AN_") ||
							FH.getID().startsWith("AF_")).
						map(FH->new InfoField(FH,ome)).
						forEach(FH->infoFields.add(FH));
					break;
				case v2_1:
					header.getInfoHeaderLines().stream().
						filter(acceptInfoTag).
						filter(FH->!FH.getID().contains("MEDIAN") ).
						filter(FH->
							FH.getID().contains("AC") ||
							FH.getID().contains("AN") ||
							FH.getID().contains("AF")).
						map(FH->new InfoField(FH,ome)).
						forEach(FH->infoFields.add(FH));
					break;
				default:
					throw new IllegalStateException("TODO "+this.gnomadVersion);
				}
			entry.close();
			}
		
	
		
		for(final InfoField f: infoFields) h2.addMetaDataLine(f.geOutputHeaderLine());
		
		if(!StringUtil.isBlank(this.inGnomadFilterName)) {
			h2.addMetaDataLine(new VCFFilterHeaderLine(this.inGnomadFilterName,"Variant CHROM/POS/REF was found in gnomad"));
			}
		
		if(!StringUtil.isBlank(this.overlapGnomadFilterName)) {
			h2.addMetaDataLine(new VCFFilterHeaderLine(this.overlapGnomadFilterName,"Gnomad Variant was found overlapping the variant"));
			}
		
		JVarkitVersion.getInstance().addMetaData(this, h2);
		out.writeHeader(h2);
		
		final ManifestEntry om2manifest[] = new ManifestEntry[]{null,null};
		while(iter.hasNext()) {
			final VariantContext ctx = progress.apply(iter.next());
			
			final Set<String> filters = new HashSet<>(ctx.getFilters());
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);

			for(final InfoField f: infoFields) vcb.rmAttribute(f.getOutputTag());
			
			if(!StringUtil.isBlank(this.inGnomadFilterName)) {
				filters.remove(this.inGnomadFilterName);
				}
			
			filters.removeAll(all_filters_from_gnomad);
			
			if(!StringUtil.isBlank(this.overlapGnomadFilterName)) {
				filters.remove(this.overlapGnomadFilterName);
				}
			
			if(this.skipFiltered && ctx.isFiltered() )
				{
				vcb.filters(filters);
				out.add(vcb.make());
				continue;
				}
			if(ctx.getContig().equals("MT") || ctx.getContig().equals("chrM") || ctx.getContig().contains("_")) {
				vcb.filters(filters);
				out.add(vcb.make());
				continue;
				}
			final List<Allele> alternateAlleles = ctx.getAlternateAlleles();
			String newid = null;
			boolean set_filter_ctx_is_in_gnomad = false;
			boolean found_gnomad_overlapping_variant = false;
			
			for(int omeIndex=0;omeIndex<2;omeIndex++)
				{
				final OmeType omeType = omeIndex==0?OmeType.exome:OmeType.genome;
				if(this.useGenomeOnly && !omeType.equals(OmeType.genome)) continue;
				if(om2manifest[omeIndex]!=null && !om2manifest[omeIndex].acceptContig(ctx.getContig()))
					{
					om2manifest[omeIndex].close();
					om2manifest[omeIndex]=null;
					}
				if(om2manifest[omeIndex]==null)
					{
					om2manifest[omeIndex] = this.manifestEntries.stream().filter(M->M.omeType.equals(omeType) && M.acceptContig(ctx.getContig())).findFirst().orElse(null);
					if(om2manifest[omeIndex]==null) continue;
					LOG.debug("Opening "+om2manifest[omeIndex].uri);
					om2manifest[omeIndex].open();
					}
				// variant overlapping 'ctx'
				final List<VariantContext> overlappingVariants = om2manifest[omeIndex].findOverlapping(ctx);
				if(!overlappingVariants.isEmpty()) found_gnomad_overlapping_variant = true;
				
				
				final List<VariantContext> gnomadVariants = overlappingVariants.
							stream().
							filter(V->V.getStart()==ctx.getStart() && V.getReference().equals(ctx.getReference())).
							collect(Collectors.toList());

				
				
				if(!gnomadVariants.isEmpty()) {
					set_filter_ctx_is_in_gnomad=true;
					
					// set new id ?
					if( newid == null) {
						newid = gnomadVariants.
								stream().
								filter(V->V.hasID()).
								map(V->V.getID()).
								findFirst().
								orElse(null);
						}
					
					// add FILTER(s)
					if(!StringUtil.isBlank(this.filteredInGnomadFilterPrefix)) {
						filters.addAll(
							gnomadVariants.
							stream().
							filter(V->V.isFiltered()).
							flatMap(V->V.getFilters().stream()).
							filter(F->!F.equals(VCFConstants.PASSES_FILTERS_v4)).
							map(F->this.filteredInGnomadFilterPrefix+"_"+omeType.name().toUpperCase()+"_"+F).
							collect(Collectors.toList())
							);
						}
					}
				
				
				// loop over each field
				for(final InfoField infoField: infoFields)
					{
					if(!infoField.ome.equals(omeType)) continue;
					if(gnomadVariants.stream().noneMatch(V->V.hasAttribute(infoField.original.getID()))) continue;
					
					
					if(infoField.getOutputLineCount().equals(VCFHeaderLineCount.INTEGER) && 
						infoField.original.getCountType().equals(VCFHeaderLineCount.INTEGER)) {
						final Set<Object> set = gnomadVariants.
									stream().
									filter(V->V.hasAttribute(infoField.original.getID())).
									map(V->infoField.parse(V.getAttribute(infoField.original.getID()))).
									collect(Collectors.toSet());
						if(set.isEmpty()) continue;
						if(set.size()==1) {
							vcb.attribute(infoField.getOutputTag(),set.iterator().next());
							}
						else if(this.ignore_info_found_twice)
							{
							LOG.warn("Found more than one value ("+set+") for "+infoField+" "+ ctx.getContig()+":"+ctx.getStart());
							vcb.attribute(infoField.getOutputTag(),set.iterator().next());
							}
						else
							{
							LOG.error("Found more than one value ("+set+") for "+infoField+" "+ ctx.getContig()+":"+ctx.getStart()+". Are you using a lift-overed gnomad ? Use option --ignore-error0 to skip this error");
							progress.close();
							return -1;
							}
						}
					else if(infoField.original.getCountType()== VCFHeaderLineCount.A)
						{
						final Object numbers[]=new Object[alternateAlleles.size()];
						Arrays.fill(numbers, infoField.getDefault());
						for(int x=0;x< alternateAlleles.size();++x)
							{
							final Allele alt=alternateAlleles.get(x);
							if(alt.equals(Allele.SPAN_DEL)) continue;
							for(final VariantContext gv:gnomadVariants)
								{
								int idx = gv.getAlternateAlleles().indexOf(alt);
								if(idx==-1) continue;
								if(!gv.hasAttribute(infoField.original.getID())) continue;
								final List<Object> array = gv.getAttributeAsList(infoField.original.getID());
								if(idx>=array.size()) continue;
								numbers[x] = infoField.parse(array.get(idx));
								}
							}
						vcb.attribute(infoField.getOutputTag(), numbers);
						}
					else if(infoField.original.isFixedCount() && infoField.original.getCount()==1) {
						final Object numbers[]=new Object[alternateAlleles.size()];
						Arrays.fill(numbers, infoField.getDefault());
						boolean something_found_for_an_allele = false;
						for(int x=0;x< alternateAlleles.size();++x)
							{
							final Allele alt=alternateAlleles.get(x);
							if(alt.equals(Allele.SPAN_DEL)) continue;
							for(final VariantContext gv:gnomadVariants)
								{
								int idx = gv.getAlternateAlleles().indexOf(alt);
								if(idx==-1) continue;
								if(!gv.hasAttribute(infoField.original.getID())) continue;
								something_found_for_an_allele = true;
								final Object val=  infoField.parse(gv.getAttribute(infoField.original.getID()));
								if(val==null) continue;
								numbers[x] = val;
								
								}
							}
						if(something_found_for_an_allele) {
							vcb.attribute(infoField.getOutputTag(), numbers);
							}
						}
					else
						{
						LOG.error("not handled "+infoField);
						progress.close();
						return -1;
						}
					}
				}
			if(set_filter_ctx_is_in_gnomad && !StringUtil.isBlank(this.inGnomadFilterName)) {
				filters.add(this.inGnomadFilterName);
				}
			
			if(found_gnomad_overlapping_variant && !StringUtil.isBlank(this.overlapGnomadFilterName)) {
				filters.add(this.overlapGnomadFilterName);
			}
			
			if(!this.doNotUpdateId && !ctx.hasID() && !StringUtil.isBlank(newid)) vcb.id(newid);
			vcb.filters(filters);
			out.add(vcb.make());
			}
		out.close();
		progress.close();
		for(int omeIndex=0;omeIndex<2;omeIndex++) CloserUtil.close(om2manifest[omeIndex]);
		return 0;
		}

	public VcfGnomad() {
		}
	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			if(this.gnomadBufferSize < 10) {
				LOG.error("buffer size is too small "+this.gnomadBufferSize);
				return -1;
				}
			if(this.manifestFile==null)
				{
				LOG.info("Building default manifest file...");
				for(final OmeType ot: OmeType.values()) {
					if(this.useGenomeOnly && !ot.equals(OmeType.genome)) continue;
					for(int i=1;i<= 24;++i) {
						final ManifestEntry entry=new ManifestEntry();
						entry.omeType=ot;
						switch(i)
							{
							case 23: entry.contig="X";break;	
							case 24: entry.contig="Y";break;
							default: entry.contig=String.valueOf(i);break;
							}
						switch(gnomadVersion)
							{
							case v2_1:
								if(ot==OmeType.genome)
									{
									if(i==24) continue;//no "chrY" for this version for genome
									entry.uri = "https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr"+entry.contig+".vcf.bgz";
									}
								else
									{
									entry.uri = "https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr"+entry.contig+".vcf.bgz";
									}
								break;
							default: {
								entry.close();
								LOG.error("Building default manifest is not available for version: "+this.gnomadVersion);
								return -1;
								}
							}
						this.manifestEntries.add(entry);
						}
					}
				LOG.info("Building default manifest file... Done");
				}
			else
				{
				try {
					final CharSplitter tab = CharSplitter.TAB;
					Files.lines(this.manifestFile.toPath()).forEach(L->{
						if(L.startsWith("#") || StringUtil.isBlank(L)) return;
						final String tokens[]=tab.split(L);
						if(tokens.length<3) throw new JvarkitException.TokenErrors("Expected 3 words",tokens);
						final ManifestEntry entry=new ManifestEntry();
						entry.omeType=OmeType.valueOf(tokens[0]);
						if(this.useGenomeOnly && !entry.omeType.equals(OmeType.genome)) {
							entry.close();
							return;
							}
						entry.contig = tokens[1].trim();
						entry.uri=tokens[2].trim();
						this.manifestEntries.add(entry);
						});
				} catch(final IOException err) {
					LOG.error(err);
					return -1;
					}
				}			
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}

public static void main(final String[] args) {
	new VcfGnomad().instanceMainWithExit(args);
	}
}
