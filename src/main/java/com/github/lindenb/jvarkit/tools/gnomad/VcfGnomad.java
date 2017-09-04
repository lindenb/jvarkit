package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.tools.vcfcmp.EqualRangeVcfIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/**
BEGIN_DOC
 
## Manifest
 
 the manifest is a tab delimited file containing 3 columns. It's used to map a contig to a URI
 
   * 1st column is a keyword 'exome' or 'genome'
   * 2d column is a contig name e.g: '1' .  Use '*' for 'any' chromosome
   * 3d column is a URL or file path where to find the data
 
 
## Example:
 
 ```
  curl -s "https://storage.googleapis.com/gnomad-public/release-170228/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz" |\
     gunzip -c | head -n 400 |\
     java  -jar ~/src/jvarkit-git/dist/vcfgnomad.jar -ac -gf IN_GNOMAD 

 (...)
 1	13595	.	AGT	A	379.68	AC0;IN_GNOMAD;RF	AB_HIST_ALL=0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_MEDIAN=1.44068e-01;AC=0;AC_AFR=0;AC_AMR=0;AC_ASJ=0;AC_EAS=0;AC_FIN=0;AC_Female=0;AC_Male=0;AC_NFE=0;AC_OTH=0;AC_POPMAX=.;AC_SAS=0;AC_raw=1;AF=0.00000e+00;AF_AFR=0.00000e+00;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=0.00000e+00;AF_FIN=0.00000e+00;AF_Female=0.00000e+00;AF_Male=0.00000e+00;AF_NFE=0.00000e+00;AF_OTH=0.00000e+00;AF_POPMAX=.;AF_SAS=0.00000e+00;AF_raw=9.99900e-06;AN=50778;AN_AFR=4986;AN_AMR=10892;AN_ASJ=1274;AN_EAS=7560;AN_FIN=694;AN_Female=24940;AN_Male=25838;AN_NFE=17556;AN_OTH=1486;AN_POPMAX=.;AN_SAS=6330;AN_raw=100010;AS_FilterStatus=RF|AC0;AS_RF=1.49748e-01;BaseQRankSum=-4.60000e-01;CSQ=-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034|YES|||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|6/6||ENST00000450305.2:n.561_562delTG||558-559||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.847_848delTG||844-845||||||1||1||deletion|1|HGNC|37102|YES||||||||||||3|||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene|||||||||||1|807|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|3/3||ENST00000515242.2:n.840_841delTG||837-838||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/4||ENST00000518655.2:n.678_679delTG||675-676||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene|||||||||||1|814|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site|||||||||||1||||deletion|1||||||||||||||||||||||||||||||||||||||||||||;ClippingRankSum=5.63000e-01;DP=2519792;DP_HIST_ALL=20921|3680|466|85|62|97|652|4365|4551|3656|2891|2039|1464|1114|954|811|688|497|352|310;DP_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;DP_MEDIAN=118;DREF_MEDIAN=3.98107e-38;FS=1.59250e+01;GC=25389,0,0;GC_AFR=2493,0,0;GC_AMR=5446,0,0;GC_ASJ=637,0,0;GC_EAS=3780,0,0;GC_FIN=347,0,0;GC_Female=12470,0,0;GC_Male=12919,0,0;GC_NFE=8778,0,0;GC_OTH=743,0,0;GC_SAS=3165,0,0;GC_raw=50004,1,0;GQ_HIST_ALL=11211|8535|2038|2055|803|203|195|95|28|49|65|37|115|64|88|117|164|34|237|23872;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;GQ_MEDIAN=99;Hom=0;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_Female=0;Hom_Male=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=0;Hom_raw=0;InbreedingCoeff=-4.37000e-02;MQ=3.15600e+01;MQRankSum=-8.97000e-01;POPMAX=.;QD=3.22000e+00;ReadPosRankSum=-1.23200e+00;SOR=1.09000e-01;VQSLOD=-1.83100e+00;VQSR_NEGATIVE_TRAIN_SITE;VQSR_culprit=QD;gnomad.exome.AC_AFR=0;gnomad.exome.AC_AMR=0;gnomad.exome.AC_ASJ=0;gnomad.exome.AC_EAS=0;gnomad.exome.AC_FIN=0;gnomad.exome.AC_Female=0;gnomad.exome.AC_Male=0;gnomad.exome.AC_NFE=0;gnomad.exome.AC_OTH=0;gnomad.exome.AC_raw=1;gnomad.exome.AN_AFR=4986;gnomad.exome.AN_AMR=10892;gnomad.exome.AN_ASJ=1274;gnomad.exome.AN_EAS=7560;gnomad.exome.AN_FIN=694;gnomad.exome.AN_Female=24940;gnomad.exome.AN_Male=25838;gnomad.exome.AN_NFE=17556;gnomad.exome.AN_OTH=1486;gnomad.exome.AN_raw=100010;gnomad.genome.AC_AFR=0;gnomad.genome.AC_AMR=0;gnomad.genome.AC_ASJ=0;gnomad.genome.AC_EAS=0;gnomad.genome.AC_FIN=0;gnomad.genome.AC_Female=0;gnomad.genome.AC_Male=0;gnomad.genome.AC_NFE=0;gnomad.genome.AC_OTH=0;gnomad.genome.AC_raw=1;gnomad.genome.AN_AFR=8680;gnomad.genome.AN_AMR=794;gnomad.genome.AN_ASJ=224;gnomad.genome.AN_EAS=1592;gnomad.genome.AN_FIN=3490;gnomad.genome.AN_Female=13274;gnomad.genome.AN_Male=16168;gnomad.genome.AN_NFE=13754;gnomad.genome.AN_OTH=908;gnomad.genome.AN_raw=30500

 
 ```

## Note to self: Another alternative with VariantAnnotator,

but I think it slower...

(javascript / Makefile generation)

```javascript
out.print(" ${java.exe} -jar ${gatk.jar} -R $(REF) -L $(addsuffix .tmp.vcf,$@) -T VariantAnnotator --variant $(addsuffix .tmp.vcf,$@) -o $(addsuffix .tmp2.vcf,$@) --resourceAlleleConcordance ");

out.print(" --resource:gnomad_exome /commun/data/pubdb/broadinstitute.org/gnomad/release-170228/vcf/exome/gnomad.exomes.r2.0.1.sites.vcf.gz ");
out.print("$(foreach A,${GFIELDS}, -E gnomad_exome.${A} ) ");

var genome="/commun/data/pubdb/broadinstitute.org/gnomad/release-170228/vcf/genome/gnomad.genomes.r2.0.1.sites."+chrom+".vcf.gz";

out.print("$(if $(realpath "+genome+"), --resource:gnomad_genome  "+genome+"  $(foreach A,${GFIELDS}, -E gnomad_genome.${A} ) )");
```



## Generating jar helper for knime
  
(for the people in my lab)
  
generate big jar

```
$ cd jvarkit
$ rm -rf tmp && mkdir tmp && echo '1.jar:2.jar:...N.jar:vcfgnomad.jar' | tr ":" "\n" | sort | uniq | while read F; do unzip -o $F -d tmp ; done && jar cvf vcfgnomad4knime.jar -C tmp . && rm -rf tmp
```
  
Open KNIME

we're going to create the following workflow : http://imgur.com/a/QcrKW

* create a new Node `java Snippet`
* in the tab 'additional libraries', add 'vcfgnomad4knime.jar'.
* in the tab 'java snippet'. Declare the following inputs: `c_CHROM,c_POS,c_REF,c_ALT`, the output string `GNOMAD`.

And insert the following code:

```java
// Your custom imports:
import  com.github.lindenb.jvarkit.tools.gnomad.VcfGnomad.KnimeAdapter;
// Enter your code here:

System.setProperty("http.proxyHost","cache.ha.univ-nantes.fr");
System.setProperty("https.proxyHost","cache.ha.univ-nantes.fr");
System.setProperty("http.proxyPort","3128");
System.setProperty("https.proxyPort","3128");

final KnimeAdapter app= new KnimeAdapter();
if(app.instanceMain(new String[]{"-ac",c_CHROM,String.valueOf(c_POS),c_REF,c_ALT})==0)
	{
	out_GNOMAD = app.getOutputString();
	}
else
	{
	out_GNOMAD = ".";
	}

// expression end

```
 

 END_DOC
 */
@Program(name="vcfgnomad",
	description="Peek annotations from gnomad",
	keywords={"vcf","annotation","gnomad"})
public class VcfGnomad extends Launcher{
	
	private static final Logger LOG = Logger.build(VcfGnomad.class).make();
	/** allele specific population in gnomad */
	private final static String POPS[]=new String[]{"AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH", "Male", "Female","SAS", "raw", "POPMAX"}; 
	/** 'ome'-type section */
	private enum OmeType {exome,genome};
	/** entries mapping chromosome/type->vcf.gz */
	private List<ManifestEntry> manifestEntries=new ArrayList<>();
	
	private final Map<String,Integer> contig2tid;
	
	
	private class ManifestEntry
		implements Closeable
		{
		OmeType omeType;
		String contig;
		String uri;
		
		/** when using Tabix reader */
		TabixVcfFileReader gnomad_tabix=null;
		/** when using vcf streaming */
		VcfIterator gnomad_vcf_iterator = null;
		EqualRangeVcfIterator gnomad_equal_range=null;
		
		int buffferChromEnd=0;
		final Map<ContigPosRef,VariantContext> buffer=new HashMap<>();
		@Override
		public void close() throws IOException {
			CloserUtil.close(gnomad_tabix);
			CloserUtil.close(gnomad_equal_range);
			CloserUtil.close(gnomad_vcf_iterator);
			this.buffer.clear();
			this.buffferChromEnd=0;
			this.gnomad_tabix=null;
			this.gnomad_vcf_iterator=null;
			}
		public void open()  throws IOException
			{
			if(VcfGnomad.this.streaming)
				{
				this.gnomad_vcf_iterator = VCFUtils.createVcfIterator(this.uri);
				final SAMSequenceDictionary dict = this.gnomad_vcf_iterator.getHeader().getSequenceDictionary();
				if(dict==null || dict.isEmpty()) throw new JvarkitException.VcfDictionaryMissing(this.uri);
				final Comparator<VariantContext> cmp = VCFUtils.createTidPosComparator(dict);
				this.gnomad_equal_range = new EqualRangeVcfIterator(this.gnomad_vcf_iterator, cmp);
				}
			else
				{
				this.gnomad_tabix=new TabixVcfFileReader(this.uri);
				}
			}
		
		
		
		/** find matching variant in tabix file, use a buffer to avoid multiple random accesses */
		VariantContext findMatching(final VariantContext userVariantCtx)
			{
			
			if( VcfGnomad.this.streaming) {
				try {
					final List<VariantContext> found = this.gnomad_equal_range.next(
							userVariantCtx
							);
					for(final VariantContext ctx:found)
						{
						if( !ctx.getReference().equals(userVariantCtx.getReference())) continue;
						if( VcfGnomad.this.filteredGnomad && ctx.isFiltered()) continue;
						if( VcfGnomad.this.noMultiAltGnomad && ctx.getAlternateAlleles().size()>1) continue;
						return ctx;
						}
					return null;
					}
				catch(final IOException err)
					{
					throw new RuntimeIOException(err);
					}
				}
			else
				{
				final ContigPosRef userCtx=new ContigPosRef(userVariantCtx);
				//past last buffer ? refill buffer
				if(this.buffferChromEnd <= userCtx.getPos())
					{
					buffer.clear();
					this.buffferChromEnd = userCtx.getPos() + VcfGnomad.this.gnomadBufferSize;
					final Iterator<VariantContext> iter=this.gnomad_tabix.iterator(
							userVariantCtx.getContig(),
							Math.max(0,userCtx.getPos()-1),
							this.buffferChromEnd
							);
					while(iter.hasNext())
						{
						final VariantContext ctx = iter.next();
						if( VcfGnomad.this.filteredGnomad && ctx.isFiltered()) continue;
						if( VcfGnomad.this.noMultiAltGnomad && ctx.getAlternateAlleles().size()>1) continue;
						final ContigPosRef key= new ContigPosRef(ctx);
						this.buffer.put(key,ctx);
						}
					CloserUtil.close(iter);
					}
				return this.buffer.get(userCtx);
				}
			}
		
		}
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-m","--manifest"},description="manifest file descibing how to map a contig to an URI . 3 columns: 1) exome|genome 2) contig 3) path or URL.")
	private File manifestFile=null;
	@Parameter(names={"-filtered","--filtered"},description="Skip Filtered User Variants")
	private boolean skipFiltered=false;
	@Parameter(names={"-filteredGnomad","--filteredGnomad"},description="[20170706] Skip Filtered GNOMAD Variants")
	private boolean filteredGnomad=false;
	@Parameter(names={"-noMultiAltGnomad","--noMultiAltGnomad"},description="[20170706] Skip Multi Allelic GNOMAD Variants")
	private boolean noMultiAltGnomad=false;
	@Parameter(names={"-gf","--gnomadFilter"},description="if defined, add this FILTER when the variant is found in nomad")
	private String inGnomadFilterName=null;
	@Parameter(names={"-ac","--alleleconcordance"},description="ALL Alt allele must be found in gnomad before setting a FILTER")
	private boolean alleleconcordance=false;
	@Parameter(names={"--noAlleleCount"},description="do Not Insert AC /Allele Count")
	private boolean doNotInsertAlleleCount=false;
	@Parameter(names={"--noAlleleNumber"},description="do Not Insert AN /Allele Number")
	private boolean doNotInsertAlleleNumber=false;
	@Parameter(names={"--noAlleleFreq"},description="do Not Insert AF /Allele Freq.")
	private boolean doNotInsertAlleleFreq=false;
	@Parameter(names={"--bufferSize"},description="When we're looking for variant in Exac, load the variants for 'N' bases instead of doing a random access for each variant")
	private int gnomadBufferSize=100000;
	@Parameter(names={"--streaming"},description="[20170707] Don't use tabix random-access (which are ok for small inputs) but you a streaming process (better to annotate a large WGS file). Assume dictionaries are sorted the same way.")
	private boolean streaming=false;
	
	
	private class InfoField
		{
		final OmeType ome;
		final String tag;
		final boolean is_AC;
		final VCFHeaderLineType lineType;
		final List<Object> attributes=new ArrayList<>();
		InfoField(String tag, OmeType ome,boolean is_AC,final VCFHeaderLineType lineType) {
			this.tag=tag;
			this.ome=ome;
			this.is_AC = is_AC;
			this.lineType=lineType;
			}
		public String getOutputTag() {
			return "gnomad_"+ this.ome.name()+"_"+this.tag;
		}
		
		VCFInfoHeaderLine makeVCFInfoHeaderLine()
			{
			if(!is_AC)
				{
				return new VCFInfoHeaderLine(
						getOutputTag(),1,
						this.lineType,
						"Field "+this.tag+" extracted from Gnomad ("+ome.name()+")"
						);
				}
			else
				{
				return new VCFInfoHeaderLine(
						getOutputTag(),VCFHeaderLineCount.A,
						this.lineType,
						"Field "+this.tag+" extracted from Gnomad ("+ome.name()+")"
						);
				}
			}
		void fill(final VariantContext ctx,final VariantContext gnomadCtx)
			{
			this.attributes.clear();
			
			if(!is_AC)
				{
				int att=gnomadCtx.getAttributeAsInt(this.tag, -9999);
				if(att>=0) {
					this.attributes.add(att);
					}
				else
					{
					this.attributes.add(null);
					}
				}
			else
				{
				final List<Allele> galts=gnomadCtx.getAlternateAlleles();
				final List<String> gatts = gnomadCtx.getAttributeAsStringList(this.tag,null);
				for(final Allele a:ctx.getAlternateAlleles())
					{
					Object found=null;
					//final int idx=gnomadCtx.getAlleleIndex(a);//non idx(REF)==0
					final int idx=galts.indexOf(a);
					
					if(idx>=0) {
						if(idx<gatts.size() && gatts.get(idx)!=null && !gatts.get(idx).equals(".")) {
							switch(this.lineType)
								{
								case Integer: found=Integer.parseInt(gatts.get(idx));break;
								case Float: found=Float.parseFloat(gatts.get(idx));break;
								default: throw new JvarkitException.ShouldNeverHappen(this.lineType.name());
								}
							}
						}
					this.attributes.add(found);
					}
				}
			}
		}
	
	private String normalizeContig(final String contig)
		{
		if(contig.startsWith("chr")) return contig.substring(3);
		return contig;
		}
	
	/** change user's chromosome for a variant if needed */
	private VariantContext normalizeVariantContig(final VariantContext userVariantCtx)
		{
		final String ensemblContig = normalizeContig(userVariantCtx.getContig());
		return ensemblContig.equals(userVariantCtx.getContig())
				? userVariantCtx
				: new VariantContextBuilder(userVariantCtx).chr(ensemblContig).make()
				;
		}
	
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VcfIterator iter,
			final VariantContextWriter out
			) {
		final ManifestEntry ome2manifest[]=new ManifestEntry[OmeType.values().length];
		Arrays.fill(ome2manifest,null);

		try {
			final List<InfoField> infoFields=new ArrayList<>();
			for(OmeType ome:OmeType.values()) {
				for(final String pop: POPS)
					{
					if(!doNotInsertAlleleCount) infoFields.add(new InfoField("AC_"+pop,ome,true,VCFHeaderLineType.Integer));
					if(!doNotInsertAlleleFreq) infoFields.add(new InfoField("AF_"+pop,ome,true,VCFHeaderLineType.Float));
					if(!doNotInsertAlleleNumber) infoFields.add(new InfoField("AN_"+pop,ome,pop.equals("POPMAX"),VCFHeaderLineType.Integer));
					}
				if(!doNotInsertAlleleCount) infoFields.add(new InfoField("AC",ome,true,VCFHeaderLineType.Integer));
				if(!doNotInsertAlleleFreq) infoFields.add(new InfoField("AF",ome,true,VCFHeaderLineType.Float));
				if(!doNotInsertAlleleNumber) infoFields.add(new InfoField("AN",ome,false,VCFHeaderLineType.Integer));
				}
			String prevContig=null;
			final VCFHeader h2=new VCFHeader(iter.getHeader());
			if(inGnomadFilterName!=null)
				{
				h2.addMetaDataLine(new VCFFilterHeaderLine(inGnomadFilterName,"Variant is in Gnomad"));
				}
			
			for(final InfoField infoField: infoFields)
				{
				h2.addMetaDataLine(infoField.makeVCFInfoHeaderLine());
				}
			out.writeHeader(h2);
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(h2).logger(LOG);
			while(iter.hasNext()) {
				final VariantContext ctx = progress.watch(iter.next());
				
				if(this.skipFiltered && ctx.isFiltered() )
					{
					out.add(ctx);
					continue;
					}
				final String ensemblContig=normalizeContig(ctx.getContig());

				/* CONTIG has changed, update the CONTIG */
				if(prevContig==null || !prevContig.equals(ctx.getContig())) {
					LOG.debug("Data for "+ctx.getContig());
					prevContig=ctx.getContig();
					for(OmeType ome: OmeType.values())
						{
						ManifestEntry newEntry = null;

						for(final ManifestEntry e: this.manifestEntries)
							{
							if(e.omeType!=ome) continue;
							if(e.contig.equals("*"))
								{
								//accept
								}
							else if(!e.contig.equals(ensemblContig)){
								continue;
								}
							newEntry=e;
							break;
							}
						if(newEntry==null)
							{
							LOG.warn("No Gnomad Data for "+ctx.getContig()+" / "+ome);
							}
						final ManifestEntry prevEntry = ome2manifest[ome.ordinal()];
						if(prevEntry==null && newEntry==null){
							ome2manifest[ome.ordinal()]=null;
							}
						else if(newEntry!=null && prevEntry!=null &&
								prevEntry.uri.equals(newEntry.uri))
							{
							// no need to re-open
							//continue with prev entry
							}
						else if(newEntry==null && prevEntry!=null)
							{
							prevEntry.close();
							ome2manifest[ome.ordinal()]=null;
							}
						else
							{
							if(prevEntry!=null) prevEntry.close();
							ome2manifest[ome.ordinal()]=newEntry;
							LOG.info("opening "+newEntry.uri);
							newEntry.open();
							}
						}
					}
				/** END UPDATE CONTIG */
				
				for(final InfoField infoField: infoFields)
					{
					infoField.attributes.clear();
					}
				
				boolean setfilter=false;
				
				// lopp over exome and genome data
				for(int i=0;i< ome2manifest.length;++i) {
					ManifestEntry entry = ome2manifest[i];
					if(entry==null) continue;
					
					final VariantContext ctx2=entry.findMatching(normalizeVariantContig(ctx));
					if(ctx2==null) continue;
					for(final InfoField infoField: infoFields)
						{
						if(infoField.ome!=entry.omeType) continue;
						infoField.fill(ctx, ctx2);
						}
					if(this.alleleconcordance)
						{
						//stream all ALT. return false if we found one ALT that is not found in Gnomad
						setfilter = !ctx.getAlternateAlleles().stream().
								filter(A->!ctx2.getAlternateAlleles().contains(A)).
								findAny().isPresent();
						}
					else
						{
						setfilter=true;
						}					
					}
				
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				if(setfilter && this.inGnomadFilterName!=null)
					{
					vcb.filter(inGnomadFilterName);
					}
				for(final InfoField infoField: infoFields)
					{
					if(infoField.attributes.isEmpty()) continue;
					if(!infoField.attributes.stream().filter(N->N!=null).findAny().isPresent()) continue;
					vcb.attribute(infoField.getOutputTag(), infoField.attributes);
					}
				out.add(vcb.make());
				}
			progress.finish();
			return 0;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		}
		finally {
			CloserUtil.close(Arrays.asList(ome2manifest));
		}
	}

	public VcfGnomad() {
		contig2tid=new HashMap<>(25);
		for(int i=1;i<=22;++i) contig2tid.put(String.valueOf(i), i);
		contig2tid.put("X",23);
		contig2tid.put("Y",24);
		
	}
	
@Override
public int doWork(final List<String> args) {
	if(this.gnomadBufferSize < 10) {
		this.gnomadBufferSize  = 10;
		}
	if(this.manifestFile==null)
		{
		LOG.info("Building default manifest file...");
		for(OmeType ot: OmeType.values()) {
			for(int i=1;i<= 23;++i) {
				final ManifestEntry entry=new ManifestEntry();
				entry.omeType=ot;
				switch(i)
					{
					default: entry.contig=String.valueOf(i);break;
					case 23: entry.contig="X";break;					
					}
				if(ot==OmeType.genome)
					{
					entry.uri = "https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites."+entry.contig+".vcf.gz";
					}
				else
					{
					entry.uri = "https://storage.googleapis.com/gnomad-public/release-170228/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz";
					}
				this.manifestEntries.add(entry);
				}
			}
		LOG.info("Building default manifest file... Done");
		}
	else
		{
		try {
		Files.lines(this.manifestFile.toPath()).forEach(L->{
			if(L.startsWith("#") || L.trim().isEmpty()) return;
			final String tokens[]=L.split("[\t]");
			if(tokens.length<3) throw new JvarkitException.TokenErrors("Expected 3 words",tokens);
			final ManifestEntry entry=new ManifestEntry();
			entry.omeType=OmeType.valueOf(tokens[0]);
			entry.contig = tokens[1].trim();
			entry.uri=tokens[2].trim();
			VcfGnomad.this.manifestEntries.add(entry);
		});
		} catch(final IOException err) {
			LOG.error(err);
			return -1;
		}
		}
	return doVcfToVcf(args, outputFile);
	}

/** For SolenaLS & Julien in knime. temporary solution */
@Deprecated
public static final class KnimeAdapter extends VcfGnomad
	{
	private String outputString="";
	private static class SingletonVcfIterator
		extends AbstractIterator<VariantContext>
		implements VcfIterator
		{
		final VCFHeader header=new VCFHeader();
		final Iterator<VariantContext> delegate;
		SingletonVcfIterator(VariantContext ctx)
			{
			this.delegate = Collections.singletonList(ctx).iterator();
			}
		@Override
		public AbstractVCFCodec getCodec() {
			throw new JvarkitException.ShouldNeverHappen("getCodec");
			}
		@Override
		public VCFHeader getHeader() {
			return header;
			}
		@Override
		protected VariantContext advance() {
			return this.delegate.hasNext()?this.delegate.next():null;
			}
		@Override
		public void close() throws IOException {
			
			}
		}
	private static class SingletonWriter
		implements VariantContextWriter
		{
		VariantContext variant=null;
		@Override
		public void writeHeader(VCFHeader arg0) { }
		@Override
		public void add(VariantContext ctx) {
			this.variant=ctx;
			}
		@Override
		public void close() { }
		@Override
		public boolean checkError() {return false;}
		}
	@Override
	protected int doVcfToVcf(final List<String> inputs,final File outorNull) {
		super.gnomadBufferSize=5;
		if(inputs.size()!=4) throw new JvarkitException.UserError("expected 4 fields : CHROM/POS/REF/ALTS");
		final List<Allele> alleles=new ArrayList<>();
		final Allele ref= Allele.create(inputs.get(2),true);
		final long start = Long.parseLong(inputs.get(1));
		alleles.add(ref);
		alleles.addAll(Arrays.stream(inputs.get(3).split("[,]")).map(
				S->Allele.create(S,false)
				).collect(Collectors.toList()));
		final VariantContext  vc=new VariantContextBuilder("knime",
				inputs.get(0),
				start,
				start + ref.length() -1,
				alleles
				).make();
		final SingletonWriter w=new SingletonWriter();
		int ret= doVcfToVcf("knime", new SingletonVcfIterator(vc), w);
		if(ret!=0 || w.variant==null) return -1;
		for(final String key:w.variant.getAttributes().keySet())
			{
			if(!(this.outputString==null || this.outputString.isEmpty()))
				{
				this.outputString+=";";
				}
			this.outputString+=key+"="+String.join(",",w.variant.getAttributeAsStringList(key, "."));
			}
		return ret;
		}
	public String getOutputString()
		{
		return this.outputString;
		}
	public static void main(String[] args) {
		final KnimeAdapter app= new KnimeAdapter();
		if(app.instanceMain(args)==0)
			{
			app.stdout().println(app.outputString);
			};
		}

	}
	
public static void main(String[] args) {
	new VcfGnomad().instanceMainWithExit(args);
	}
}
