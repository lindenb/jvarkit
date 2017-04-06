package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/**
 BEGIN_DOC
 
 ## Manifest
 
 the manifest is a tab delimited file containing 3 columns. It's used to map a contig to a URI
 
 1st column is a keyword 'exome' or 'genome'
 2d column is a contig name e.g: '1' .  Use '*' for 'any' chromosome
 3d column is a URL or file path where to find the data
 
 
 ## Example:
 
 ```
  curl -s "https://storage.googleapis.com/gnomad-public/release-170228/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz" |\
     gunzip -c | head -n 400 |\
     java  -jar ~/src/jvarkit-git/dist/vcfgnomad.jar -ac -gf IN_GNOMAD 

 (...)
 1	13595	.	AGT	A	379.68	AC0;IN_GNOMAD;RF	AB_HIST_ALL=0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_MEDIAN=1.44068e-01;AC=0;AC_AFR=0;AC_AMR=0;AC_ASJ=0;AC_EAS=0;AC_FIN=0;AC_Female=0;AC_Male=0;AC_NFE=0;AC_OTH=0;AC_POPMAX=.;AC_SAS=0;AC_raw=1;AF=0.00000e+00;AF_AFR=0.00000e+00;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=0.00000e+00;AF_FIN=0.00000e+00;AF_Female=0.00000e+00;AF_Male=0.00000e+00;AF_NFE=0.00000e+00;AF_OTH=0.00000e+00;AF_POPMAX=.;AF_SAS=0.00000e+00;AF_raw=9.99900e-06;AN=50778;AN_AFR=4986;AN_AMR=10892;AN_ASJ=1274;AN_EAS=7560;AN_FIN=694;AN_Female=24940;AN_Male=25838;AN_NFE=17556;AN_OTH=1486;AN_POPMAX=.;AN_SAS=6330;AN_raw=100010;AS_FilterStatus=RF|AC0;AS_RF=1.49748e-01;BaseQRankSum=-4.60000e-01;CSQ=-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034|YES|||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|6/6||ENST00000450305.2:n.561_562delTG||558-559||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.847_848delTG||844-845||||||1||1||deletion|1|HGNC|37102|YES||||||||||||3|||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene|||||||||||1|807|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|3/3||ENST00000515242.2:n.840_841delTG||837-838||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/4||ENST00000518655.2:n.678_679delTG||675-676||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene|||||||||||1|814|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site|||||||||||1||||deletion|1||||||||||||||||||||||||||||||||||||||||||||;ClippingRankSum=5.63000e-01;DP=2519792;DP_HIST_ALL=20921|3680|466|85|62|97|652|4365|4551|3656|2891|2039|1464|1114|954|811|688|497|352|310;DP_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;DP_MEDIAN=118;DREF_MEDIAN=3.98107e-38;FS=1.59250e+01;GC=25389,0,0;GC_AFR=2493,0,0;GC_AMR=5446,0,0;GC_ASJ=637,0,0;GC_EAS=3780,0,0;GC_FIN=347,0,0;GC_Female=12470,0,0;GC_Male=12919,0,0;GC_NFE=8778,0,0;GC_OTH=743,0,0;GC_SAS=3165,0,0;GC_raw=50004,1,0;GQ_HIST_ALL=11211|8535|2038|2055|803|203|195|95|28|49|65|37|115|64|88|117|164|34|237|23872;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;GQ_MEDIAN=99;Hom=0;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_Female=0;Hom_Male=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=0;Hom_raw=0;InbreedingCoeff=-4.37000e-02;MQ=3.15600e+01;MQRankSum=-8.97000e-01;POPMAX=.;QD=3.22000e+00;ReadPosRankSum=-1.23200e+00;SOR=1.09000e-01;VQSLOD=-1.83100e+00;VQSR_NEGATIVE_TRAIN_SITE;VQSR_culprit=QD;gnomad.exome.AC_AFR=0;gnomad.exome.AC_AMR=0;gnomad.exome.AC_ASJ=0;gnomad.exome.AC_EAS=0;gnomad.exome.AC_FIN=0;gnomad.exome.AC_Female=0;gnomad.exome.AC_Male=0;gnomad.exome.AC_NFE=0;gnomad.exome.AC_OTH=0;gnomad.exome.AC_raw=1;gnomad.exome.AN_AFR=4986;gnomad.exome.AN_AMR=10892;gnomad.exome.AN_ASJ=1274;gnomad.exome.AN_EAS=7560;gnomad.exome.AN_FIN=694;gnomad.exome.AN_Female=24940;gnomad.exome.AN_Male=25838;gnomad.exome.AN_NFE=17556;gnomad.exome.AN_OTH=1486;gnomad.exome.AN_raw=100010;gnomad.genome.AC_AFR=0;gnomad.genome.AC_AMR=0;gnomad.genome.AC_ASJ=0;gnomad.genome.AC_EAS=0;gnomad.genome.AC_FIN=0;gnomad.genome.AC_Female=0;gnomad.genome.AC_Male=0;gnomad.genome.AC_NFE=0;gnomad.genome.AC_OTH=0;gnomad.genome.AC_raw=1;gnomad.genome.AN_AFR=8680;gnomad.genome.AN_AMR=794;gnomad.genome.AN_ASJ=224;gnomad.genome.AN_EAS=1592;gnomad.genome.AN_FIN=3490;gnomad.genome.AN_Female=13274;gnomad.genome.AN_Male=16168;gnomad.genome.AN_NFE=13754;gnomad.genome.AN_OTH=908;gnomad.genome.AN_raw=30500

 
  ```
  
 END_DOC
 */
@Program(name="vcfgnomad",description="Peek annotations from gnomad",keywords={"vcf","annotation","gnomad"})
public class VcfGnomad extends Launcher{
	
	private static final Logger LOG = Logger.build(VcfGnomad.class).make();
	private static final int BUFFER_SIZE = 100000;
	/** allele specific population in gnomad */
	private final static String POPS[]=new String[]{"AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH", "Male", "Female", "raw", "POPMAX"}; 
	/** 'ome'-type section */
	private enum OmeType {exome,genome};
	/** entries mapping chromosome/type->vcf.gz */
	private List<ManifestEntry> manifestEntries=new ArrayList<>();
	
	
	private static class ManifestEntry
		implements Closeable
		{
		OmeType omeType;
		String contig;
		String uri;
		TabixVcfFileReader tabix=null;
		int buffferChromEnd=0;
		final List<VariantContext> buffer=new ArrayList<>(BUFFER_SIZE);
		@Override
		public void close() throws IOException {
			CloserUtil.close(tabix);
			buffer.clear();
			this.buffferChromEnd=0;
			tabix=null;
			}
		/** find matching variant in tabix file, use a buffer to avoid multiple random accesses */
		VariantContext findMatching(final String normContig,final int pos,Allele ref)
			{
			//past last buffer ? refill buffer
			if(this.buffferChromEnd<=pos)
				{
				buffer.clear();
				this.buffferChromEnd=pos+BUFFER_SIZE;
				final Iterator<VariantContext> iter=this.tabix.iterator(
						normContig, Math.max(0,pos-1),
						this.buffferChromEnd
						);
				while(iter.hasNext())
					{
					this.buffer.add(iter.next());
					}
				CloserUtil.close(iter);
				}
			
			for(VariantContext ctx:this.buffer)
				{
				if(!ctx.getContig().equals(normContig)) continue;
				if(ctx.getStart()<pos) continue;
				if(ctx.getStart()>pos) break;
				if(!ctx.getReference().equals(ref)) continue;
				return ctx;
				}
			return null;
			}
		
		}
	
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;
	@Parameter(names={"-m","--manifest"},description="manifest file descibing how to map a contig to an URI . 3 columns: 1) exome|genome 2) contig 3) path or URL.")
	private File manifestFile=null;
	@Parameter(names={"-filtered","--filtered"},description="Skip Filtered")
	private boolean skipFiltered=false;

	@Parameter(names={"-gf","--gnomadFilter"},description="if defined, add this FILTER when the variant is found in nomad")
	private String inGnomadFilterName=null;
	@Parameter(names={"-ac","--alleleconcordance"},description="ALL Alt allele must be found in gnomad before setting a FILTER")
	private boolean alleleconcordance=false;
	
	private class InfoField
		{
		final OmeType ome;
		final String tag;
		final List<Integer> attributes=new ArrayList<>();
		InfoField(String tag, OmeType ome) {
			this.tag=tag;
			this.ome=ome;
			}
		public String getOutputTag() {
			return "gnomad_"+ this.ome.name()+"_"+this.tag;
		}
		VCFInfoHeaderLine makeVCFInfoHeaderLine()
			{
			return new VCFInfoHeaderLine(
					getOutputTag(),VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer,
					"Field "+this.tag+" extracted from Gnomad ("+ome.name()+")"
					);
			}
		void fill(final VariantContext ctx,final VariantContext gnomadCtx)
			{
			this.attributes.clear();
			final List<Allele> galts=gnomadCtx.getAlternateAlleles();
			final List<String> gatts = gnomadCtx.getAttributeAsStringList(this.tag,null);
			for(final Allele a:ctx.getAlternateAlleles())
				{
				Integer found=null;
				//final int idx=gnomadCtx.getAlleleIndex(a);//non idx(REF)==0
				final int idx=galts.indexOf(a);
				
				if(idx>=0) {
					if(idx<gatts.size() && gatts.get(idx)!=null && !gatts.get(idx).equals(".")) {
						found=Integer.parseInt(gatts.get(idx));
						}
					}
				this.attributes.add(found);
				}
			}
		}
	
	private String normalizeContig(String contig)
		{
		if(contig.startsWith("chr")) contig=contig.substring(3);
		return contig;
		}
	
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator iter, VariantContextWriter out) {
		final ManifestEntry ome2manifest[]=new ManifestEntry[OmeType.values().length];
		Arrays.fill(ome2manifest,null);

		try {
			final List<InfoField> infoFields=new ArrayList<>();
			for(OmeType ome:OmeType.values()) {
				for(final String pop: POPS)
					{
					infoFields.add(new InfoField("AC_"+pop,ome));
					infoFields.add(new InfoField("AN_"+pop,ome));
					}
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
			while(iter.hasNext()) {
				final VariantContext ctx=iter.next();
				if(this.skipFiltered && ctx.isFiltered() )
					{
					out.add(ctx);
					continue;
					}
				final String ensemblContig=normalizeContig(ctx.getContig());

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
							newEntry.tabix=new TabixVcfFileReader(newEntry.uri);
							}
						}
					}
				
				
				for(final InfoField infoField: infoFields)
					{
					infoField.attributes.clear();
					}
				
				boolean setfilter=false;
				
				// lopp over exome and genome data
				for(int i=0;i< ome2manifest.length;++i) {
					ManifestEntry entry = ome2manifest[i];
					if(entry==null) continue;
					
					final VariantContext ctx2=entry.findMatching(ensemblContig, ctx.getStart(), ctx.getReference());
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
			return 0;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		}
		finally {
			CloserUtil.close(Arrays.asList(ome2manifest));
		}
	}

@Override
public int doWork(List<String> args) {
	if(this.manifestFile==null)
		{
		LOG.info("Building default manifest file...");
		for(OmeType ot: OmeType.values()) {
			for(int i=1;i<= 23;++i) {
				ManifestEntry entry=new ManifestEntry();
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
		}
	else
		{
		try {
		Files.lines(this.manifestFile.toPath()).forEach(L->{
			if(L.startsWith("#") || L.trim().isEmpty()) return;
			final String tokens[]=L.split("[\t]");
			if(tokens.length<3) throw new JvarkitException.TokenErrors("Expected 3 words",tokens);
			ManifestEntry entry=new ManifestEntry();
			entry.omeType=OmeType.valueOf(tokens[0]);
			entry.contig = tokens[1].trim();
			entry.uri=tokens[2].trim();
			VcfGnomad.this.manifestEntries.add(entry);
		});
		} catch(IOException err) {
			LOG.error(err);
			return -1;
		}
		}
	return doVcfToVcf(args, outputFile);
	}
	
public static void main(String[] args) {
	new VcfGnomad().instanceMainWithExit(args);
	}
}
