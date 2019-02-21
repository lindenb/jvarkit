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
package com.github.lindenb.jvarkit.tools.gnomadpext;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
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
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Context

https://www.biorxiv.org/content/10.1101/554444v1

> "Transcript expression-aware annotation improves rare variant discovery and interpretation" 
> Here, we develop a transcript-level annotation metric, the proportion expressed across transcripts (pext), which summarizes isoform quantifications for variants. We calculate this metric using 11,706 tissue samples from the Genotype Tissue Expression project (GTEx) and show that it clearly differentiates between weakly and highly evolutionarily conserved exons, a proxy for functional importance. 

## Example

```
# bgzip if needed and index the database with tabix
$ bgzip data.tsv
$ tabix -f -b 2 -e 2 -s 1 -c 'c' data.tsv.gz

# annotate
$ java -jar dist/vcfgnomadpext.jar -d data.tsv.gz input.vcf.gz

(...)
1	135803	.	C	G	.	.	AC=8;AF=0.500;AN=16;DP=113;GNOMAD_PEXT=ensg:ENSG00000237683|csq:splice_acceptor_variant|symbol:AL627309.1|lof:LC|lof_flag:NON_CAN_SPLICE|Spleen:1.0|Brain_FrontalCortex_BA9_:1.0|SmallIntestine_TerminalIleum:1.0|Skin_SunExposed_Lowerleg_:1.0|Artery_Coronary:1.0|Brain_Hippocampus:1.0|Esophagus_Muscularis:1.0|Brain_Nucleusaccumbens_basalganglia_:1.0|Artery_Tibial:1.0|Brain_Hypothalamus:1.0|Adipose_Visceral_Omentum_:1.0|Nerve_Tibial:1.0|Brain_CerebellarHemisphere:1.0|Breast_MammaryTissue:1.0|Liver:1.0|Skin_NotSunExposed_Suprapubic_:1.0|AdrenalGland:1.0|Pancreas:1.0|Lung:1.0|Pituitary:1.0|Muscle_Skeletal:1.0|Colon_Transverse:1.0|Artery_Aorta:1.0|Heart_AtrialAppendage:1.0|Adipose_Subcutaneous:1.0|Esophagus_Mucosa:1.0|Heart_LeftVentricle:1.0|Brain_Cerebellum:1.0|Brain_Cortex:1.0|Thyroid:1.0|Stomach:1.0|WholeBlood:1.0|Brain_Anteriorcingulatecortex_BA24_:1.0|Brain_Putamen_basalganglia_:1.0|Brain_Caudate_basalganglia_:1.0|Colon_Sigmoid:1.0|Esophagus_GastroesophagealJunction:1.0|Brain_Amygdala:1.0|mean_proportion:1.0	GT:DP	1/10/0:39	1/1	1/0:25	0/0:11	1/0	0/1:38	0/1
1	135804	.	G	A	.	.	AC=6;AF=0.375;AN=16;DP=82;GNOMAD_PEXT=ensg:ENSG00000237683|csq:splice_acceptor_variant|symbol:AL627309.1|lof:LC|lof_flag:NON_CAN_SPLICE|Spleen:1.0|Brain_FrontalCortex_BA9_:1.0|SmallIntestine_TerminalIleum:1.0|Skin_SunExposed_Lowerleg_:1.0|Artery_Coronary:1.0|Brain_Hippocampus:1.0|Esophagus_Muscularis:1.0|Brain_Nucleusaccumbens_basalganglia_:1.0|Artery_Tibial:1.0|Brain_Hypothalamus:1.0|Adipose_Visceral_Omentum_:1.0|Nerve_Tibial:1.0|Brain_CerebellarHemisphere:1.0|Breast_MammaryTissue:1.0|Liver:1.0|Skin_NotSunExposed_Suprapubic_:1.0|AdrenalGland:1.0|Pancreas:1.0|Lung:1.0|Pituitary:1.0|Muscle_Skeletal:1.0|Colon_Transverse:1.0|Artery_Aorta:1.0|Heart_AtrialAppendage:1.0|Adipose_Subcutaneous:1.0|Esophagus_Mucosa:1.0|Heart_LeftVentricle:1.0|Brain_Cerebellum:1.0|Brain_Cortex:1.0|Thyroid:1.0|Stomach:1.0|WholeBlood:1.0|Brain_Anteriorcingulatecortex_BA24_:1.0|Brain_Putamen_basalganglia_:1.0|Brain_Caudate_basalganglia_:1.0|Colon_Sigmoid:1.0|Esophagus_GastroesophagealJunction:1.0|Brain_Amygdala:1.0|mean_proportion:1.0	GT:DP	0/00/0	0/1:15	0/0:5	1/1:15	0/1:17	0/0:8	1/1:22
1	137619	.	G	A	.	.	AC=7;AF=0.438;AN=16;DP=193;GNOMAD_PEXT=ensg:ENSG00000237683|csq:splice_donor_variant|symbol:AL627309.1|lof:LC|lof_flag:NON_CAN_SPLICE|Spleen:1.0|Brain_FrontalCortex_BA9_:1.0|SmallIntestine_TerminalIleum:1.0|Skin_SunExposed_Lowerleg_:1.0|Artery_Coronary:1.0|Brain_Hippocampus:1.0|Esophagus_Muscularis:1.0|Brain_Nucleusaccumbens_basalganglia_:1.0|Artery_Tibial:1.0|Brain_Hypothalamus:1.0|Adipose_Visceral_Omentum_:1.0|Nerve_Tibial:1.0|Brain_CerebellarHemisphere:1.0|Breast_MammaryTissue:1.0|Liver:1.0|Skin_NotSunExposed_Suprapubic_:1.0|AdrenalGland:1.0|Pancreas:1.0|Lung:1.0|Pituitary:1.0|Muscle_Skeletal:1.0|Colon_Transverse:1.0|Artery_Aorta:1.0|Heart_AtrialAppendage:1.0|Adipose_Subcutaneous:1.0|Esophagus_Mucosa:1.0|Heart_LeftVentricle:1.0|Brain_Cerebellum:1.0|Brain_Cortex:1.0|Thyroid:1.0|Stomach:1.0|WholeBlood:1.0|Brain_Anteriorcingulatecortex_BA24_:1.0|Brain_Putamen_basalganglia_:1.0|Brain_Caudate_basalganglia_:1.0|Colon_Sigmoid:1.0|Esophagus_GastroesophagealJunction:1.0|Brain_Amygdala:1.0|mean_proportion:1.0	GT:DP	0/1:25	0/1:10	0/1:46	0/0:32	0/0	1/0	1/1:31	0/1:49
```



END_DOC
 */
@Program(name="vcfgnomadpext",
	description="Peek annotations from gnomadpext",
	keywords={"vcf","annotation","gnomad"},
	creationDate="20190220",
	modificationDate="20190220"
	)
public class VcfGnomadPext extends Launcher{
	
	private static final Logger LOG = Logger.build(VcfGnomadPext.class).make();
	

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--database","--pext"},description="Pext file. tab delimited :(chrom\\tpos\\tref\\talt\\ttx_annotation). Bgziped and indexed with tabix.",required=true)
	private String pextDatabasePath=null;
	@Parameter(names={"--bufferSize"},description="When we're looking for variant in Gnomad, load the variants for 'N' bases instead of doing a random access for each variant. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int gnomadBufferSize= 10_000;
	@Parameter(names={"-filtered","--filtered"},description="Skip Filtered User Variants")
	private boolean skipFiltered=false;
	@Parameter(names={"-T","--tissues"},description="Restrict to those tissues.")
	private Set<String> restrictTissues=new HashSet<>();

	
	
	private class PextEntry
		implements Locatable
		{
		final String contig;
		final int pos;
		final Allele ref;
		final Allele alt;
		final String jsonStr;
		PextEntry(final String tokens[])
			{
			this.contig = tokens[0];
			this.pos = Integer.parseInt(tokens[1]);
			this.ref = Allele.create(tokens[2], true);
			this.alt = Allele.create(tokens[3], false);
			this.jsonStr = tokens[4];
			}
		
		@Override
		public String getContig() {
			return this.contig;
			}
		@Override
		public int getStart() {
			return this.pos;
			}
		@Override
		public int getEnd() {
			return this.pos;
			}
		@Override
		public String toString() {
			return contig+":"+pos+"("+ref.getDisplayString()+"/"+alt.getDisplayString()+")="+this.jsonStr;
			}
		}
		
		
	private final ContigNameConverter ensemblCtgConvert = ContigNameConverter.createConvertToEnsembl();
	private final List<PextEntry> buffer = new ArrayList<>();
	private Interval lastInterval = null;
	
	/** find matching variant in tabix file, use a buffer to avoid multiple random accesses */
	final List<PextEntry> findOverlapping(final TabixFileReader tabix,final VariantContext ctx)
		{
		
		final String normContig = this.ensemblCtgConvert.apply(ctx.getContig());
		if(StringUtil.isBlank(normContig)) return Collections.emptyList();
		
		if(!buffer.isEmpty() && !buffer.get(0).contig.equals(normContig)) {
			this.buffer.clear();
			}
		
		
		if(this.lastInterval==null ||
			!this.lastInterval.getContig().equals(normContig) ||
			!CoordMath.encloses(lastInterval.getStart(), lastInterval.getEnd(), ctx.getStart(), ctx.getEnd())
			)
			{
			final CharSplitter tab = CharSplitter.TAB;
			this.buffer.clear();
			this.lastInterval = new Interval(
					normContig,
					Math.max(0, ctx.getStart()-10),
					ctx.getEnd()+ VcfGnomadPext.this.gnomadBufferSize
				);
			final Iterator<String> iter= tabix.iterator(
					this.lastInterval.getContig(),
					this.lastInterval.getStart(),
					this.lastInterval.getEnd()
					);
			while(iter.hasNext())
				{
				final String tokens[]=tab.split(iter.next());
				this.buffer.add(new PextEntry(tokens));
				}
			CloserUtil.close(iter);
			}
		
		return this.buffer.stream().
				filter(V->
						V.getStart()==ctx.getStart() && 
						V.getEnd()==ctx.getEnd() &&
						V.getContig().equals(normContig) && 
						V.ref.equals(ctx.getReference()) &&
						ctx.getAlleles().contains(V.alt)).
				collect(Collectors.toList());
		}
		
		
		
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter out
			)
		{
		final JsonParser jsonParser = new JsonParser();

		final String standard_pext_header[]=new String[] {"chrom","pos","ref","alt","tx_annotation"};
		final VCFHeader h0 = iter.getHeader();
		if(!SequenceDictionaryUtils.isGRCh37(h0)) {
			LOG.error("Input is NOT GRCh37 ");
			return -1;
			}
		
		TabixFileReader gextFileReader = null;
		final CharSplitter tab =  CharSplitter.TAB;
		try {
			gextFileReader  = new TabixFileReader(this.pextDatabasePath);
			final String line1= gextFileReader.readLine();
			if(StringUtils.isBlank(line1)) {
				LOG.error("Cannot read first line of "+this.pextDatabasePath);
				return -1;
				}
			if(!Arrays.equals(tab.split(line1),standard_pext_header))
				{
				LOG.error("Bad header in "+this.pextDatabasePath+" expected " +
						String.join("(tab)",standard_pext_header)+" but got " +
						line1.replace("\t", "(tab)"));
				return -1;
				}
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(h0).logger(LOG).build();
		
			final VCFHeader h2 = new VCFHeader(h0);
			final VCFInfoHeaderLine pexInfo = new VCFInfoHeaderLine(
					"GNOMAD_PEXT",
					VCFHeaderLineCount.A,
					VCFHeaderLineType.String,
					"Gnomad Data from "+this.pextDatabasePath
					);
			h2.addMetaDataLine(pexInfo);
			
			JVarkitVersion.getInstance().addMetaData(this, h2);
			out.writeHeader(h2);
		
			while(iter.hasNext()) {
				final VariantContext ctx = progress.apply(iter.next());
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				
				vcb.rmAttribute(pexInfo.getID());
				
				
				if(!ctx.isVariant() || (this.skipFiltered && ctx.isFiltered()) )
					{
					out.add(vcb.make());
					continue;
					}
				
				final List<PextEntry> entries = findOverlapping(gextFileReader,ctx);
				if(entries.isEmpty()) {
					out.add(vcb.make());
					continue;
					}
				
			
				final List<String> altInfo=new ArrayList<>(ctx.getAlleles().size());
				for(int allele_idx=1/* 0 is ref */;allele_idx< ctx.getAlleles().size();allele_idx++)
					{
					final Allele alt = ctx.getAlleles().get(allele_idx);
					final PextEntry entry = entries.stream().filter(E->E.alt.equals(alt)).findFirst().orElse(null);
					if(entry==null)
						{
						altInfo.add(".");
						}
					else
						{
						final JsonElement e= jsonParser.parse(entry.jsonStr);
						if(!e.isJsonArray()) throw new IllegalStateException("not an array: "+entry.jsonStr);
						final JsonArray array=e.getAsJsonArray();
						final StringBuilder sb=new StringBuilder();
						for(int x=0;x<array.size();++x)
							{
							if(x>0) sb.append("&");
							final StringBuilder sb2=new StringBuilder();
							final JsonObject obj = array.get(x).getAsJsonObject();
							for(final Map.Entry<String,JsonElement> kv:obj.entrySet())
								{
								String key=kv.getKey();
								// "Brain_FrontalCortex_BA9_": 1.0,
								if(key.endsWith("_")) key=key.substring(0,key.length()-1);
								// as far as I can see, tissues start with a uppercase
								if(Character.isUpperCase(key.charAt(0)) && 
									!this.restrictTissues.isEmpty() &&
									!this.restrictTissues.contains(key)
									) continue;
								final JsonElement v = kv.getValue();
								if(v.isJsonNull()) continue;
								if(v.getAsJsonPrimitive().isString()) {
									final String strv = v.getAsString();
									if(strv.equals("NaN")) continue;
									}
								
								if(sb2.length()>0) sb2.append("|");
								sb2.append(key);
								sb2.append(":");
								sb2.append(kv.getValue().getAsString());
								}
							
							sb.append(sb2.toString());
							}
						if(sb.length()==0) sb.append(".");
						altInfo.add(sb.toString());
						}
					}
				//at least none is not '.'
				if(altInfo.stream().anyMatch(S->!S.equals("."))) {
					vcb.attribute(pexInfo.getID(), altInfo);
				}
			
				out.add(vcb.make());
				}
			out.close();
			progress.close();
			gextFileReader.close();gextFileReader=null;
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(gextFileReader);
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.gnomadBufferSize < 10) {
			LOG.error("buffer size is too small "+this.gnomadBufferSize);
			return -1;
			}
		return doVcfToVcf(args,this.outputFile);
		}

public static void main(final String[] args) {
	new VcfGnomadPext().instanceMainWithExit(args);
	}
}
