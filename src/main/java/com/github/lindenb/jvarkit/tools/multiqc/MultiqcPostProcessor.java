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
package com.github.lindenb.jvarkit.tools.multiqc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.util.Maps;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

/**
BEGIN_DOC

## Motivation

Enhances multiqc output by reading the data folder and producing new plots (eg. boxplot per population

If no group is defined and the tool can find the file `dragen_ploidy.txt`, this file is used to create a group with male and females

Doesn't work for now ? https://github.com/MultiQC/MultiQC/issues/2689

## Example

```
# run first time
multiqc --force --file-list input.list
java -jar jvarkit.jar multiqcpostproc --sample2collection sample2collection.tsv multiqc_data -o OUTDIR2
find OUTDIR2 --type f -name "*.json" >> input.list
# run second time
multiqc --force --file-list input.list
```

END_DOC
*/
@Program(name="multiqcpostproc",
description="Enhances multiqc output by reading the data folder and producing new plots (eg. boxplot per population.",
keywords={"multiqc"},
creationDate="20240708",
modificationDate="20240730",
jvarkit_amalgamion = true,
menu="Utilities"
)
public class MultiqcPostProcessor extends Launcher {
	private static final Logger LOG = Logger.of(MultiqcPostProcessor.class);
	private final char GE_UNICODE ='\u2265';
	@Parameter(names={"-o","--output"},description="output directory",required = true)
	private Path outputDirectory = null;
	@Parameter(names={"--sample2collection"},description=SamplePopulation.OPT_DESC)
	private Path sample2collectionPath=null;
	@Parameter(names={"--custom"},description="custom mapping file (undocumented")
	private Path customMapping=null;
	@Parameter(names={"--beeswarm"},description="use plot_type=beeswarm instead of boxplot.")
	private boolean use_beeswarm = false;
	@Parameter(names={"--title"},description="main section title")
	private String main_section_title = "";
	@Parameter(names={"--description"},description="main section description")
	private String main_section_description = "";


	private final SamplePopulation sampleCollection = new SamplePopulation();
	
	private static class SampleValue {
		final String sample;
		final double value;
		SampleValue(final String sample,double value) {
			this.sample = sample;
			this.value = value;
		}
	}
	
	private static class FileContent {
		FileHeader fileHeader;
		List<FileHeader.RowMap> rows;
	}
	
	private abstract class Handler {
		final Map<String,String> properties=new HashMap<>();
		Handler(final Map<String,String> prop) {
			this.properties.putAll(prop);
			}
		
		public boolean isHandlerForFile(String fname) {
			final String[] s=  CharSplitter.COMMA.split(this.properties.getOrDefault("filename", ""));
			return Arrays.stream(s).
					map(S->{
						if(S.startsWith("glob:")) {
							final PathMatcher pm = FileSystems.getDefault().getPathMatcher(S);
							return new Predicate<String>() {
								public boolean test(final String s) {
									boolean b= pm.matches(Paths.get(s).getFileName());
									LOG.debug("glob matcher for "+s+"/"+S+" is "+b);
									return b;
									}
								};
							}
						else
							{
							final String x1=S;
							return new Predicate<String>() {
								public boolean test(final String s) {
									return x1.equals(s);
									}
								};
							}
						}).
						anyMatch(S->S.test(fname)
						);
			}
		protected String fixSampleName(final String sn) {
			final String idxstats=".idxstat";
			if(sn.endsWith(idxstats)) {
				return fixSampleName(sn.substring(0,sn.length()-idxstats.length()));
				}
			return sn;
			}
		protected Set<String> getDataColumns() {
			final String s = this.properties.getOrDefault("data", "");
			return Arrays.asList(CharSplitter.COMMA.split(s)).stream().filter(S->!StringUtils.isBlank(S)).collect(Collectors.toSet());
			}
		protected String getSectionName(String title) {
			return this.properties.getOrDefault("section_name", "${title} per population").replace("${title}",title);
			}
		protected String getSectionDescription(String title) {
			return this.properties.getOrDefault("section_description",getSectionName(title)).replace("${title}",title);
			}
		
		protected void saveAs(final Set<Path> filetSet,Collection<SampleValue> sample_values,final String titleRaw) throws IOException {
			if(sample_values.isEmpty()) return;

			//remove unicode
			final String title = titleRaw.replaceAll(String.valueOf(GE_UNICODE),">=");

			
			// make a copy of SamplePopulation
			final SamplePopulation snpop =new SamplePopulation(MultiqcPostProcessor.this.sampleCollection);
			snpop.retainSamples(sample_values.stream().
					map(ROW-> ROW.sample).
					collect(Collectors.toSet())
					);
			
			/* is it worth plotting ?*/
			boolean found_deviation=false;
			final Map<String,List<SampleValue>> pop2data=new TreeMap<>();
			
			for(final SamplePopulation.Population pop:snpop.getPopulations()) {
				final List<SampleValue> values= sample_values.stream().
					filter(SV->pop.containsKey(SV.sample)).
					sorted((A,B)->Double.compare(A.value, B.value)).
					collect(Collectors.toList());
				/* is it worth plotting ?*/
				if(values.size()>1 && values.get(0).value < values.get(values.size()-1).value) {
					found_deviation=true;
					}
				pop2data.put(pop.getName(), values);
				}
			if(pop2data.isEmpty() || !found_deviation) {
				LOG.info("no data to plot for "+filetSet+"/"+title);
				return;
				}
			
			final Path inputFile = filetSet.iterator().next();
			final Path outfile = outputDirectory.resolve(IOUtils.getFilenameWithoutCommonSuffixes(inputFile)+"_"+ title.replaceAll("[^A-Za-z0-9_]+", "_").toLowerCase()+"_mqc.json");
			LOG.info("writing "+outfile);
			try(PrintWriter pw = IOUtils.openPathForPrintWriter(outfile)) {
				final String id = StringUtils.md5( inputFile.getFileName().toString()+title);
				try(JsonWriter w = new JsonWriter(pw)) {
					w.beginObject();
					
					w.name("parent_id");
					w.value("p"+StringUtils.md5(StringUtils.ifBlank(MultiqcPostProcessor.this.main_section_title,MultiqcPostProcessor.class.getName())));
					w.name("parent_name");
					w.value(StringUtils.ifBlank(MultiqcPostProcessor.this.main_section_title,"Per Population"));
					w.name("parent_description");
					w.value(StringUtils.ifBlank(MultiqcPostProcessor.this.main_section_description,"MULTIQC data grouped by population using jvarkit "+MultiqcPostProcessor.class.getSimpleName()));
					
					w.name("id");
					w.value("section_id"+id);
					w.name("section_name");
					w.value(getSectionName(title));
					w.name("section_description");
					w.value(getSectionDescription(title)+" from "+inputFile.getFileName());
					
					
					w.name("plot_type");
					w.value(use_beeswarm?"beeswarm":"box");
					w.name("pconfig");
					w.beginObject();
					w.name("id");
					w.value("plot__id"+id);
					w.name("title");
					w.value(title+ " per population");
					w.name("xlab");
					w.value(title);
					w.name("ylab");
					w.value("collection");
					if(use_beeswarm) {
						w.name("xmin");
						w.value(pop2data.values().stream().flatMap(it->it.stream()).mapToDouble(X->X.value).min().orElse(0.0));
						w.name("xmax");
						w.value(pop2data.values().stream().flatMap(it->it.stream()).mapToDouble(X->X.value).max().orElse(1.0));
					}
					
					w.endObject();
					
					w.name("data");
					w.beginObject();
					if(use_beeswarm) {
						for(final String sn:pop2data.values().stream().flatMap(L->L.stream()).map(K->K.sample).collect(Collectors.toSet())) {
							w.name(sn);
							w.beginObject();
							for(final String popName:pop2data.keySet()) {
								final SampleValue snv = pop2data.get(popName).stream().filter(it->it.sample.equals(sn)).findFirst().orElse(null);
								if(snv==null) continue;
								final SamplePopulation.Population pop = snpop.getPopulationByName(popName);
								w.name(pop.getName()+" (N="+pop.size()+")");
								w.value(snv.value);
								}
							w.endObject();
							}
						}
					else
						{
						for(final String popName:pop2data.keySet()) {
							final SamplePopulation.Population pop = snpop.getPopulationByName(popName);
							w.name(pop.getName()+" (N="+pop.size()+")");
							w.beginArray();
							for(SampleValue sv: pop2data.get(popName)) {
								w.value(sv.value);
								}
							w.endArray();
							}
						}
					w.endObject();
					
					w.endObject();
					w.flush();
					}
				pw.flush();
				}
			}
		
		public abstract void apply(final Set<Path> pathSet);
		};
	
	private abstract class AbstractHandler extends Handler {
		AbstractHandler(final Map<String,String> prop) {
			super(prop);
			}
		
		protected String getSampleColumn() {
			return this.properties.getOrDefault("sample", "Sample");
			}
		
	
		
		/** for some tools multiqc normalize by x. Eg: samtools stats results is 'x' in MB (0.0001 ) */
		protected double getFactor() {
			return Double.parseDouble(this.properties.getOrDefault("factor", "1.0"));
			}
		}
	
	private class JsonHandler extends Handler {
		JsonHandler(final Map<String,String> prop) {
			super(prop);
			}
		private void scan(JsonElement root,String dataCol,final Map<String,SampleValue> sample_values) {
			if(root.isJsonArray()) {
				JsonArray L= root.getAsJsonArray();
				for(int i=0;i< L.size();++i) {
					scan(L.get(i),dataCol,sample_values);
					}
				}
			else if(root.isJsonObject()) {
				final JsonObject hash = root.getAsJsonObject();
				boolean found=false;
				for(Map.Entry<String, JsonElement> kv:hash.entrySet()) {
					final JsonElement jsE = hash.get(kv.getKey());
					if(!jsE.isJsonObject()) continue;
					final JsonObject hash2 = jsE.getAsJsonObject();
					if(!hash2.has(dataCol)) continue;
					final JsonElement jsV = hash2.get(dataCol);
					if(!jsV.isJsonPrimitive()) continue;
					if(!jsV.getAsJsonPrimitive().isNumber()) continue;
					final double v;
					try {
						v= jsV.getAsJsonPrimitive().getAsDouble();
						}
					catch(Throwable err) {
						continue;
						}
					final String sn = fixSampleName(kv.getKey());
					if(!MultiqcPostProcessor.this.sampleCollection.hasSample(sn)) continue;
					if(sample_values.containsKey(sn)) {
						LOG.warn("duplicate key for sample "+sn+" "+dataCol);
						continue;
						}
					sample_values.put(sn,new SampleValue(sn, v));
					found=true;
					}
				if(!found) {
					for(Map.Entry<String, JsonElement> kv:hash.entrySet()) {
						scan(kv.getValue(),dataCol,sample_values);
						}
					}
				
				}
			}
		
		
		@Override
		public void apply(final Set<Path> pathSet) {
			try {
				final JsonParser parser=new JsonParser();

				for(final String dataCol: getDataColumns()) {
					final Map<String,SampleValue> sample_values=new HashMap<>();
					for(final Path f:pathSet) {
						final JsonElement root;
						try(Reader r=Files.newBufferedReader(f)) {
							root = parser.parse(r);
							}
						scan(root,dataCol,sample_values);
						}
					if(sample_values.isEmpty()) {
						LOG.warn("no data "+dataCol+" in "+pathSet);
						continue;
						}
					saveAs(pathSet,sample_values.values(),dataCol);
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				}
			}
		}
	
	@SuppressWarnings("unused")
	private class DistributionHandler extends AbstractHandler {
		DistributionHandler(final Map<String,String> prop) {
			super(prop);
			}
		
		private OptionalDouble parseKeyValue(String s,final String key) {
			if(StringUtils.isBlank(s)) return OptionalDouble.empty();
			
			try {
				if(!(s.startsWith("(") && s.endsWith(")") && s.contains(","))) return OptionalDouble.empty();
				s=s.substring(1, s.length()-1);
				final String[] tokens = CharSplitter.COMMA.split(s);
				if(tokens.length!=2) return OptionalDouble.empty();
				tokens[0]= tokens[0].trim();
				if(StringUtils.isBlank(tokens[0])) return OptionalDouble.empty();
				if(!tokens[0].equals(key)) return OptionalDouble.empty();
				tokens[1]= tokens[1].trim();
				if(StringUtils.isBlank(tokens[1])) return OptionalDouble.empty();
				final double v= Double.parseDouble(tokens[1]);
				if(Double.isNaN(v)) return OptionalDouble.empty();
				if(Double.isInfinite(v)) return OptionalDouble.empty();
				return OptionalDouble.of(v);
			} catch(final NumberFormatException f) {
				return OptionalDouble.empty();
			}
		}
		@Override
		public void apply(final Set<Path> pathSet) {
			try {
				
				final int sampleColumn =0;
				
				final double factor = getFactor();
				for(final String dataCol:  getDataColumns()) {
					final Map<String,SampleValue> sample_values=new HashMap<>();
					for(Path f: pathSet) {
						final FileContent fc = readFileContent(f);
						fc.fileHeader.assertColumn(getSampleColumn(),sampleColumn);
						
						fc.rows.stream().
							map(ROW->{
								final List<String> L=ROW.asList();
								final String sn= fixSampleName(L.get(sampleColumn));
								final OptionalDouble v = L.subList(1,L.size()).
										stream().
										map(it->parseKeyValue(it,dataCol)).
										filter(it->it.isPresent()).
										mapToDouble(it->it.getAsDouble()).
										findFirst();
								if(!v.isPresent()) return null;
										
								return new SampleValue(sn,v.getAsDouble()*factor);
								}
							).
							filter(SV->SV!=null).
							filter(SV->MultiqcPostProcessor.this.sampleCollection.hasSample(SV.sample)).
							forEach(SV->sample_values.put(SV.sample, SV));
						}
					saveAs(pathSet,sample_values.values(),dataCol);
					} //end loop over data columns
				}
			catch(final IOException err) {
				LOG.error(err);
				}
			}
		
	}
	
	private class BoxPlotHandler extends AbstractHandler {
		BoxPlotHandler(Map<String,String> prop) {
			super(prop);
			}
		
		
		private OptionalDouble parseDouble(final String s) {
			if(StringUtils.isBlank(s)) return OptionalDouble.empty();
			
			try {
				final double v= Double.parseDouble(s);
				if(Double.isNaN(v)) return OptionalDouble.empty();
				if(Double.isInfinite(v)) return OptionalDouble.empty();
				return OptionalDouble.of(v);
			} catch(final NumberFormatException f) {
				return OptionalDouble.empty();
			}
		}
		@Override
		public void apply(final Set<Path> pathSet) {
			try {
				final double factor = getFactor();
				for(final String dataCol: getDataColumns()) {
					if(StringUtils.isBlank(dataCol)) continue;
					final Map<String,SampleValue> sample_values=new HashMap<>();
					
					for(Path f: pathSet) {
						final FileContent fc = readFileContent(f);
						fc.fileHeader.assertColumnExists(getSampleColumn());
						
						if(!fc.fileHeader.containsKey(dataCol)) {
							LOG.warn("no column "+dataCol+" in "+f+" "+fc.fileHeader);
							continue;
							}
						fc.fileHeader.assertColumnExists(dataCol);
	
						fc.rows.stream().
							filter(ROW->parseDouble(ROW.get(dataCol)).isPresent()).
							map(ROW->new SampleValue(fixSampleName(ROW.get(getSampleColumn())), factor * parseDouble(ROW.get(dataCol)).getAsDouble())).
							filter(SV->MultiqcPostProcessor.this.sampleCollection.hasSample(SV.sample)).
							forEach(SV->sample_values.put(SV.sample, SV));
						}
					saveAs(pathSet,sample_values.values(),dataCol);
					} //end loop over data columns
				}
			catch(final IOException err) {
				LOG.error(err);
				}
			}
		}
		
	
	private FileContent readFileContent(final Path f) throws IOException {
		final FileContent fc = new FileContent();
		try(BufferedReader br = IOUtils.openPathForBufferedReading(f)) {
			String line=br.readLine();
			if(line==null) throw new IOException("cannot read first line of "+f);
			fc.fileHeader = new FileHeader(line, CharSplitter.TAB);
			fc.rows=fc.fileHeader.readAll(br);
			return fc;
			}
		
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			
			
			
			final List<Handler> handlers = new ArrayList<>();
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","singletons",
					"filename",String.join(",",
							"bcftools-stats-singletons.txt",
							"mqc_bcftools-stats-singletons__STDIN_.txt"
							)
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","depth",
					"filename",String.join(",",
							"bcftools-stats-sequencing-depth.txt",
							"glob:mqc_bcftools-stats-depth*.txt"
							)
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","nSNPs,nIndels",
					"filename",String.join(",",
							"bcftools-stats-sites.txt",
							"glob:mqc_bcftools-stats-sites*.txt"
							)
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","tstv",
					"filename",String.join(",",
							"bcftools-stats-tstv.txt",
							"glob:mqc_bcftools-stats-tstv*.txt"
							)
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","Total sequences,Mapped &amp; paired,Properly paired,Duplicated,QC Failed,Reads MQ0,Mapped bases (CIGAR),Bases Trimmed,Duplicated bases,Diff chromosomes,Other orientation,Inward pairs,Outward pairs",
					"factor","1000000",
					"filename","samtools-stats-dp.txt"
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","PCT of genome with coverage [  20x: inf),PCT of genome with coverage [  50x: inf),Average chr X coverage over genome,Average chr Y coverage over genome,Average mitochondrial coverage over genome,Average autosomal coverage over genome,Median autosomal coverage over genome,Mean/Median autosomal coverage ratio over genome,Aligned reads in genome pct",
					"filename","dragen_cov_metrics.txt"
					)));
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "Total","Biallelic pct","Multiallelic pct","Ti/Tv ratio","Percent Callability"),
					"filename","dragen_vc_metrics.txt"
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "Variants","Biallelic","Multiallelic","SNP","Indel","Ti/Tv","Callability","Autosome callability"),
					"filename","dragen-vc-metrics-table.txt"
					)));
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","Total input reads,Total input reads pct,Number of duplicate marked reads,Number of duplicate marked reads pct,Number of unique reads (excl. duplicate marked reads) pct,Reads with mate sequenced,Reads with mate sequenced pct,Reads without mate sequenced,Reads without mate sequenced pct,QC-failed reads,QC-failed reads pct,Mapped reads,Mapped reads pct,Mapped reads adjusted for filtered mapping,Mapped reads adjusted for filtered mapping pct,Mapped reads R1,Mapped reads R1 pct,Mapped reads R2,Mapped reads R2 pct,Number of unique & mapped reads (excl. duplicate marked reads),Number of unique & mapped reads (excl. duplicate marked reads) pct,Unmapped reads,Unmapped reads pct,Unmapped reads adjusted for filtered mapping,Unmapped reads adjusted for filtered mapping pct,Adjustment of reads matching non-reference decoys,Adjustment of reads matching non-reference decoys pct,Singleton reads (itself mapped; mate unmapped),Singleton reads (itself mapped; mate unmapped) pct,Paired reads (itself & mate mapped),Paired reads (itself & mate mapped) pct,Properly paired reads,Properly paired reads pct,Not properly paired reads (discordant),Not properly paired reads (discordant) pct,Paired reads mapped to different chromosomes,Paired reads mapped to different chromosomes pct,Paired reads mapped to different chromosomes (MAPQ>=10),Paired reads mapped to different chromosomes (MAPQ>=10) pct,Reads with MAPQ [40:inf),Reads with MAPQ [40:inf) pct,Reads with MAPQ [30:40),Reads with MAPQ [30:40) pct,Reads with MAPQ [20:30),Reads with MAPQ [20:30) pct,Reads with MAPQ [10:20),Reads with MAPQ [10:20) pct,Reads with MAPQ [ 0:10),Reads with MAPQ [ 0:10) pct,Reads with MAPQ NA (Unmapped reads),Reads with MAPQ NA (Unmapped reads) pct,Reads with indel R1,Reads with indel R1 pct,Reads with indel R2,Reads with indel R2 pct,Total bases,Total bases R1,Total bases R2,Mapped bases,Mapped bases R1,Mapped bases R2,Soft-clipped bases,Soft-clipped bases pct,Soft-clipped bases R1,Soft-clipped bases R1 pct,Soft-clipped bases R2 pct,Hard-clipped bases,Hard-clipped bases pct,Hard-clipped bases R1,Hard-clipped bases R1 pct,Hard-clipped bases R2,Hard-clipped bases R2 pct,Mismatched bases R1,Mismatched bases R1 pct,Mismatched bases R2 pct,Mismatched bases R1 (excl. indels),Mismatched bases R1 (excl. indels) pct,Mismatched bases R2 (excl. indels) pct,Q30 bases,Q30 bases pct,Q30 bases R1,Q30 bases R1 pct,Q30 bases R2,Q30 bases R2 pct,Q30 bases (excl. dups & clipped bases),Total alignments,Secondary alignments,Supplementary (chimeric) alignments,Estimated read length,Bases in reference genome,Insert length: mean,Insert length: median,Insert length: standard deviation,DRAGEN mapping rate [mil. reads/second],Secondary alignments pct,Q30 bases (excl. dups & clipped bases) pct,Mapped bases R1 pct,Mapped bases R2 pct",
					"filename","dragen_map_metrics.txt"
					)));
			
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","Total input reads,Total input bases,Total input bases R1,Total input bases R2,Average input read length,Total trimmed reads,Total trimmed bases,Average bases trimmed per read,Average bases trimmed per trimmed read,Remaining poly-G K-mers R1 3prime,Remaining poly-G K-mers R2 3prime,Poly-G soft trimmed reads unfiltered R1 3prime,Poly-G soft trimmed reads unfiltered R2 3prime,Poly-G soft trimmed reads filtered R1 3prime,Poly-G soft trimmed reads filtered R2 3prime,Poly-G soft trimmed bases unfiltered R1 3prime,Poly-G soft trimmed bases unfiltered R2 3prime,Poly-G soft trimmed bases filtered R1 3prime,Poly-G soft trimmed bases filtered R2 3prime,Total filtered reads,Reads filtered for minimum read length R1,Reads filtered for minimum read length R2",
					"filename","dragen-trimmer-metrics-table.txt"
					)));
			
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","Aligned bases,Aligned bases in genome,Aligned bases in genome pct,Average alignment coverage over genome,Uniformity of coverage (PCT > 0.2*mean) over genome,Uniformity of coverage (PCT > 0.4*mean) over genome,PCT of genome with coverage [1500x:inf),PCT of genome with coverage [1000x:inf),PCT of genome with coverage [500x:inf),PCT of genome with coverage [100x:inf),PCT of genome with coverage [50x:inf),PCT of genome with coverage [20x:inf),PCT of genome with coverage [15x:inf),PCT of genome with coverage [10x:inf),PCT of genome with coverage [3x:inf),PCT of genome with coverage [1x:inf),PCT of genome with coverage [0x:inf),PCT of genome with coverage [1000x:1500x),PCT of genome with coverage [500x:1000x),PCT of genome with coverage [100x:500x),PCT of genome with coverage [50x:100x),PCT of genome with coverage [20x:50x),PCT of genome with coverage [15x:20x),PCT of genome with coverage [10x:15x),PCT of genome with coverage [3x:10x),PCT of genome with coverage [1x:3x),PCT of genome with coverage [0x:1x),Average chr X coverage over genome,Average chr Y coverage over genome,Average mitochondrial coverage over genome,Average autosomal coverage over genome,Median autosomal coverage over genome,Mean/Median autosomal coverage ratio over genome,Aligned reads,Aligned reads in genome,Aligned reads in genome pct",
					"filename","dragen_wgs_cov_metrics.txt"
					)));
			
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","Variants,Biallelic,Multiallelic,SNP,Indel,Ins,Del,Hom ins,Het ins,Hom del,Het del,Het indel,X/Y SNP ratio,Ti/Tv,Het,Hom,Het/Hom,Callability,Autosome callability",
					"filename","dragen-gvcf-metrics-table.txt"
					)));
			
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","Total,Total pct,Biallelic,Biallelic pct,Multiallelic,Multiallelic pct,SNPs,SNPs pct,Insertions (Hom),Insertions (Hom) pct,Insertions (Het),Insertions (Het) pct,Deletions (Hom),Deletions (Hom) pct,Deletions (Het),Deletions (Het) pct,Indels (Het),Indels (Het) pct,Chr X number of SNPs over genome,Chr Y number of SNPs over genome,(Chr X SNPs)/(chr Y SNPs) ratio over genome,SNP Transitions,SNP Transversions,Ti/Tv ratio,Heterozygous,Homozygous,Het/Hom ratio,In dbSNP,In dbSNP pct,Not in dbSNP,Not in dbSNP pct,Percent Callability,Percent Autosome Callability,Insertions,Deletions,Indels,Insertions pct,Deletions pct,Indels pct",
					"filename","dragen_gvcf_metrics.txt"
					)));
			
			
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","Unique,Duplicated,Unmapped",
					"filename","mapping_dup_percentage_plot_Unique_vs_duplicated_vs_unmapped.txt"
					)));
			
			
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "Aln reads","Depth","Med aut cov",String.valueOf(GE_UNICODE)+"20x",String.valueOf(GE_UNICODE)+"50x"),
					"filename","dragen-cov-metrics-own-section-wgs-table.txt"
					)));
			
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "Duplicates","Singletons","Mate mapped to diff chr"),
					"filename","samtools-flagstat-dp_Percentage_of_total.txt"
					)));
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "MEDIAN_INSERT_SIZE"),
					"filename","multiqc_picard_insertSize.txt"
					)));
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "total_deduplicated_percentage","%GC","Sequence length"),
					"filename","multiqc_fastqc.txt"
					)));
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "raw_total_sequences","filtered_sequences","reads_mapped_and_paired","reads_unmapped"),
					"filename","glob:multiqc_samtools_stats*.txt"
					)));
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "total_passed","total_failed","reads_mapped_and_paired","reads_unmapped"),
					"filename","glob:multiqc_samtools_flagstat*.txt"
					)));
			
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "Mean Read Length"),
					"filename","glob:mqc_picard_alignment_readlength_plot_*.txt"
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "percent_trimmed"),
					"filename","multiqc_cutadapt.txt"
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data",String.join(",", "Intergenic","TTS","Unassigned","exon","intron","promoter-TSS"),
					"filename","multiqc_mlib_peak_annotation-plot.txt"
					)));
			
			handlers.add(new JsonHandler(Maps.of(
					"data",String.join(",",
							"Unmapped reads pct",
							"Median autosomal coverage over genome",
			                "Average chr X coverage over genome",
			                "Average chr Y coverage over genome",
			                "Singleton reads (itself mapped; mate unmapped) pct",
							"Biallelic pct",
							"Multiallelic pct",
			                "Deletions (Het) pct",
			                "Deletions (Hom) pct",
			                "Deletions pct",
			                "Hard-clipped bases pct",
			                "In dbSNP pct",
			                "Indels pct",
			                "Soft-clipped bases pct",
			                "Soft-clipped bases R1 pct",
			                "Soft-clipped bases R2 pct",
			                "Hard-clipped bases R1 pct",
			                "Hard-clipped bases R2 pct",
							"Aligned bases in genome pct",
			                "PCT of genome with coverage [  50x: inf)",
			                "PCT of genome with coverage [  20x: inf)",
			                "PCT of genome with coverage [  15x: inf)",
			                "PCT of genome with coverage [  10x: inf)",
							"10_x_pc","30_x_pc","50_x_pc",
			                "Number of duplicate marked reads pct",
			                "Mismatched bases R1 pct",
			                "Mismatched bases R2 pct",
			                "Paired reads mapped to different chromosomes pct",
			                "Properly paired reads pct",
							"mean_coverage","median_coverage",
							"reads_mapped_percent","reads_duplicated_percent","insert_size_average","error_rate","average_length","average_quality"),
					"filename","multiqc_data.json"
					)));
			/*
			handlers.add(new DistributionHandler(Maps.of(
					"data",String.join(",","10,20,30"),
					"section_name","Cumulative Coverage. % of bases mapped more than ${title}",
					"filename","mosdepth-cumcoverage-dist-id.txt"
					)));
			*/
			
			if(customMapping!=null) {
				try(BufferedReader br = IOUtils.openPathForBufferedReading(customMapping)) {
					final List<Map.Entry<String, String>> rows= new ArrayList<>();
					for(;;) {
						String line = br.readLine();
						if(StringUtils.isBlank(line)) {
							final Map<String,String> hash = new HashMap<>();
							for(String key : rows.stream().map(KV->KV.getKey()).collect(Collectors.toSet())) {
								hash.put(key, rows.stream().filter(KV->KV.getKey().equals(key)).map(KV->KV.getValue()).collect(Collectors.joining(",")));
								}
							handlers.add(new BoxPlotHandler(hash));
							if(line==null) break;
							rows.clear();
							continue;
							}
						rows.add(new AbstractMap.SimpleEntry<>(
							StringUtils.substringBefore(line,":").trim(),
							StringUtils.substringAfter(line,":").trim()
							));
					}
				}
			}
			
			
			final List<Path> filesForInput =  IOUtils.unrollPaths(args).
					stream().flatMap(F->{
				try {
					if(Files.isDirectory(F)) {
						IOUtils.assertDirectoryIsReadable(F);
						return Files.list(F);
						}
					return Stream.of(F);
					}
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
				}).collect(Collectors.toList());
			
			
			
			LOG.info(filesForInput);
			IOUtil.assertDirectoryIsWritable(outputDirectory);
			
			if(sample2collectionPath!=null) {
				this.sampleCollection.load(this.sample2collectionPath);
				}
			else {
				final Path dragenPath = filesForInput.stream().filter(F->F.endsWith("dragen_ploidy.txt") || F.endsWith("dragen_ploidy_table.txt")).findFirst().orElse(null);
				if(dragenPath!=null) {
					final FileContent fc=readFileContent(dragenPath);
					final String dragenPloidy =fc.fileHeader.containsKey("Ploidy estimation")?"Ploidy estimation":"Sex";

					
					if(fc.fileHeader.containsKey(dragenPloidy) && fc.fileHeader.containsKey("Sample")) {
						LOG.warning("no sample2collection defined. Using ploidy from dragen: "+dragenPloidy+" in "+dragenPath);
						for(FileHeader.RowMap row:fc.rows) {
							if(row.get("Sample").isEmpty()) continue;
							if(row.get(dragenPloidy).isEmpty()) continue;
							this.sampleCollection.insert(row.get("Sample"),row.get(dragenPloidy));
							}
						}
					}
				}
			
			if(this.sampleCollection.isEmpty()) {
				LOG.error("no sample/collection defined");
				return -1;
				}
			
			for(final Handler handler : handlers) {
				final Set<Path> pathSet = filesForInput.stream().
					filter(PATH->handler.isHandlerForFile(PATH.getFileName().toString())).
					collect(Collectors.toSet());
				if(pathSet.isEmpty()) {
					continue;
					}
				LOG.info("running "+pathSet);
				handler.apply(pathSet);
				}
			return 0;
		} catch (Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
public static void main(String[] args) {
	new MultiqcPostProcessor().instanceMainWithExit(args);
	}
}
