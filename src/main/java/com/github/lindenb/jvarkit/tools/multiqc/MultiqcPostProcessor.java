/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.util.IOUtil;

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
java -jar jeter.jar multiqc_data -o OUTDIR2
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
modificationDate="20240708",
jvarkit_amalgamion = true,
menu="Utilities"
)
public class MultiqcPostProcessor extends Launcher {
	private static final Logger LOG = Logger.build(MultiqcPostProcessor.class).make();
	private final char GE_UNICODE ='\u2265';
	@Parameter(names={"-o","--output"},description="output directory",required = true)
	private Path outputDirectory = null;
	@Parameter(names={"--sample2collection"},description=SamplePopulation.OPT_DESC)
	private Path sample2collectionPath=null;
	@Parameter(names={"--custom"},description="custom mapping file (undocumented")
	private Path customMapping=null;

	private final SamplePopulation sampleCollection = new SamplePopulation();
	
	private static class FileContent {
		FileHeader fileHeader;
		List<FileHeader.RowMap> rows;
	}
	
	private abstract class Handler {
		public abstract boolean isHandlerForFile(final String filename);
		public abstract void apply(Path f);
		};
	
		private class BoxPlotHandler extends Handler {
			final Map<String,String> properties=new HashMap<>();
			BoxPlotHandler(Map<String,String> prop) {
				this.properties.putAll(prop);
				}
			@Override
			public boolean isHandlerForFile(String fname) {
				String[] s=  CharSplitter.COMMA.split(this.properties.getOrDefault("filename", ""));
				return Arrays.stream(s).anyMatch(S->S.equals(fname));
				}
			protected String getSampleColumn() {
				return this.properties.getOrDefault("sample", "Sample");
				}
			public String[] getDataColumns() {
				String s = this.properties.getOrDefault("data", "");
				return CharSplitter.COMMA.split(s);
				}
			protected String fixSampleName(final String sn) {
				return sn;
				}
			/** for some tools multiqc normalize by x. Eg: samtools stats results is 'x' in MB (0.0001 ) */
			protected double getFactor() {
				return Double.parseDouble(this.properties.getOrDefault("factor", "1.0"));
				}
			private boolean isDouble(String s) {
				if(StringUtils.isBlank(s)) return false;
				try {
					Double.parseDouble(s);
					return true;
				} catch(NumberFormatException f) {
					return false;
				}
			}
			@Override
			public void apply(final Path f) {
				try {
					final FileContent fc = readFileContent(f);
					fc.fileHeader.assertColumnExists(getSampleColumn());
					final String[] dataColumns = getDataColumns();
					final double factor = getFactor();
					for(final String dataCol: dataColumns) {
						if(StringUtils.isBlank(dataCol)) continue;
						if(!fc.fileHeader.containsKey(dataCol)) {
							LOG.warn("no column "+dataCol+" in "+f+" "+fc.fileHeader);
							continue;
							}
						//remove unicode
						final String datColPlain = dataCol.replaceAll(String.valueOf(GE_UNICODE),">=");
						// make a copy of SamplePopulation
						final SamplePopulation snpop =new SamplePopulation(MultiqcPostProcessor.this.sampleCollection);
						snpop.retainSamples(fc.rows.stream().
								map(ROW-> fixSampleName(ROW.get(getSampleColumn()))).
								collect(Collectors.toSet())
								);
						
						fc.fileHeader.assertColumnExists(dataCol);
						/* is it worth plotting ?*/
						boolean found_deviation=false;
						final Map<String,double[]> pop2data=new TreeMap<>();
						for(final SamplePopulation.Population pop:snpop.getPopulations()) {
							final double[] values= fc.rows.stream().
								filter(ROW->pop.containsKey(fixSampleName(ROW.get(getSampleColumn())))).
								map(ROW->ROW.get(dataCol)).
								filter(V->isDouble(V)).
								mapToDouble(V->Double.parseDouble(V)*factor).
								sorted().
								toArray();
							/* is it worth plotting ?*/
							if(values.length>1 && values[0] < values[values.length-1]) {
								found_deviation=true;
								}
							pop2data.put(pop.getName(), values);
							}
						if(pop2data.isEmpty() || !found_deviation) {
							LOG.info("no data to plot");
							return;
						}
						final Path outfile = outputDirectory.resolve(IOUtils.getFilenameWithoutCommonSuffixes(f)+"_"+datColPlain.replaceAll("[^A-Za-z0-9_]+", "_").toLowerCase()+"_mqc.json");
						LOG.info("writing "+outfile);
						try(PrintWriter pw = IOUtils.openPathForPrintWriter(outfile)) {
							final String id = StringUtils.md5( f.getFileName().toString()+dataCol);
							try(JsonWriter w = new JsonWriter(pw)) {
								w.beginObject();
								
								w.name("parent_id");
								w.value("p"+StringUtils.md5(MultiqcPostProcessor.class.getName()));
								w.name("parent_name");
								w.value("Per Population");
								w.name("parent_description");
								w.value("MULTIQC data grouped by population using jvarkit "+MultiqcPostProcessor.class.getSimpleName());
								
								w.name("id");
								w.value("section_id"+id);
								w.name("section_name");
								w.value(datColPlain+ " per population");
								w.name("description");
								w.value(datColPlain+ " per population from "+f.getFileName());
								
								
								w.name("plot_type");
								w.value("box");
								w.name("pconfig");
								w.beginObject();
								w.name("id");
								w.value("plot__id"+id);
								w.name("title");
								w.value(datColPlain+ " per population");
								w.name("xlab");
								w.value(datColPlain);
								w.name("ylab");
								w.value("collection");
								w.endObject();
								
								w.name("data");
								w.beginObject();
								for(final String popName:pop2data.keySet()) {
									final SamplePopulation.Population pop = snpop.getPopulationByName(popName);
									w.name(pop.getName()+" (N="+pop.size()+")");
									w.beginArray();
									for(double v: pop2data.get(popName)) {
										w.value(v);
										}
									w.endArray();
									}
								w.endObject();
								
								w.endObject();
								w.flush();
								}
							pw.flush();
							}
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
					"filename","bcftools-stats-singletons.txt"
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","depth",
					"filename","bcftools-stats-sequencing-depth.txt"
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","nSNPs,nIndels",
					"filename","bcftools-stats-sites.txt"
					)));
			handlers.add(new BoxPlotHandler(Maps.of(
					"data","tstv",
					"filename","bcftools-stats-tstv.txt"
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
					"data",String.join(",", "Aln reads","Depth","Med aut cov",String.valueOf(GE_UNICODE)+"20x",String.valueOf(GE_UNICODE)+"50x"),
					"filename","dragen-cov-metrics-own-section-wgs-table.txt"
					)));
			
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
			
			
			final String input = super.oneAndOnlyOneFile(args);
			final Path inputDirectory = Paths.get(input);
			IOUtils.assertDirectoryIsReadable(inputDirectory);
			
			

			final List<Path> fileInDir = Files.list(inputDirectory).
					collect(Collectors.toList());
			LOG.info(fileInDir);
			IOUtil.assertDirectoryIsWritable(outputDirectory);
			
			if(sample2collectionPath!=null) {
				this.sampleCollection.load(this.sample2collectionPath);
				}
			else {
				final Path dragenPath = fileInDir.stream().filter(F->F.endsWith("dragen_ploidy.txt") || F.endsWith("dragen_ploidy_table.txt")).findFirst().orElse(null);
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
			
			for(Path f : fileInDir ) {
				LOG.info("test "+f);
				final Handler handler = handlers.stream().
					filter(H->H.isHandlerForFile(f.getFileName().toString())).
					findFirst().
					orElse(null);
				if(handler==null) continue;
				LOG.info("running "+f);
				handler.apply(f);
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
