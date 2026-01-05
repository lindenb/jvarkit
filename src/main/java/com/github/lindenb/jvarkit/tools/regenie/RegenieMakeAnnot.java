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
package com.github.lindenb.jvarkit.tools.regenie;
import java.io.BufferedReader;
import java.io.Closeable;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Random;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFReader;
/**
## mask

example mask for snpeff

```
3_prime_UTR_truncation	 0.5	UTR,UTR3
3_prime_UTR_variant	 0.1	UTR,UTR3
5_prime_UTR_premature_start_codon_gain_variant	 0.2	UTR,UTR5
5_prime_UTR_truncation	 0.5	UTR,UTR5
5_prime_UTR_variant	 0.2	UTR,UTR5
bidirectional_gene_fusion	 1.0	.
conservative_inframe_deletion	 0.1	protein_altering
conservative_inframe_insertion	 0.3	protein_altering
disruptive_inframe_deletion	 0.2	protein_altering
disruptive_inframe_insertion	 0.2	protein_altering
downstream_gene_variant	 0.1	downstream,updownstream
exon_loss_variant	 1.0	protein_altering
exon_loss_variant	 1.0	protein_altering
frameshift_variant	 0.4	.
gene_fusion	 0.9	protein_altering
initiator_codon_variant	0.5	protein_altering
intergenic_region	 0.001	.
intragenic_variant	 0.01	.
intron_variant	 0.05	intronic
missense_variant	 0.9	protein_altering
non_coding_transcript_exon_variant	 0.1	non_coding
non_coding_transcript_variant	 0.1	non_coding
splice_acceptor_variant	 0.5	protein_altering,splice
splice_donor_variant	 0.5	protein_altering,splice
splice_region_variant	 0.5	protein_altering,splice
start_lost	 0.6	protein_altering
start_retained_variant	0.1	synonymous
stop_gained	 0.9	protein_altering
stop_lost	 0.6	protein_altering
stop_retained_variant	 0.2	synonymous
synonymous_variant	 0.1	synonymous
upstream_gene_variant	 0.1	upstream,updownstream
```


## Example

```
bcftools view -O v 'chr4.bcf' |\
	java -Djava.io.tmpdir=TMP -jar "jvarkit.jar" regenieslidingannot \
		--window-size "3000" \
		--window-shift "1500" \
		-f 0.01,0.001 |\
	java -Djava.io.tmpdir=TMP -jar "jvarkit.jar" regeniemakeannot \
		--prefix "chr4_3000_1500_chunk" \
		-o ${PWD}/OUT \
		--gzip \
		-N 5000

$ find OUT/ | sort | cut -d '/' -f 14- | tail
OUT/chr4_3000_1500_chunk097_000.setfile.txt.gz
OUT/chr4_3000_1500_chunk098_000.aaf.txt.gz
OUT/chr4_3000_1500_chunk098_000.annot.txt.gz
OUT/chr4_3000_1500_chunk098_000.mask.txt.gz
OUT/chr4_3000_1500_chunk098_000.setfile.txt.gz
OUT/chr4_3000_1500_chunk099_000.aaf.txt.gz
OUT/chr4_3000_1500_chunk099_000.annot.txt.gz
OUT/chr4_3000_1500_chunk099_000.mask.txt.gz
OUT/chr4_3000_1500_chunk099_000.setfile.txt.gz
OUT/manifest.tsv

```

 */
@Program(name="regeniemakeannot",
description="Create annotation files for regenie from a TSV input file",
keywords={"vcf","regenie","burden"},
creationDate="20250311",
modificationDate="20250901",
generate_doc = true,
jvarkit_amalgamion = true
)
public class RegenieMakeAnnot extends Launcher {
	private static final Logger LOG = Logger.of(RegenieMakeAnnot.class);

	@Parameter(names={"-o","--output"},description="output dir.",required = true)
	private Path outputDir=null;
	@Parameter(names={"--prefix"},description="prefix for output files")
	private String filePrefix = "chunk";
	@Parameter(names={"-m","--masks"},description="mask file. TSV file. no header. 3 columns prediction_name/score/comma-separated-mask_names. if undefined, will produce one mask per prediction")
	private Path masksFile = null;
	@Parameter(names={"-N"},description="max item per chunck. -1: no limit ")
	private int max_items=10_000;
	@Parameter(names={"--vcf"},description="keep variant id only if it's ID was found in this VCF file")
	private Path externalVCFPath = null;
	@Parameter(names={"--gzip","-Z"},description="compress output files")
	private boolean gzip_files=false;
	@Parameter(names={"--reserve"},description="reserve 'n' output files of non-overlaping gene/target")
	private int reserve_output = 20;
	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

	private final Map<String,Prediction> predictions_hash = new HashMap<>();//tmpFile = Files.createTempFile("regenie", ".tmp");

	
	 
	
	private class Output implements Closeable{
		private final int index;
		private int count_written=0;
		private int batch_size=0;
		PrintWriter annot = null;
		PrintWriter setfile = null;
		PrintWriter mask = null;
		PrintWriter aaf = null;
		final Set<String> seen_predictions = new HashSet<>();
		final IntervalTreeMap<Target> treeMap = new IntervalTreeMap<>();
		Output(int index) {
			this.index=index;
			}
		void prepare(final PrintWriter manifest) throws IOException {
			if(RegenieMakeAnnot.this.max_items>0 && count_written>=RegenieMakeAnnot.this.max_items ) {
				close();
				count_written=0;
				batch_size++;
				}
			if(count_written==0) {
				final String prefix = String.format("%s%03d_%03d.",filePrefix,index,batch_size);
				final String suffix = ".txt"+(gzip_files?".gz":"");
				final Path p1 = outputDir.resolve(prefix+"annot"+suffix);
				final Path p2 = outputDir.resolve(prefix+"setfile"+suffix);
				final Path p3 = outputDir.resolve(prefix+"mask"+suffix);
				final Path p4 = outputDir.resolve(prefix+"aaf"+suffix);
				this.annot = IOUtils.openPathForPrintWriter(p1);
				this.setfile = IOUtils.openPathForPrintWriter(p2);
				this.mask = IOUtils.openPathForPrintWriter(p3);
				this.aaf = IOUtils.openPathForPrintWriter(p4);
				manifest.print(p1.toString());
				manifest.print("\t");
				manifest.print(p2.toString());
				manifest.print("\t");
				manifest.print(p3.toString());
				manifest.print("\t");
				manifest.print(p4.toString());
				manifest.println();
				}
			count_written++;
			}
		@Override
		public void close() {
			if(count_written>0) 
				{
				final Set<String> distinct_masks = RegenieMakeAnnot.this.predictions_hash.values().
						stream().
						flatMap(SET->SET.masks.stream()).
						collect(Collectors.toSet());
				
				for (final String mask_name : distinct_masks) {
					final Set<String> predset = this.seen_predictions.
							stream().
							filter(P->RegenieMakeAnnot.this.predictions_hash.get(P).masks.contains(mask_name)).
							collect(Collectors.toSet())
							;
					if(predset.isEmpty()) continue;
					this.mask.print(mask_name);
					this.mask.print("\t");
					this.mask.print(String.join(",",predset));
					this.mask.println();
					}
				
				this.annot.flush();
				this.setfile.flush();
				this.mask.flush();
				this.aaf.flush();
				this.annot.close();
				this.setfile.close();
				this.mask.close();
				this.aaf.close();
				this.seen_predictions.clear();
				}
			}
		}
	
	/** min max of a gene */
	private static class Target {
		String contig;
		int start;
		int end;
		Output output=null;
		}

	
	private static class Variation implements Locatable {
		String contig;
		int pos;
		String id;
		String gene;
		String prediction;
		double score;
		double cadd;
		double frequency;
		byte is_singleton;
		@Override public String getContig() { return contig; }
		@Override public int getStart() { return pos; }
		@Override public int getEnd() { return pos; }
		@Override
		public String toString() {
			return contig+":"+pos+":"+id+":"+gene;
			}
		}

	
	private static class Prediction {
		final String name;
		OptionalDouble score=OptionalDouble.empty();
		final Set<String> masks = new HashSet<>();

		Prediction(final String name) {
			this.name  = name;
			}
		}


	private static class RowCodec extends AbstractDataCodec<Variation> {
		@Override
		public void encode(DataOutputStream dos, final Variation variant) throws IOException {
			dos.writeUTF(variant.contig);
			dos.writeInt(variant.pos);
			dos.writeUTF(variant.id);
			dos.writeUTF(variant.gene);
			dos.writeUTF(variant.prediction);
			dos.writeDouble(variant.score);
			dos.writeDouble(variant.cadd);
			dos.writeDouble(variant.frequency);
			dos.writeByte(variant.is_singleton);
		}

		@Override
		public Variation decode(DataInputStream dis) throws IOException {
			final Variation v = new Variation();
			try {
				v.contig = dis.readUTF();
			} catch (EOFException err) {
				return null;
			}
			v.pos = dis.readInt();
			v.id = dis.readUTF();
			v.gene = dis.readUTF();
			v.prediction = dis.readUTF();
			v.score = dis.readDouble();
			v.cadd = dis.readDouble();
			v.frequency = dis.readDouble();
			v.is_singleton = dis.readByte();
			return v;
			}

		@Override
		public RowCodec clone() {
			return new RowCodec();
		}
	}
	
	protected String fixContig(final String ctg) {
		//if(ctg.equals("X") || ctg.equals("chrX")) return "23";
		//if(ctg.equals("Y") || ctg.equals("chrY")) return "24";
		if(ctg.startsWith("chr")) return ctg.substring(3);
		return ctg;
		}

	
	private boolean isEmpty(final String s) {
		return StringUtils.isBlank(s) || s.equals(".");
	}

	private static int compareVariant(Variation t1, Variation t2, int level) {
		int i = t1.contig.compareTo(t2.contig);// CHROMOSOME
		if (i != 0)
			return i;
		i = t1.gene.compareTo(t2.gene);// GENE
		if (i != 0)
			return i;
		if (level == 1)
			return 0;
		return Integer.compare(t1.pos, t2.pos);
	}

	
	
	
	@Override
	public int doWork(final List<String> args) {
		//Path tmpFile=null;
		BufferedVCFReader bufferedVCFReader = null;
		VCFReader vcfReader = null;
		ContigNameConverter vcfCtgConverter=null;
		final List<Output>  outputs = new ArrayList<>();
		try {
			final String input = oneFileOrNull(args);
			final Pattern spaces_regex = Pattern.compile("[ \t]");
			final Map<String, Target> target_hash = new HashMap<>(50_000);
			
			if(this.externalVCFPath!=null) {
				vcfReader  = VCFReaderFactory.makeDefault().open(this.externalVCFPath,true);
				final SAMSequenceDictionary dict = new SequenceDictionaryExtractor().extractRequiredDictionary(vcfReader.getHeader());
				vcfCtgConverter = ContigNameConverter.fromOneDictionary(dict);
				bufferedVCFReader = new BufferedVCFReader(vcfReader, 10000);
				bufferedVCFReader.setSimplifier(VC->new VariantContextBuilder(VC).attributes(Collections.emptyMap()).noGenotypes().make());
				}
 			
			if(this.masksFile!=null) {
				try(BufferedReader br= IOUtils.openPathForBufferedReading(masksFile)) {
					for(;;) {
						String line = br.readLine();
						if(line==null) break;
						if(line.startsWith("#")) continue;
						final String[] tokens = spaces_regex.split(line);
						if(tokens.length!=3) throw new JvarkitException.TokenErrors(3, tokens);
						if(predictions_hash.containsKey(tokens[0])) {
							throw new IllegalArgumentException("duplicate prediction "+tokens[0]+" in "+masksFile);
							}
						final Prediction pred = new Prediction(tokens[0]);
						if(!isEmpty(tokens[1])) {
							pred.score = OptionalDouble.of(Double.parseDouble(tokens[1]));
							}
						if(pred.name.equals(RegenieFunctionalAnnot.FIRST_INTRON)) {
							pred.masks.add(RegenieFunctionalAnnot.FIRST_INTRON);
							}
						else
							{
							pred.masks.addAll(Arrays.asList(CharSplitter.COMMA.split(tokens[2])));
							}
						

						pred.masks.remove("");
						pred.masks.remove(".");
						if(pred.masks.isEmpty()) {
							throw new IllegalArgumentException("no mask in "+line+" in "+masksFile);
							}
						predictions_hash.put(pred.name, pred);
					}
				}
			}
			
			final SortingCollection<Variation> sorter = SortingCollection.newInstance(Variation.class, new RowCodec(),
					(A, B) -> compareVariant(A, B, 2), writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths());
			sorter.setDestructiveIteration(true);

			try(BufferedReader br= (input==null?IOUtils.openStdinForBufferedReader():IOUtils.openPathForBufferedReading(Paths.get(input)))) {
				String line = br.readLine();
				if(line==null) throw new IOException("cannot read header line");
				if(line.startsWith("#")) line=line.substring(1);
				final FileHeader fileHeader = new FileHeader(line, S->Arrays.asList(spaces_regex.split(S)));
				int col_chrom = fileHeader.containsKey("CHROM")?fileHeader.getColumnIndex("CHROM"):fileHeader.getColumnIndex("CONTIG");
				int col_id = fileHeader.getColumnIndex("ID");
				int col_pos = fileHeader.getColumnIndex("POS");
				int col_gene = fileHeader.getColumnIndex("GENE");
				int col_prediction = fileHeader.getColumnIndex("ANNOTATION");
				int col_score = fileHeader.containsKey("SCORE")?fileHeader.getColumnIndex("SCORE"):-1;
				int col_cadd = fileHeader.containsKey("CADD")?fileHeader.getColumnIndex("CADD"):-1;
				int col_singleton = fileHeader.getColumnIndex("SINGLETON");
				int col_freq = fileHeader.getColumnIndex("FREQ");

					
				while((line=br.readLine())!=null) {
					final FileHeader.RowMap row = fileHeader.toMap(line);
					final Variation v = new Variation();
					v.id = row.at(col_id);
					v.contig = fixContig(row.at(col_chrom));
					v.gene =row.at(col_gene);
					v.prediction = row.at(col_prediction);
					if(this.masksFile==null) {
						if(!this.predictions_hash.containsKey(v.prediction)) {
							final Prediction p  = new Prediction(v.prediction);
							p.masks.add(v.prediction);
							p.score=OptionalDouble.of(1.0);
							this.predictions_hash.put(p.name,p);
							}
						}
					else
						{
						if(!this.predictions_hash.containsKey(v.prediction)) {
							throw new IllegalArgumentException("got prediction "+v.prediction+" in "+line+" but no mask was defined : " + this.predictions_hash.keySet());
							}
						}
					final String scoreStr=col_score<0?"":row.at(col_score);
					if(isEmpty(scoreStr)) {
						final Prediction p  = this.predictions_hash.get(v.prediction);
						if(p==null) throw new IllegalArgumentException("got prediction "+v.prediction+" in "+line+" but no mask was defined : " + this.predictions_hash.keySet());
						if(!p.score.isPresent()) throw new IllegalArgumentException("no score was defined form prediction "+v.prediction);
						v.score = p.score.getAsDouble();
						}
					else
						{
						v.score = Double.parseDouble(scoreStr);
						}
					v.pos = Integer.parseInt(row.at(col_pos));
					
					final String caddStr=col_cadd<0?"":row.at(col_cadd);
					if(isEmpty(caddStr)) {
						v.cadd = 0.0;
						}
					else
						{
						v.cadd = Double.parseDouble(caddStr);
						}
					
					v.is_singleton = Byte.parseByte(row.at(col_singleton));
					v.frequency = Double.parseDouble(row.at(col_freq));

					
					/* check this ID is present in the VCF file */
					if(bufferedVCFReader!=null) {
						final String ctg=vcfCtgConverter.apply(v.contig);
						if(StringUtils.isBlank(ctg)) continue;
						if(bufferedVCFReader.query(ctg, v.pos, v.pos).stream().
							noneMatch(VC->VC.hasID() && VC.getID().equals(v.id))) continue;
						}
					
					
					final String gene_id = v.contig+"~"+v.gene;
					Target t = target_hash.get(gene_id);
					if(t==null) {
						t =  new Target();
						t.contig = v.contig;
						t.start = v.pos;
						t.end = v.pos;
						target_hash.put(gene_id, t);
						}
					t.start = Math.min(t.start, v.pos);
					t.end = Math.max(t.end, v.pos);
					sorter.add(v);
					}
				}
			sorter.doneAdding();
			
			final Random choose_any_output =  new Random(0L);
			while(outputs.size()<Math.max(2, reserve_output)) {
				outputs.add(new Output(outputs.size()));
				}
			for(Target t:target_hash.values()) {
				final Interval loc = new Interval(t.contig, t.start, t.end);
				
				final List<Output> available_output = outputs.stream().
						filter(O->!O.treeMap.containsOverlapping(loc)).
						collect(Collectors.toList());
				Output output = null;
				if(available_output.isEmpty()) {
					output = new Output(outputs.size());
					outputs.add(output);
					}
				else
					{
					output = available_output.get(choose_any_output.nextInt(available_output.size()));
					}
				if(output.treeMap.containsOverlapping(loc) ) {
					throw new IllegalStateException("found overlap ?");
					}
				output.treeMap.put(loc, t);
				t.output = output;
				}

			try(PrintWriter manifestW = IOUtils.openPathForPrintWriter(this.outputDir.resolve("manifest.tsv"))) {
			try (CloseableIterator<Variation> iter0 = sorter.iterator()) {
				try (EqualRangeIterator<Variation> iter = new EqualRangeIterator<>(iter0, (A, B) -> compareVariant(A, B, 1))) {
					while (iter.hasNext()) {
						final List<Variation> gene_variants = iter.next();
						// sort and remove duplicates, if any
						Collections.sort(gene_variants,(A,B)->{
							int i = A.contig.compareTo(B.contig);
							if(i!=0) return i;
							i = Integer.compare(A.pos, B.pos);
							if(i!=0) return i;
							i = A.id.compareTo(B.id);
							return i;
							});
						int x=0;
						while(x+1< gene_variants.size()) {
							final Variation v1 = gene_variants.get(x  );
							final Variation v2 = gene_variants.get(x+1);
							if(v1.id.equals(v2.id)) {
								gene_variants.remove(x+1);
							} else {
								x++;
								}
							}
						final Variation first = gene_variants.get(0);
						final String gene_id = first.contig+"~"+first.gene;
						final Target target = target_hash.get(gene_id);
						if(target==null) throw new IllegalStateException(gene_id);
						final Output output =  target.output;
						output.prepare(manifestW);
						
						
						for (Variation v : gene_variants) {
							output.annot.print(v.id);
							output.annot.print(" ");
							output.annot.print(v.gene);
							output.annot.print(" ");
							output.annot.print(v.prediction);
							output.annot.print(" ");
							output.annot.print(v.score);
							output.annot.print(" ");
							output.annot.print(v.cadd);
							output.annot.println();
							
							output.seen_predictions.add(v.prediction);
							
							
							output.aaf.print(v.id);
							output.aaf.print(" ");
							output.aaf.print(v.frequency);
							output.aaf.print(" ");
							output.aaf.print(v.is_singleton);
							output.aaf.println();
							}
						output.setfile.print(first.gene);// gene
						output.setfile.print("\t");
						output.setfile.print(first.contig);// contig
						output.setfile.print("\t");
						output.setfile.print((int) gene_variants.stream().mapToInt(it -> it.pos).average().getAsDouble());
						output.setfile.print("\t");
						output.setfile.print(gene_variants.stream().map(it -> it.id).collect(Collectors.joining(",")));
						output.setfile.println();
						}
					}
				}
			manifestW.flush();
			}
			sorter.cleanup();
			
			if(vcfReader!=null) {
				vcfReader.close();
				vcfReader=null;
				}
			if(bufferedVCFReader!=null) {
				bufferedVCFReader.close();
				bufferedVCFReader=null;
				}		
			
			for(Output o: outputs) {
				o.close();
				}
			return 0;
		} catch (Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
			if(vcfReader!=null) {
				try { vcfReader.close(); } catch(IOException err) {}
				}
			if(bufferedVCFReader!=null) {
				try { bufferedVCFReader.close(); } catch(IOException err) {}
				}			
			//if(tmpFile!=null) Files.deleteIfExists(tmpFile);
			for(Output o: outputs) {
				o.close();
			}
		}
	}

	public static void main(final String[] args) {
		new RegenieMakeAnnot().instanceMainWithExit(args);
	}
}
