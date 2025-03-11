package com.github.lindenb.jvarkit.tools.regenie;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public abstract class AbstractRegenieAnnot extends Launcher {
	private static final String CADD_PHRED = "CADD_PHRED";
	private static final String GNOMAD_AF = "gnomad_genome_AF_NFE";
	
	@Parameter(names = "-f", description = "comma separated of Allele frequencies , I will use the highest to discard frequent variants.")
	private String freqStr="0.01";
	@Parameter(names = "-a", description = "output annotation file", required = true)
	private Path annotationFileOut = null;
	@Parameter(names = "-s", description = "set list file output", required = true)
	private Path setListFileOut = null;
	@Parameter(names = "-m", description = "mask file output", required = true)
	private Path maskFileOut = null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

	
	
	protected static class Variation {
		String contig;
		int pos;
		String id;
		String gene;
		String prediction;
		double score;
		double cadd;
		}

	
	protected static class Prediction {
		final String name;
		final double score;
		boolean found = false;
		final Set<String> masks = new HashSet<>();

		Prediction(String name,double score) {
			this.name  = name;
			this.score  = score;
			}
		}


	private static class RowCodec extends AbstractDataCodec<Variation> {
		@Override
		public void encode(DataOutputStream dos, Variation variant) throws IOException {
			dos.writeUTF(variant.contig);
			dos.writeInt(variant.pos);
			dos.writeUTF(variant.id);
			dos.writeUTF(variant.gene);
			dos.writeUTF(variant.prediction);
			dos.writeDouble(variant.score);
			dos.writeDouble(variant.cadd);
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
			return v;
			}

		@Override
		public RowCodec clone() {
			return new RowCodec();
		}
	}

	protected final Map<String,Prediction> name2prediction = new HashMap<>();
	
	protected Prediction makeScore(String pred, double score,String...masks) {
		final Prediction p = new Prediction(pred,score);
		p.masks.addAll(Arrays.asList(masks));
		if(name2prediction.containsKey(p.name)) throw new IllegalArgumentException("duplicate "+p.name);
		name2prediction.put(pred, p);
		return p;
		}

	protected AbstractRegenieAnnot() {

		}

	protected abstract void dump(final SortingCollection<Variation> sorter,final VariantContext ctx) throws Exception;
		
	
	protected String fixContig(final String ctg) {
		//if(ctg.equals("X") || ctg.equals("chrX")) return "23";
		//if(ctg.equals("Y") || ctg.equals("chrY")) return "24";
		if(ctg.startsWith("chr")) return ctg.substring(3);
		return ctg;
		}

	private boolean keepVariant(final double max_freq,final VariantContext ctx) {
		if(ctx.hasAttribute(GNOMAD_AF) && ctx.getAttributeAsDouble(GNOMAD_AF, 0.0) > max_freq) return false;
		if(ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) && ctx.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0) > max_freq) return false;
		return true;
		}

	protected double getCaddScore(final VariantContext ctx) {
		Double cadd_phred = null;
		if (ctx.hasAttribute(CADD_PHRED)) {
				final String s = ctx.getAttributeAsString(CADD_PHRED, ".");
				if (!(s.equals(".") || StringUtils.isBlank(s))) {
					cadd_phred = Double.valueOf(s);
				}
			}
		return cadd_phred==null?0.0:cadd_phred.doubleValue();
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

	
	protected VCFHeader initVcfHeader(VCFHeader h) {
		return h;
	}
	
	protected abstract Logger getLogger();
	
	@Override
	public int doWork(List<String> args) {
		try {
			final double freq = Arrays.stream(this.freqStr.trim().split("[,]")).mapToDouble(S->Double.parseDouble(S)).max().orElse(0.01);
			
			
			final SortingCollection<Variation> sorter = SortingCollection.newInstance(Variation.class, new RowCodec(),
					(A, B) -> compareVariant(A, B, 2), writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths());
			sorter.setDestructiveIteration(true);


			try (VCFIterator iter = new VCFIteratorBuilder().open(System.in)) {
				initVcfHeader(iter.getHeader());
				
				while (iter.hasNext()) {
					final VariantContext vc = iter.next();
					if (vc.getNAlleles() != 2)
						throw new IOException(vc.getContig() + ":" + vc.getStart() + ":" + vc.getAlleles());
					if(!keepVariant(freq,vc)) continue;
					dump(sorter, vc);
				} // end while
			}
			sorter.doneAdding();

			
			final Set<String> seen_predictions=new HashSet<>();
			try (CloseableIterator<Variation> iter0 = sorter.iterator()) {
				try (EqualRangeIterator<Variation> iter = new EqualRangeIterator<>(iter0, (A, B) -> compareVariant(A, B, 1))) {

					try (PrintWriter annotOut = IOUtils.openPathForPrintWriter(this.annotationFileOut)) {
						try (PrintWriter setFileOut = IOUtils.openPathForPrintWriter(this.setListFileOut)) {

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
									Variation v1 = gene_variants.get(x  );
									Variation v2 = gene_variants.get(x+1);
									if(v1.id.equals(v2.id)) {
										gene_variants.remove(x+1);
									} else {
										x++;
									}
								}
								
								final Variation first = gene_variants.get(0);

								for (Variation v : gene_variants) {
									annotOut.print(v.id);
									annotOut.print(" ");
									annotOut.print(v.gene);
									annotOut.print(" ");
									annotOut.print(v.prediction);
									annotOut.print(" ");
									annotOut.print(v.score);
									annotOut.print(" ");
									annotOut.print(v.cadd);
									annotOut.println();
									
									seen_predictions.add(v.prediction);
									}
								setFileOut.print(first.gene);// gene
								setFileOut.print("\t");
								setFileOut.print(first.contig);// contig
								setFileOut.print("\t");
								setFileOut.print((int) gene_variants.stream().mapToInt(it -> it.pos).average().getAsDouble());
								setFileOut.print("\t");
								setFileOut.print(gene_variants.stream().map(it -> it.id).collect(Collectors.joining(",")));
								setFileOut.println();

							}
							setFileOut.flush();
						}
						annotOut.flush();
					}
				}

			}

			sorter.cleanup();
			
			try (PrintWriter maskOut = IOUtils.openPathForPrintWriter(this.maskFileOut)) {
				final Set<String> distinct_masks = name2prediction.values().stream().flatMap(P->P.masks.stream()).collect(Collectors.toSet());
				
				for (final String mask_name : distinct_masks) {
					final Set<String> predset = seen_predictions.
							stream().
							filter(P->name2prediction.get(P).masks.contains(mask_name)).
							map(P->name2prediction.get(P).name).
							collect(Collectors.toSet())
							;
					if(predset.isEmpty()) continue;
					maskOut.print(mask_name);
					maskOut.print("\t");
					maskOut.print(String.join(",",predset));
					maskOut.println();
					}
				}
			return 0;
		} catch (Throwable err) {
			getLogger().error(err);
			return -1;
		}
	}

}
