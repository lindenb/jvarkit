package com.github.lindenb.jvarkit.tools.rvtests;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;

public class RVTestsMakeSetFiles extends Launcher {
	private static final Logger LOG = Logger.of(RVTestsMakeSetFiles.class);

	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();


	
	private static class Variation {
		int tid;
		int start;
		int end;
		Allele ref;
		Allele alt;
		String key;

		
		int compareVariantOnKey(Variation v) {
			return this.key.compareTo(v.key);
			}
		
		int compareVariantOnKeyAndPosition(Variation v) {
			int i = compareVariantOnKey(v);
			if (i != 0) return i;
			i =Integer.compare(this.tid, v.tid);// CHROMOSOME
			if (i != 0) return i;
			i = Integer.compare(this.start, v.start);// CHROMOSOME
			if (i != 0) return i;
			return Integer.compare(this.end, v.end);
		 	}
		
		@Override
		public String toString() {
			return tid+":"+start+":"+ref.getDisplayString()+":"+alt.getDisplayString();
			}
		}



	private static class RowCodec extends AbstractDataCodec<Variation> {
		@Override
		public void encode(DataOutputStream dos, final Variation variant) throws IOException {
			dos.writeInt(variant.tid);
			dos.writeInt(variant.start);
			dos.writeInt(variant.end);
			dos.writeUTF(variant.ref.getDisplayString());
			dos.writeUTF(variant.alt.getDisplayString());
			dos.writeUTF(variant.key);
			}
		
		@Override
		public Variation decode(DataInputStream dis) throws IOException {
			final Variation v = new Variation();
			try {
				v.tid = dis.readInt();
			} catch (EOFException err) {
				return null;
			}
			v.start = dis.readInt();
			v.end = dis.readInt();
			v.ref = Allele.create( dis.readUTF(),true);
			v.alt = Allele.create( dis.readUTF(),false);
			v.key = dis.readUTF();
			return v;
			}

		@Override
		public RowCodec clone() {
			return new RowCodec();
		}
	}
	
	protected String fixContig(final String ctg) {
		if(ctg.startsWith("chr")) return ctg.substring(3);
		return ctg;
		}

	
	public RVTestsMakeSetFiles() {
		}


	

	
	@Override
	public int doWork(List<String> args) {
		try {
			final SortingCollection<Variation> sorter = SortingCollection.newInstance(Variation.class, new RowCodec(),
					(A, B) -> A.compareVariantOnKeyAndPosition(B),
					writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths());
			sorter.setDestructiveIteration(true);

			
			try (CloseableIterator<Variation> iter0 = sorter.iterator()) {
				try (EqualRangeIterator<Variation> iter = new EqualRangeIterator<>(iter0, (A, B) -> A.compareVariantOnKey(B))) {
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
			
			return 0;
			}
		catch(Throwable err ) {
			LOG.error(err);
			return -1;
			}
		}

	public static void main(String[] args) {
		new RVTestsMakeSetFiles().instanceMainWithExit(args);
	}
}
