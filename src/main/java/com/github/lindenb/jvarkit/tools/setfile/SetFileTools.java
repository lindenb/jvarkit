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
package com.github.lindenb.jvarkit.tools.setfile;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.util.IntervalExtender;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.setfile.SetFileReaderFactory;
import com.github.lindenb.jvarkit.setfile.SetFileRecord;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC


## action = combine
```
$ echo -e "w RF01:150-200\nx RF01:1-100,RF02:1-100\ny RF01:90-200,RF03:1-100\nz RF03:50-150,RF04:100-200" | java -jar dist/setfiletools.jar -R src/test/resources/rotavirus_rf.fa combine
w_x	RF01:1-100,RF01:150-200,RF02:1-100
w_y	RF01:90-200,RF03:1-100
w_z	RF01:150-200,RF03:50-150,RF04:100-200
x_y	RF01:1-200,RF02:1-100,RF03:1-100
x_z	RF01:1-100,RF02:1-100,RF03:50-150,RF04:100-200
y_z	RF01:90-200,RF03:1-150,RF04:100-200
```

## action = bedbed

creates a setfile from two bed files. First bed is "peaks.bed" second bed is "genes.bed".
The output is a set file . Each record in the output setFile is a 'gene' where all items are the overlapping peaks. 
```
gunzip -c in.gtf.gz|\
awk -F '\t' '($3=="gene") {B=int($4)-1;X=100;printf("%s\t%d\t%s\n",$1,(B<X?0:B-X),int($5)+X);}' |\
sort -t $'\t' -k1,1 -k2,2n  |\
java -jar dist/setfiletools.jar -R ref.dict  in.bed.gz -
```

## action = intersectbed

print  whole setrecords overlapping bed file, there is to trimming

```
java -jar dist/setfiletools.jar -R ref.dict  in.bed.gz in.setfile
```


END_DOC

**/
@Program(name="setfiletools",
description="Utilities for the setfile format",
creationDate="20210125",
modificationDate="20220426",
keywords={"setfile","bed"},
jvarkit_amalgamion = true,
menu="BED Manipulation"
)
public class SetFileTools extends Launcher {
	private static final Logger LOG = Logger.build(SetFileTools.class).make();
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxRef = null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT+ ". For action=cluster, output is: "+ArchiveFactory.OPT_DESC)
	protected Path outputFile=null;
	@Parameter(names={"-t","--trim-chr"},description="Remove chr prefix in chromosome names on output.")
	protected boolean trim_chr_prefix = false;
	@Parameter(names={"-U","--remove-unused-interval"},description="Remove ")
	protected boolean remove_unused_interval = false;
	@Parameter(names={"--bed"},description="Restrict input to this bed file.")
	protected Path intersectBedPath = null;
	@Parameter(names={"--vcf","--vcfs"},description="Restrict input to thoses vcf file(s). A file with the '.list' suffix is interpreted as a list of paths to the vcfs.")
	protected List<String> intersectVcfPath = new ArrayList<>();
	@Parameter(names={"--stringency"},description="Validation Stringency")
	protected ValidationStringency validationStringency = ValidationStringency.LENIENT;
	@Parameter(names={"--min-variants-per-setfile"},description="when using vcf, only keep the setfile is there are at least 'x' overlapping variants.")
	private int min_variant_per_setfile = 1;
	@Parameter(names={"--extend"},description="Extends each interval. "+IntervalExtender.OPT_DESC)
	private String extendStr = null;
	@Parameter(names={"--disable-uniq"},description="disable unique-name checker.")
	private boolean disable_uniq = false;

	
	/** global fai */
	private SAMSequenceDictionary theDict = null;
	/** global sorter */
	private Comparator<Locatable> theSorter;
	
	/** check uniq names */
	private class UniqNameIterator extends AbstractCloseableIterator<SetFileRecord> {
		final CloseableIterator<SetFileRecord> delegate;
		final Set<String> names = new HashSet<>(10_000);
		
		UniqNameIterator(final CloseableIterator<SetFileRecord> delegate) {
			this.delegate = delegate;
			}
		
		@Override
		protected final SetFileRecord advance() {
			if(!this.delegate.hasNext()) return null;
			final SetFileRecord rec = this.delegate.next();
			if(!this.names.add(rec.getName())) {
				throw new RuntimeIOException("setfile contains duplicate name \""+rec.getName()+"\".");
				}
			return rec;
			}

		@Override
		public void close() {
			delegate.close();
			}
		}
	
	/** iterator extending+merging intervals */
	private class ExtenderIterator extends AbstractCloseableIterator<SetFileRecord> {
		final IntervalExtender xtender;
		final CloseableIterator<SetFileRecord> delegate;
		
		ExtenderIterator(final CloseableIterator<SetFileRecord> delegate,final IntervalExtender xtender) {
			this.delegate = delegate;
			this.xtender = xtender;
			if(this.xtender.isShriking()) throw new IllegalArgumentException("shrinking is not supported:"+ extendStr);
			}
		
		@Override
		protected final SetFileRecord advance() {
			while(this.delegate.hasNext()) {
				final SetFileRecord rec = this.delegate.next();
				final List<Locatable> L = new ArrayList<>(rec.size());
				for(int i=0;i< rec.size();i++) {
					final Locatable loc0 = rec.get(i);
					final Locatable loc = this.xtender.apply(loc0);
					L.add(loc);
					}
				return SetFileRecord.create(rec.getName(),sortAndMerge(L));
				}
			return null;
			}

		
		@Override
		public void close() {
			delegate.close();
			}
		}
	
	/** iterator filtering on VCF */
	private class IntersectVcfIterator extends AbstractCloseableIterator<SetFileRecord> {
		final List<BufferedVCFReader> vcfReaders = new ArrayList<>();
		final List<UnaryOperator<String>> contigConverters = new ArrayList<>();
		final CloseableIterator<SetFileRecord> delegate;
		IntersectVcfIterator(final CloseableIterator<SetFileRecord> delegate,final List<Path> paths) {
			this.delegate = delegate;
			final VCFReaderFactory vrf = VCFReaderFactory.makeDefault();
			final UnaryOperator<VariantContext> simplifier = V->new VariantContextBuilder(V).
					noGenotypes().
					noID().
					unfiltered().
					attributes(Collections.emptyMap()).
					make();
			for(final Path path:paths) {
				@SuppressWarnings("resource")
				final BufferedVCFReader vr=new BufferedVCFReader(vrf.open(path,true),10_000).setSimplifier(simplifier);
				this.vcfReaders.add(vr);
				final VCFHeader h = vr.getHeader();
				this.contigConverters.add(ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(h)));
				}
			}
		
		@Override
		protected final SetFileRecord advance() {
			while(this.delegate.hasNext()) {
				final SetFileRecord rec = this.delegate.next();
				final List<Locatable> L = new ArrayList<>(rec.size());
				long n = 0L;
				for(int i=0;i< rec.size();i++) {
					final Locatable loc = rec.get(i);
					long n2 = 0L;
					for(int k=0;k< this.vcfReaders.size();++k) {
						final BufferedVCFReader vr = this.vcfReaders.get(k);
						final String ctg = this.contigConverters.get(k).apply(loc.getContig());
						if(StringUtil.isBlank(ctg)) continue;
						try(CloseableIterator<VariantContext> iter2 = vr.query(new SimpleInterval(ctg,loc.getStart(),loc.getEnd()))) {
							n2 += iter2.stream().filter(V->V.overlaps(loc)).count();
							}
						}
					if(n2==0L && remove_unused_interval) {
						continue;
						}
					L.add(loc);
					n+=n2;
					}
				if(L.isEmpty() || n< min_variant_per_setfile) continue;
				return SetFileRecord.create(rec.getName(), L);
				}
			return null;
			}

		
		@Override
		public void close() {
			delegate.close();
			for(BufferedVCFReader b:this.vcfReaders) try {b.close();}
			catch(final IOException err) {}
			}
		}
	
	/** open path or null(stdin) */
	private BedLineReader openBedLineReader(final Path pathOrNull) {
		try
			{
			BedLineReader br;
			if(pathOrNull!=null) {
				br = new BedLineReader(pathOrNull);
				}
			else {
				br = new BedLineReader(stdin(),"stdin");
				}
			br.setContigNameConverter(ContigNameConverter.fromOneDictionary(SetFileTools.this.theDict));
			br.setValidationStringency(validationStringency);
			return br;
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	
	private class IntersectBedIterator extends AbstractCloseableIterator<SetFileRecord> {
		final IntervalTreeMap<Boolean> intervalTreeMap;
		final CloseableIterator<SetFileRecord> delegate;
		IntersectBedIterator(final CloseableIterator<SetFileRecord> delegate,final Path bedPath) {
			this.delegate = delegate;
			try(BedLineReader br = openBedLineReader(bedPath)) {
				this.intervalTreeMap = br.toIntervalTreeMap(X->Boolean.TRUE);
				}
			}
		
		@Override
		protected final SetFileRecord advance() {
			while(this.delegate.hasNext()) {
				final SetFileRecord rec = this.delegate.next();
				final List<Locatable> L = new ArrayList<>(rec.size());
				for(int i=0;i< rec.size();i++) {
					final Locatable loc = rec.get(i);
					boolean keep = this.intervalTreeMap.containsOverlapping(loc);
					if(!keep && remove_unused_interval) {
						continue;
						}
					L.add(loc);
					}
				if(L.isEmpty()) continue;
				return SetFileRecord.create(rec.getName(), L);
				}
			return null;
			}

		
		@Override
		public void close() {
			delegate.close();
			}
		}
	

	private enum Action {
		view,
		combine,
		bedbed,
		intersectbed,
		stats,
		};
		
	private String noChr(final String contig) {
		if(trim_chr_prefix && contig.toLowerCase().startsWith("chr")) {
			return contig.substring(3);
		}
		return contig;
	}
	
	private List<Locatable> sortAndMerge(final List<? extends Locatable> L0) {
			final List<Locatable> L = new ArrayList<>(L0);
			// merge overlapping
			Collections.sort(L, theSorter);
			int i=0;
			while(i+1 < L.size()) {
				final Locatable xi = L.get(i  );
				final Locatable xj = L.get(i+1);
				if(xi.overlaps(xj)) {
					L.set(i, new SimpleInterval(
							xi.getContig(),
							Math.min(xi.getStart(),xj.getStart()),
							Math.max(xi.getEnd(),xj.getEnd())
							));
					L.remove(i+1);
				} else {
					i++;
				}
			}
		return L;
		}
	
	private CloseableIterator<SetFileRecord> openSetFileIterator(final List<String> args) throws IOException {
		CloseableIterator<SetFileRecord> iter  = null;
		final String input = oneFileOrNull(args);
		final SetFileReaderFactory srf  = new SetFileReaderFactory(this.theDict);
		if(input==null) {
			iter = srf.open(IOUtils.openStdinForBufferedReader());
		} else {
			iter  =srf.open(IOUtils.openURIForBufferedReading(input));
		}
		if(!StringUtils.isBlank(this.extendStr)) {
			final IntervalExtender xtExtender =  IntervalExtender.of(this.theDict, this.extendStr);
			if(xtExtender.isChanging()) {
				iter  = new ExtenderIterator(iter, xtExtender);
				}
		}
		
		
		if(intersectBedPath!=null) {
			iter  = new IntersectBedIterator(iter, this.intersectBedPath);
			}
		if(!intersectVcfPath.isEmpty()) {
			iter = new IntersectVcfIterator(iter, IOUtils.unrollPaths(this.intersectVcfPath));
			}
		if(!disable_uniq) {
			iter  = new UniqNameIterator(iter);
			}
		return iter;
		}
	
	/** print SetFileRecord to pw */
	private void print(PrintWriter pw,SetFileRecord setfile) {
		if(setfile.isEmpty()) return;
		pw.write(setfile.getName());
		for(int i=0;i< setfile.size();i++) {
			final Locatable rec = setfile.get(i);
			pw.write(i==0?"\t":",");
			pw.write(noChr(rec.getContig()));
			pw.write(":");
			pw.write(String.valueOf(rec.getStart()));
			pw.write("-");
			pw.write(String.valueOf(rec.getEnd()));
			}
		pw.write("\n");
	}
	

	
	private int view(final List<String> args) throws IOException {
		try(CloseableIterator<SetFileRecord> iter = openSetFileIterator(args)) {
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				while(iter.hasNext()) {
					final SetFileRecord rec = iter.next();
					print(pw,rec);
					}
				pw.flush();
				}
			}
		return 0;
		}
	
	private int combine(final List<String> args) throws IOException {
	try(CloseableIterator<SetFileRecord> iter = openSetFileIterator(args)) {
		final List<SetFileRecord> L2 = iter.stream().collect(Collectors.toList());
		try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
			for(int i=0;i+1< L2.size();i++) {
				final SetFileRecord r1 = L2.get(i);
				for(int j=i+1;j< L2.size();j++) {
					final SetFileRecord r2 = L2.get(j);
					final List<Locatable> L = new ArrayList<>(r1.size()+r2.size());
					L.addAll(r1);
					L.addAll(r2);
					final String name = String.join("_", r1.getName(),r2.getName());
					print(pw,SetFileRecord.create(name, sortAndMerge(L)));
					}
				}
			pw.flush();
			}
		}
	return 0;
	}
	
	
	
	
	
	
	
	/** print  whole setrecords overlapping bed file, there is to trimming */
	private int doInterBed(final List<String> args) throws IOException {
		if(this.intersectBedPath!=null) {
			LOG.info("intersectBedPath shouldn' be specified");
			return -1;
		}
		if(!this.intersectVcfPath.isEmpty()) {
			LOG.info("intersectVcfPath shouldn't be specified");
			return -1;
		}
		if(args.size()!=2) {
			LOG.error("expected 2 files but got "+args.size()+" "+args);
			return -1;
		}
		if(args.get(0).equals("-") && args.get(1).equals("-")) {
			LOG.error("cannot use both files on stdin");
			return -1;
		}
		
		final IntervalTreeMap<BedLine> peaksTreeMap;
		try(BedLineReader blr = openBedLineReader(args.get(0).equals("-")?null:Paths.get(args.get(0)))) {
			peaksTreeMap = blr.toIntervalTreeMap();
			}
		
		try(CloseableIterator<SetFileRecord> iter = openSetFileIterator((args.get(1).equals("-")?Collections.emptyList():args.subList(1,2)))) {
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				while(iter.hasNext()) {
					final SetFileRecord rec = iter.next();
					if(rec.stream().noneMatch(B->peaksTreeMap.containsOverlapping(B))) continue;
					print(pw,rec);
					}
				pw.flush();
				}
			}
		return 0;
		}
	
	/** statistics for setFile*/
	private int doStats(final List<String> args) throws IOException {
		final Counter<String> chrom2count=new Counter<>();
		final DiscreteMedian<Integer> d_size = new DiscreteMedian<>();
		final DiscreteMedian<Integer> d_nitems = new DiscreteMedian<>();
		final DiscreteMedian<Integer> d_distance = new DiscreteMedian<>();
		final DiscreteMedian<Integer> d_item_size = new DiscreteMedian<>();
		for(final SAMSequenceRecord ssr: this.theDict.getSequences()) {
			chrom2count.initializeIfNotExists(noChr(ssr.getSequenceName()));
			}
		chrom2count.initializeIfNotExists("*multiple*");
		chrom2count.initializeIfNotExists("*empty*");

		try(CloseableIterator<SetFileRecord> iter = openSetFileIterator(args)) {
			while(iter.hasNext()) {
				final SetFileRecord rec = iter.next();
				final Set<String> chroms = rec.getChromosomes();
				switch(chroms.size()) {
					case 0: chrom2count.incr("*empty*"); break;
					case 1: chrom2count.incr(noChr(chroms.iterator().next())); break;
					default: chrom2count.incr("*multiple*"); break;
					}
				int len = rec.stream().mapToInt(B->B.getLengthOnReference()).sum();
				d_size.add(len);
				d_nitems.add(rec.size());
				if(rec.size()>0) {
					len = len/rec.size();
					d_item_size.add(len);
				}
				
				if(rec.size()>1 && chroms.size()==1) {
					int d = 0;
					final List<Locatable> L = sortAndMerge(rec);
					for(int i=0;i+1<L.size();i++) {
						d+= (rec.get(i+1).getStart() - rec.get(i).getEnd());
						}
					d=d/(L.size()-1);
					d_distance.add(d);
					}
				}
			}
		try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
			for(final String key: chrom2count.keySetDecreasing()) {
				pw.println("C\trecords-per-chrom\t"+key+"\t"+chrom2count.count(key));
			}
			pw.println("AS\taverage-size\t"+(d_size.isEmpty()?".":String.valueOf(d_size.getAverage().orElse(0.0))));
			pw.println("MS\tmedian-size\t"+(d_size.isEmpty()?".":String.valueOf(d_size.getMedian().orElse(0.0))));
			
			pw.println("AIS\taverage-item-size\t"+(d_item_size.isEmpty()?".":String.valueOf(d_item_size.getAverage().orElse(0.0))));
			pw.println("MIS\tmedian-item-size\t"+(d_item_size.isEmpty()?".":String.valueOf(d_item_size.getMedian().orElse(0.0))));

			
			pw.println("AN\taverage-nitems\t"+(d_nitems.isEmpty()?".":String.valueOf(d_nitems.getAverage().orElse(0.0))));
			pw.println("MN\tmedian-nitems\t"+(d_nitems.isEmpty()?".":String.valueOf(d_nitems.getMedian().orElse(0.0))));
			
			pw.println("AD\taverage-distance-between-items\t"+(d_distance.isEmpty()?".":String.valueOf(d_distance.getAverage().orElse(0.0))));
			pw.println("MD\tmedian-distance-between-items\t"+(d_distance.isEmpty()?".":String.valueOf(d_distance.getMedian().orElse(0.0))));
			
			pw.flush();	
			}
		return 0;
		}

	
	/** create setfile each line is a BED (file2) containing the components of overlapping from BED (file1) */
	private int doBedBed(final List<String> args) throws IOException {
		if(this.intersectBedPath!=null) {
			LOG.info("intersectBedPath shouldn' be specified");
			return -1;
		}
		if(!this.intersectVcfPath.isEmpty()) {
			LOG.info("intersectVcfPath shouldn't be specified");
			return -1;
		}
		if(args.size()!=2) {
			LOG.error("expected 2 files but got "+args.size()+" "+args);
			return -1;
		}
		if(args.get(0).equals("-") && args.get(1).equals("-")) {
			LOG.error("cannot use both files on stdin");
			return -1;
		}
		
		final IntervalTreeMap<BedLine> peaksTreeMap;
		try(BedLineReader blr = openBedLineReader(args.get(0).equals("-")?null:Paths.get(args.get(0)))) {
			peaksTreeMap = blr.toIntervalTreeMap();
			}
		
		long id_generator=0L;
		try(BedLineReader iter = openBedLineReader(args.get(1).equals("-")?null:Paths.get(args.get(1)))) {
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				while(iter.hasNext()) {
					final BedLine rec= iter.next();
					final List<Locatable> L = sortAndMerge(peaksTreeMap.
							getOverlapping(rec).stream().
							sorted(theSorter).
							collect(Collectors.toList()));
					if(L.isEmpty()) continue;
					final String id = String.format("%s_%d_%d_%d",
						noChr(rec.getContig()),
						rec.getStart(),
						rec.getEnd(),
						++id_generator
						);
					print(pw,SetFileRecord.create(id, L));
					}
				pw.flush();
				}
			}

		
		return 0;
	}
	
	@Override
	public int doWork(final List<String> args0) {
		if(args0.isEmpty()) {
			LOG.error("action parameter is missing.");
			return -1;
			}
		try {
			this.theDict= SequenceDictionaryUtils.extractRequired(this.faidxRef);
			this.theSorter = new ContigDictComparator(this.theDict).createLocatableComparator();
			
			final Action action = Action.valueOf(args0.get(0));
			final List<String> args = args0.subList(1, args0.size());
			switch(action) {
				case view: return view(args);
				case combine: return combine(args);
				case bedbed: return doBedBed(args);
				case intersectbed:  return doInterBed(args);
				case stats:  return doStats(args);
				default: LOG.error("not implemented "+action);return -1;
				}
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new SetFileTools().instanceMainWithExit(args);
	}

}
