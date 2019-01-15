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
package com.github.lindenb.jvarkit.tools.lumpysv;

import java.io.BufferedReader;
import java.io.File;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;


/**
BEGIN_DOC

## Input

input is a set of VCF file or a file ending with '.list' and containing one line per VCF path.

## Example

```
$ find DIR -name "*.vcf" > vcf.list
$ mkdir -p BDB
$ java -jar dist/lumpysort.jar --bdb BDB/  vcf.list > merged.ved
```

END_DOC
 */
@Program(name="lumpysort",
description="sort and merge a set of Lumpy-SV VCF files.",
keywords={"lumpy","vcf","sort"},
generate_doc=false
)
public class LumpySort 
	 extends Launcher {
	
	private static final Logger LOG = Logger.build(LumpySort.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-f","--fraction"},description="Required Overlap fraction between two intervals.")
	private double fraction_overlap = 0.8;
	@Parameter(names={"-slop","--slop"},description="slop for 'x' bases in both directions.")
	private int slop_size = 0;
	@Parameter(names={"-secondary","--secondary"},description="keep SECONDARY record. \"Secondary breakend in a multi-line variants\".")
	private boolean keep_secondary = false;
	@Parameter(names={"-vf","--variant-filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = JexlVariantPredicate.create("");
	@Parameter(names={"-dm","--do-not-merge"},description="Do not merge genotypes, just sort")
	private boolean do_not_merge_ctx = false;
	@Parameter(names={"-gt","--genotype"},description="genotypes with SU>0 will be called with ALT (haploid)")
	private boolean do_genotype = false;
	@Parameter(names={"-B","--bed"},description="restrict to variants overlapping this BED file.")
	private File bedFile = null;
	@Parameter(names={"-bdb","--bdb"},description="Berkeley DB working directory. It must exists and must be writeable.",required=true)
	private File bdbHomeDir = null;
	
	/** encoder for VariantCtx -> line */
	private VCFEncoder vcfEncoder = null;
	/** encoder for line -> VariantCtx */
	private final VCFCodec vcfCodec = new VCFCodec();
	
	
	
	/** berkeleyDB LumpyVar encoder */
	private  class LumpyVarBinding extends TupleBinding<LumpyVar>
		{
		@Override
		public LumpyVar entryToObject(final TupleInput in) {
			long id=in.readLong();
			int L=in.readInt();
			byte array[]=new byte[L];
			in.readFast(array);
			final LumpyVar v= new LumpyVar(linetoVariantContext(new String(array)),id);
			return v;
			}
		@Override
		public void objectToEntry(final LumpyVar v,final TupleOutput out) {
			out.writeLong(v.id);
			final String line = variantContextToLine(v.ctx);
			final byte array[]=line.getBytes();
			out.writeInt(array.length);
			out.write(array);
			}
		}
	
	

	/** berkeley DB 'key' object */
	private static class KeySorter
		{
		StructuralVariantType st;
		String contig1;
		String contig2;/* for BnD  variant */
		int start;
		int end;
		long id;
		public int compare1(final KeySorter o) {
			int i= st.compareTo(o.st);
			if(i!=0) return i;
			//
			i = this.contig1.compareTo(o.contig1);
			if(i!=0) return i;
			//
			if(st.equals(StructuralVariantType.BND)) {
				i = this.contig2.compareTo(o.contig2);
				if(i!=0) return i;
				}
			return 0;
			}
		public int compare2(final KeySorter o) {
			int i= compare1(o);
			if(i!=0) return i;
			//
			i = this.start -o.start;
			if(i!=0) return i;
			//
			i = this.end -o.end;
			if(i!=0) return i;
			return Long.compare(id, o.id);
			}
		}
	
	/** berkeley DB 'key' comparator */
	public static class KeySorterComparator implements Comparator<byte[]>
		{
		private final KeySorterBinding binding = new KeySorterBinding();
		private KeySorter make(final byte[] o1)
			{
			final TupleInput in = new TupleInput(o1);
			final KeySorter ks = this.binding.entryToObject(in);
			try{in.close();}catch(Exception err) {}
			return ks;
			}
		@Override
		public int compare(final byte[] o1, final byte[] o2) {
			return make(o1).compare2(make(o2));
			}
		}
	
	/** berkeley DB 'key' bdinding -> DataEntry */
	private static class KeySorterBinding extends TupleBinding<KeySorter>
		{
		private final StructuralVariantType sttypes[]=StructuralVariantType.values();
		@Override
		public void objectToEntry(final KeySorter key, final TupleOutput out) {
			out.writeByte((byte)key.st.ordinal());
			out.writeString(key.contig1);
			if(key.st.equals(StructuralVariantType.BND)) {
				out.writeString(key.contig2);
				}
			out.writeInt(key.start);
			out.writeInt(key.end);
			out.writeLong(key.id);
			}
		@Override
		public KeySorter entryToObject(final TupleInput in) {
			final KeySorter ks= new KeySorter();
			ks.st = this.sttypes[(int)in.readByte()];
			ks.contig1= in.readString();
			if(ks.st.equals(StructuralVariantType.BND)) {
				ks.contig2= in.readString();
				}
			else
				{
				ks.contig2= ks.contig1;
				}
			ks.start = in.readInt();
			ks.end = in.readInt();
			ks.id = in.readLong();
			return ks;
			}	
		}
	
	
	/** VariantContext wrapper */
	private class LumpyVar
		{
		private final VariantContext ctx;
		private final long id;
		private Interval _interval=null;
		private Interval _bndinterval=null;
		private KeySorter _sortKey = null;
		
		LumpyVar(final VariantContext ctx,final long id)
			{
			this.ctx = ctx;
			this.id=id;
			}
		private Set<String> getGenotypedSamples() {
			return ctx.getGenotypes().stream().
					filter(G->isAvailableGenotype(G)).
					map(G->G.getSampleName()).
					collect(Collectors.toSet());
			}	
		
		private KeySorter getSortKey() {
			if(_sortKey==null) {
				final Function<String,String> normalize=C->C.startsWith("chr")?C.substring(3):C;
				final Interval rgn = getInterval() ;
				_sortKey = new KeySorter();
				_sortKey.st = this.ctx.getStructuralVariantType();
				_sortKey.contig1 = normalize.apply(rgn.getContig());
				if(_sortKey.st.equals(StructuralVariantType.BND)) {
					_sortKey.contig2 = normalize.apply(LumpyConstants.getBnDContig(this.ctx.getAlternateAllele(0).getDisplayString()));
					}
				else
					{
					_sortKey.contig2 = _sortKey.contig1;
					}
				_sortKey.start = rgn.getStart();
				_sortKey.end = rgn.getEnd();
				_sortKey.id=this.id;
				}
			return _sortKey;
			}
		
		private Interval getInterval() {
			if(this._interval==null) {
				if(!ctx.hasAttribute("CIPOS")) throw new IllegalArgumentException("No CIPOS in "+ctx);
				final List<Integer> ciposL= ctx.getAttributeAsIntList("CIPOS",0);
				if(ciposL.size()!=2) throw new IllegalArgumentException("len(CIPOS)!=2 in "+ctx);
				if(!ctx.hasAttribute("CIEND")) throw new IllegalArgumentException("No CIEND in "+ctx);
				final List<Integer> ciendL= ctx.getAttributeAsIntList("CIEND",0);
				if(ciendL.size()!=2) throw new IllegalArgumentException("len(CIEND)!=2 in "+ctx);
	
				this._interval =  new Interval(
					ctx.getContig(),
					Math.max(0,ctx.getStart() + ciposL.get(0) - LumpySort.this.slop_size),
					ctx.getEnd() + ciendL.get(1) + LumpySort.this.slop_size
					);
				}
			return this._interval;
			}
		private Interval getBndInterval() {
			if(this._bndinterval==null) {
				if(!ctx.hasAttribute("CIPOS")) throw new IllegalArgumentException("No CIPOS in "+ctx);
				final List<Integer> ciposL= ctx.getAttributeAsIntList("CIPOS",0);
				if(ciposL.size()!=2) throw new IllegalArgumentException("len(CIPOS)!=2 in "+ctx);
				if(!ctx.hasAttribute("CIEND")) throw new IllegalArgumentException("No CIEND in "+ctx);
				final List<Integer> ciendL= ctx.getAttributeAsIntList("CIEND",0);
				if(ciendL.size()!=2) throw new IllegalArgumentException("len(CIEND)!=2 in "+ctx);
				
				String cL;
				int pL;
				if(ctx.getStructuralVariantType()==StructuralVariantType.BND) {
					final  Map.Entry<String,Integer> entry = LumpyConstants.getBnDContigAndPos(ctx.getAlternateAllele(0).getDisplayString());
					cL = entry.getKey();
					pL = entry.getValue();
					} 
				else
					{
					cL = ctx.getContig();
					pL = ctx.getEnd();
					}
			
				this._bndinterval =  new Interval(
					cL,
					Math.max(0,pL+ciposL.get(0) - LumpySort.this.slop_size),
					pL+ciendL.get(1) + LumpySort.this.slop_size
					);
				}
			return this._bndinterval;
			}
		
		boolean canMerge(final LumpyVar o)
			{
			// we cannot have common available variants between two ctx
			final Set<String> commonSamples= new HashSet<String>(this.getGenotypedSamples());
			commonSamples.retainAll(o.getGenotypedSamples());
			
			
			if(!commonSamples.isEmpty()) {
				return false;
			}
			
			Interval L1 = this.getInterval();
			Interval L2 = o.getInterval();
			if(!LumpySort.this.overlap(L1,L2)) return false;
			if(this.ctx.getStructuralVariantType()==StructuralVariantType.BND) {
				L1 = this.getBndInterval();
				L2 = o.getBndInterval();
				if(!LumpySort.this.overlap(L1,L2)) return false;
				}
			
			return true;
			}
		
		
		}
	
	
	/** variant decoder */
	private VariantContext linetoVariantContext(final String line) {
		return this.vcfCodec.decode(line);
	}
	
	/** encoder variant decoder */
	private  String variantContextToLine(final VariantContext ctx) {
		return this.vcfEncoder.encode(ctx);
	}
	
	/** returns true two interval overlap with fraction_overlap  */
	private boolean overlap(final Interval i1,final Interval i2)
		{
		if(!i1.intersects(i2)) return false;
		final int L1 = i1.length();
		final int L2 = i2.length();
		final int  L3 = i1.getIntersectionLength(i2);
		if(L3< (int)(this.fraction_overlap*L1)) return false;
		if(L3< (int)(this.fraction_overlap*L2)) return false;
		return true;
		}
	
	/** returns true if there is a SU greater than 0 */
	private boolean isAvailableGenotype(final Genotype g)
		{
		@SuppressWarnings("deprecation")
		final int suv  = g.getAttributeAsInt("SU", 0);
		if(suv<=0)
			{
			return false;
			}
		return true;
		}
	
	@Override
	public int doWork(final List<String> args) {
	VariantContextWriter vcw = null;
	LineIterator vcfIn= null;
	Environment environment = null;
	Database variantsDb1=null;
	final List<File> inputs = IOUtil.unrollFiles(
			args.stream().map(S->new File(S)).collect(Collectors.toList()),
			".vcf",".vcf.gz");
	if(inputs.isEmpty()) {
		LOG.error("empty vcf list");
		return -1;
		}
	try {
		IOUtil.assertDirectoryIsWritable(this.bdbHomeDir);

		final Set<VCFHeaderLine> metaData = new HashSet<>();
		final Set<String> sampleNames = new TreeSet<>();
		final IntervalTreeMap<Boolean> intervalTreeMapBed;
		if(this.bedFile!=null)
			{
			intervalTreeMapBed = new IntervalTreeMap<>();
			final BedLineCodec bedLineCodec = new BedLineCodec();
			final BufferedReader br = IOUtils.openFileForBufferedReading(this.bedFile);
			br.lines().
				map(L->bedLineCodec.decode(L)).
				filter(L->L!=null).
				forEach(B->intervalTreeMapBed.put(B.toInterval(),true));
			br.close();
			}	
		else
			{
			intervalTreeMapBed = null;
			}
		
		for(int idx=0;idx< inputs.size();++idx)
			{
			final File vcfFile = inputs.get(idx);
			LOG.info("Read header "+(idx+1)+"/"+inputs.size());
			final VCFFileReader r  = new VCFFileReader(vcfFile,false);
			final VCFHeader header = r.getFileHeader();
			if(!LumpyConstants.isLumpyHeader(header))
				{
				LOG.error("doesn't look like a Lumpy-SV vcf header "+vcfFile);
				r.close();
				return -1;
				}
			
			if(!header.hasGenotypingData()) {
				LOG.error("No sample in "+vcfFile);
				r.close();
				return -1;
				}
			for(final String sampleName : header.getSampleNamesInOrder())
				{
				if(sampleNames.contains(sampleName)) {
					LOG.error("Sample found twice "+sampleName+" in "+vcfFile);
					r.close();
					return -1;
					}
				sampleNames.add(sampleName);
				}
			metaData.addAll(
					header.getMetaDataInInputOrder().
					stream().filter(H->!H.getKey().equals("fileDate")).
						collect(Collectors.toSet())
					);
			r.close();
			}
		final VCFInfoHeaderLine nSampleInfoHeaderLine = new VCFInfoHeaderLine("NSAMPLES", 1, VCFHeaderLineType.Integer,"Number of affected samples.");
		metaData.add(nSampleInfoHeaderLine);
		final VCFFormatHeaderLine chromStartFormatHeaderLine = new VCFFormatHeaderLine(
				"CB", 1, VCFHeaderLineType.Integer,"Original Variant POS");
		metaData.add(chromStartFormatHeaderLine);
		final VCFFormatHeaderLine chromEndFormatHeaderLine = new VCFFormatHeaderLine(
				"CE", 1, VCFHeaderLineType.Integer,"Original Variant END");
		metaData.add(chromEndFormatHeaderLine);

		
		
		final VCFHeader outHeader = new VCFHeader(
			metaData,
			sampleNames
			);
		final VCFHeaderVersion versions[]=VCFHeaderVersion.values();
		this.vcfEncoder = new VCFEncoder(outHeader, false, true);
		this.vcfCodec.setVCFHeader(
				outHeader,
				versions[versions.length-1]
				);
		
		
		/* open BDB env */
		final Transaction txn=null;
		environment = new Environment(this.bdbHomeDir, 
				new EnvironmentConfig().
					setAllowCreate(true).
					setReadOnly(false)
			);
		
		variantsDb1 = environment.openDatabase(txn,"variants1",
			new DatabaseConfig().
				setBtreeComparator(KeySorterComparator.class).
				setAllowCreate(true).
				setReadOnly(false).
				setTemporary(true)
			);
		
		long total_variants = 0L;
		
		final LumpyVarBinding lumpVarBinding = new LumpyVarBinding();
		final KeySorterBinding keySorterBinding = new KeySorterBinding();

		for(int idx=0;idx< inputs.size();++idx)
			{
			final long millisecstart = System.currentTimeMillis();
			final File vcfFile = inputs.get(idx);
			int nVariant = 0;
			final VCFFileReader r  = new VCFFileReader(vcfFile,false);
			
			final List<Genotype> missing =new ArrayList<>(sampleNames.size());
			for(final String sn:sampleNames)
				{
				if(r.getFileHeader().getSampleNamesInOrder().contains(sn)) continue;
				missing.add(GenotypeBuilder.createMissing(sn, 2));
				}
			
			final CloseableIterator<VariantContext> iter = r.iterator();
			while(iter.hasNext()) {
				VariantContext ctx = iter.next();
				if(!this.keep_secondary) {
					if(ctx.hasAttribute("SECONDARY")) continue;
					}
				if(!this.variantFilter.test(ctx)) continue;
				
				if(intervalTreeMapBed!=null &&
					!intervalTreeMapBed.containsOverlapping(ctx)) continue;
					
				
				final List<Genotype> gtList  = new ArrayList<>(ctx.getGenotypes());
				
				for(int gi=0;gi< gtList.size();gi++)
					{
					Genotype g= gtList.get(gi);
					final GenotypeBuilder gb;
					
					if(this.do_genotype && isAvailableGenotype(g))
						{
						gb = new GenotypeBuilder(g.getSampleName(), ctx.getAlternateAlleles());
						gb.attributes(g.getExtendedAttributes());
						}
					else
						{
						gb = new GenotypeBuilder(g);
						}
					gb.attribute(chromStartFormatHeaderLine.getID(), ctx.getStart());
					gb.attribute(chromEndFormatHeaderLine.getID(), ctx.getEnd());
					gtList.set(gi, gb.make());
					}
				
				
				gtList.addAll(missing);
				
				ctx = new VariantContextBuilder(ctx).
						genotypes(gtList).
						rmAttribute("PRPOS").
						make();
				
				final LumpyVar lvar = new LumpyVar(ctx,total_variants);
				final DatabaseEntry key = new DatabaseEntry();
				final DatabaseEntry data = new DatabaseEntry();
				
				lumpVarBinding.objectToEntry(lvar, data);
				keySorterBinding.objectToEntry(lvar.getSortKey(), key);
				if(variantsDb1.put(txn, key, data)!=OperationStatus.SUCCESS)
					{
					r.close();
					LOG.error("insertion failed");
					return -1;
					}
				nVariant++;
				total_variants++;
				}
			iter.close();
			r.close();

			LOG.info("Read  "+(idx+1)+"/"+inputs.size()+" variants of "+vcfFile+" N="+nVariant+
					" Total:"+total_variants + 
					" That took: " + Duration.ofMillis(System.currentTimeMillis() -millisecstart )
					);
			System.gc();
			}
		
		if(intervalTreeMapBed!=null) intervalTreeMapBed.clear();
		System.gc();
		
		LOG.info("Writing output");
		final List<Allele> ALLELES_NO_CALLS=
				this.do_genotype
				? Collections.singletonList(Allele.NO_CALL)
				: Arrays.asList(Allele.NO_CALL,Allele.NO_CALL)
				;
		final Cursor cursor = variantsDb1.openCursor(txn, null);


		vcw = super.openVariantContextWriter(this.outputFile);
		vcw.writeHeader(outHeader);
		
		for(;;)
			{
			final DatabaseEntry key = new DatabaseEntry();
			final DatabaseEntry data = new DatabaseEntry();
			OperationStatus status = cursor.getNext(key, data, LockMode.DEFAULT);
			if(!status.equals(OperationStatus.SUCCESS)) break;
			final LumpyVar first = lumpVarBinding.entryToObject(data);
			if(this.do_not_merge_ctx)
				{
				vcw.add(first.ctx);
				continue;
				}

			final KeySorter keySorter1 = keySorterBinding.entryToObject(key);

			
			final List<LumpyVar> buffer = new ArrayList<>();
			buffer.add(first);
			
			final DatabaseEntry key2 = new DatabaseEntry();
			final DatabaseEntry data2 = new DatabaseEntry();

			final Cursor cursor2=cursor.dup(true);
			for(;;)
				{
				status = cursor2.getNext(key2, data2, LockMode.DEFAULT);
				if(!status.equals(OperationStatus.SUCCESS)) break;
				final KeySorter keySorter2 = keySorterBinding.entryToObject(key2);

				if(keySorter1.compare1(keySorter2)!=0) 
					{
					break;
					}
				
				final LumpyVar lv = lumpVarBinding.entryToObject(data2);
				if(lv.ctx.getStart()>first.ctx.getEnd()) {
					break;
					}
				if(first.canMerge(lv))
					{
					buffer.add(lv);
					cursor2.delete();
					}
				}
			cursor2.close();
			
			
			cursor.delete();//delete 'first'
			
			
			
			final int variantStartA = buffer.stream().
					mapToInt(V->V.ctx.getStart()).
					min().getAsInt();
			final int variantStartB = (int)buffer.stream().
					mapToInt(V->V.ctx.getStart()).
					average().getAsDouble();
			final int variantStartC = buffer.stream().
					mapToInt(V->V.ctx.getStart()).
					max().getAsInt();
			
			final int variantEndA = buffer.stream().
					mapToInt(V->V.ctx.getEnd()).
					min().getAsInt();
			final int variantEndB = (int)buffer.stream().
					mapToInt(V->V.ctx.getEnd()).
					average().getAsDouble();
			final int variantEndC = buffer.stream().
					mapToInt(V->V.ctx.getEnd()).
					max().getAsInt();
			
			final VariantContextBuilder vcb = new VariantContextBuilder(
					"lumpymerge",
					first.ctx.getContig(),
					variantStartB,
					variantEndB,
					first.ctx.getAlleles()
					);
			vcb.attribute("END", variantEndB);
			vcb.attribute("SVTYPE", first.ctx.getAttribute("SVTYPE"));
			vcb.attribute("SVLEN", (int)Percentile.median().evaluate(buffer.stream().mapToInt(V->V.ctx.getEnd()-V.ctx.getStart())));
			vcb.attribute("CIPOS",Arrays.asList(variantStartB-variantStartA,variantStartC-variantStartB));
			vcb.attribute("CIEND",Arrays.asList(variantEndB-variantEndA,variantEndC-variantEndB));
			vcb.attribute("SU",buffer.stream().flatMap(V->V.ctx.getGenotypes().stream()).mapToInt(G->G.getAttributeAsInt("SU", 0)).sum());
			vcb.attribute("SR",buffer.stream().flatMap(V->V.ctx.getGenotypes().stream()).mapToInt(G->G.getAttributeAsInt("SR", 0)).sum());
			vcb.attribute("PE",buffer.stream().flatMap(V->V.ctx.getGenotypes().stream()).mapToInt(G->G.getAttributeAsInt("PE", 0)).sum());

			
			
			final Map<String,Genotype> sample2genotype = new HashMap<>(sampleNames.size());
			
			
			buffer.stream().flatMap(V->V.ctx.getGenotypes().stream()).
				filter(G->isAvailableGenotype(G)).
				forEach(G->{
				sample2genotype.put(G.getSampleName(), G);
			});
			
			vcb.attribute(nSampleInfoHeaderLine.getID(), sample2genotype.size());
			
			for(final String sn: sampleNames)
				{
				if(!sample2genotype.containsKey(sn))
					{
					sample2genotype.put(sn, new GenotypeBuilder(sn,ALLELES_NO_CALLS).
							attribute("SU",0).
							attribute("SR",0).
							attribute("PE",0).
							make());
					}	
				}
			
			vcb.genotypes(sample2genotype.values());
			vcw.add(vcb.make());
			}
		cursor.close();
		vcw.close();vcw=null;
		
		variantsDb1.close();variantsDb1=null;
		environment.close();environment=null;
		return 0;
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(vcfIn);
		CloserUtil.close(vcw);
		CloserUtil.close(variantsDb1);
		CloserUtil.close(environment);
		} 
	}
	 
	public static void main(final String[] args) {
		new LumpySort().instanceMainWithExit(args);
	}
}
