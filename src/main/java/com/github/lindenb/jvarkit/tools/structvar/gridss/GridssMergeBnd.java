/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar.gridss;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

Input is a set of VCF files or a file with the '.list' suffix containing the paths to the vcfs

END_DOC
 */
@Program(name="gridssmergebnd",
	generate_doc=false,
	description="Merge BND results from gridss",
	keywords= {"cnv","indel","sv"},
	creationDate="20200517",
	modificationDate="20200518"
	)
public class GridssMergeBnd extends Launcher{
	private static final Logger LOG = Logger.build(GridssMergeBnd.class).make();
	private static long ID_GENERATOR = 0L;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-d","--distance"},description="within distance. Two bnd are considered the same they're withing that distance.")
	private int withinDistance = 0;
	@Parameter(names={"-a","--attributes"},description="add VCF attributes AC/AN/AF")
	private boolean add_attributes = false;
	@Parameter(names={"-C","--contig"},description="limit to that contig")
	private String onlyContig = null;

	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	@ParametersDelegate
	private WritingSortingCollection sortingDelegate = new WritingSortingCollection();

	private SAMSequenceDictionary theDict = null;
	
	private class BreakPoint implements Locatable {
		final int tid;
		final int start;
		final int end;
		final int sample_id;
		final long id;
		
		BreakPoint(int tid,int start,int end,int sample_id,long id) {
			this.tid=tid;
			this.start=start;
			this.end=end;
			this.sample_id = sample_id;
			this.id = id;
		}
		
		@Override
		public String getContig()
			{
			return theDict.getSequence(this.tid).getSequenceName();
			}
		@Override
		public int getStart()
			{
			return start;
			}
		@Override
		public int getEnd()
			{
			return end;
			}
		
	    public boolean withinDistance(final BreakPoint other, int distance) {
	        return this.tid == other.tid && 
	                CoordMath.overlaps(this.start, this.end, other.start-distance, other.end+distance);
	    	}

		
		public int compare(final BreakPoint o)
			{
			int i= Integer.compare(this.tid, o.tid);
			if(i!=0) return i;
			i= Integer.compare(this.start, o.start);
			if(i!=0) return i;
			return Long.compare(this.id,o.id);
			}
		
		@Override
		public String toString() {
			return getContig()+":"+getStart()+"-"+getEnd()+" "+sample_id;
			}
		}
	
	private class BreakPointCodec extends AbstractDataCodec<BreakPoint> {
		@Override
		public BreakPoint decode(DataInputStream dis) throws IOException
			{
			int tid;
			try {
				tid = dis.readInt();
				}
			catch(EOFException err) {
				return null;
				}
			int start = dis.readInt();
			int end = dis.readInt();
			int sample_id = dis.readInt();
			long id = dis.readLong();
			return new BreakPoint(tid,start,end,sample_id,id);
			}
	
		@Override
		public void encode(DataOutputStream dos, final BreakPoint o)
				throws IOException
			{
			dos.writeInt(o.tid);
			dos.writeInt(o.start);
			dos.writeInt(o.end);
			dos.writeInt(o.sample_id);
			dos.writeLong(o.id);
			}
		
		@Override
		public AbstractDataCodec<BreakPoint> clone()
			{
			return new BreakPointCodec();
			}
		}
	
	@Override
	public int doWork(List<String> args)
		{
		int onlyContigTid = -1;
		try {
			final Map<String,Integer> sample2index = new HashMap<>();			
			final List<String> idx2sample = new ArrayList<>();
			
			SortingCollection<BreakPoint> sorter = SortingCollection.newInstance(
					BreakPoint.class,
					new BreakPointCodec(),(A,B)->A.compare(B),
					this.sortingDelegate.maxRecordsInRam,
					this.sortingDelegate.getTmpPaths()
					);
			
			sorter.setDestructiveIteration(true);
			
			for(final Path path: IOUtils.unrollPaths(args)) {
				LOG.info("Reading "+path);
				
				try(VCFIterator iter= new BcfIteratorBuilder().open(path)) {
					final VCFHeader header= iter.getHeader();
					
					if(header.getNGenotypeSamples()!=1) {
						LOG.error("Expected one sample in "+path);
						return -1;
						}
					final String sn = header.getSampleNamesInOrder().get(0);
					
					if(sample2index.containsKey(sn)) {
						LOG.error("Duplicate sample "+sn+" from "+path);
						return -1;
						}
					final int sample_id = idx2sample.size();
					idx2sample.add(sn);
					sample2index.put(sn, sample_id);
						
					
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
					if(this.theDict==null) {
						this.theDict = dict;
						if(this.onlyContig!=null) {
							onlyContigTid = this.theDict.getSequenceIndex(this.onlyContig);
							if(onlyContigTid<0) throw new JvarkitException.ContigNotFoundInDictionary(this.onlyContig,this.theDict);
							}
						}
					else
						{
						SequenceUtil.assertSequenceDictionariesEqual(dict, this.theDict);
						}
					
					
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						final int tid = this.theDict.getSequenceIndex(ctx.getContig());
						if(tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(ctx.getContig(), this.theDict);
						if(onlyContigTid!=-1 && onlyContigTid!=tid) continue;
						if(onlyContigTid!=-1 && onlyContigTid<tid) break;
						int start = ctx.getStart();
						int end = ctx.getStart();//yes start
						if(ctx.hasAttribute("CIPOS")) {
							final List<Integer> cipos = ctx.getAttributeAsIntList("CIPOS", 0);
							if(cipos.size()>0) start += cipos.get(0);
							if(cipos.size()>1) end += cipos.get(1);
							if(start>end) {
								final int tmp = start;
								start = end;
								end = tmp;
								}
							}
						final QueryInterval qi1 = new QueryInterval(tid, start, end);
						final QueryInterval qintervals[];
						
						if(!ctx.getAttributeAsString(VCFConstants.SVTYPE,"undefined").equals("BND") &&
							ctx.getLengthOnReference()>this.withinDistance) {
							int start2 =  ctx.getEnd();
							int end2 =  ctx.getEnd();
							if(ctx.hasAttribute("CIEND")) {
								final List<Integer> ciend = ctx.getAttributeAsIntList("CIEND", 0);
								if(ciend.size()>0) start2 += ciend.get(0);
								if(ciend.size()>1) end2 += ciend.get(1);
								if(start2>end2) {
									final int tmp = start2;
									start2 = end2;
									end2 = tmp;
									}
								}
							final QueryInterval qi2 = new QueryInterval(tid, start2, end2);
							qintervals = new QueryInterval[]{qi1,qi2};
							}
						else
							{
							qintervals = new QueryInterval[]{qi1};
							}
						for(QueryInterval qi: qintervals) {
							sorter.add(new BreakPoint(qi.referenceIndex, qi.start, qi.end, sample_id,++ID_GENERATOR));
							}
						}//end while
					
					}
				}
			sorter.doneAdding();
			LOG.info("Done adding");
			boolean debug=false;
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.END_KEY,VCFConstants.ALLELE_COUNT_KEY,VCFConstants.ALLELE_FREQUENCY_KEY,VCFConstants.ALLELE_NUMBER_KEY);
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.GENOTYPE_KEY);
			final VCFHeader header = new VCFHeader(metaData, new TreeSet<>(sample2index.keySet()));
			header.setSequenceDictionary(this.theDict);
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			
			
			long count=0L;
			final Allele ref_allele = Allele.create("N", true);
			final Allele alt_allele = Allele.create("<B>", false);
			final List<Allele> alleles = Arrays.asList(ref_allele,alt_allele);
			try(CloseableIterator<BreakPoint> iter0 = sorter.iterator();
				VariantContextWriter vcw = this.writingVariantsDelegate.
						dictionary(this.theDict).
						open(this.outputFile);
				)  {
				vcw.writeHeader(header);
				final PeekIterator<BreakPoint> iter = new PeekIterator<>(iter0);
				final Map<Integer,BreakPoint> sample2bnd = new HashMap<>();

				while(iter.hasNext()) {
					final BreakPoint first = iter.next();
					//debug = first.getStart() >= 32_552_334 && first.getStart()< 50_000_000 && first.getContig().equals("chr2");
					if(debug) System.err.println("A "+first);
					sample2bnd.clear();
					sample2bnd.put(first.sample_id,first);
					while(iter.hasNext()) {
						if(debug) System.err.println("B "+first);
						final BreakPoint bp2 = iter.peek();
						if(debug) System.err.println("C "+bp2);

						if(sample2bnd.values().stream().noneMatch(BP->BP.withinDistance(bp2,this.withinDistance))) {
							if(debug) System.err.println("D "+first+" "+bp2);
							break;
							}
						if(debug) System.err.println("E "+first+" "+bp2+" ");

						sample2bnd.put(bp2.sample_id,iter.next());
						}

					final int end = sample2bnd.values().stream().mapToInt(BP->BP.end).max().orElse(first.end);
					
					if(end-first.getStart() > 1000) continue;

					final VariantContextBuilder vcb = new VariantContextBuilder(null,first.getContig(),first.getStart(),end,alleles);
					if(first.getStart()!=end) vcb.attribute(VCFConstants.END_KEY,end);
					if(this.add_attributes) {
						vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, sample2bnd.size());
						vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY, sample2index.size()*2);
						vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, sample2bnd.size()/(2.0*sample2index.size()));
						}
					vcb.genotypes(
						sample2bnd.keySet().
						stream().
						map(SIDX->new GenotypeBuilder(idx2sample.get(SIDX), alleles).make()).
						collect(Collectors.toList()));
					if(debug) System.err.println("F "+first+" "+count);

					vcw.add(vcb.make());
					count++;
					if(count%10_000==0) {
						LOG.info("N="+count+" last:"+first);
						if(vcw.checkError()) {
							LOG.error("Cannot write!!! "+count+"/"+first);
							return -1;
							}
						}
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			}
		}
	public static void main(final String[] args) {
		new GridssMergeBnd().instanceMainWithExit(args);
		}
	}
