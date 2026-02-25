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
package com.github.lindenb.jvarkit.tools.minigenotyper;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.SampleAndBamFactory;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/*
BEGIN_DOC

## Input

input is a set of indexed bam or a file with the paths ending with '.list'

## Example

```
find src/test/resources/ -type f -name "S*.bam" > jeter.list
java -jar dist/jvarkit.jar minigenotyper \
 	-V src/test/resources/rotavirus_rf.vcf.gz  \
 	-R src/test/resources/rotavirus_rf.fa \
 	jeter.list --skip-illegal-variant > out.vcf
 ```


END_DOC
 */
@Program(name="minigenotyper",
	description="Simple and Stupid Variant Genotyper",
	keywords={"bam","sam","calling","vcf"},
	modificationDate="20260221",
	creationDate="20260221",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class MiniGenotyper extends Launcher   {
	private static final Logger LOG = Logger.of(MiniGenotyper.class);
	private enum Mode {random_access,streaming};
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-V","--variant"},description="Variants to Genotype",required=true)
	private Path variantPath;
	@Parameter(names={"-R","--reference"},description= INDEXED_FASTA_REFERENCE_DESCRIPTION+". This can be used multiple times if there is more than one REF.")
	private List<String> otherReferences = new ArrayList<>();
	@Parameter(names={"-mapq","--mapq"},description="min mapping quality")
	private int mapq  = 1;
	@Parameter(names={"--min-base-quality"},description="min base quality")
	private int min_base_quality = -1;
	@Parameter(names={"--min-gt-depth"},description="min genotype DP")
	private int min_genotype_depth = 5;
	@Parameter(names={"--bad-ad-ratio"},description="Ignore Alt if x < AD ratio < (1-x)")
	private double bad_ad_ratio = 0.1;
	@Parameter(names={"--disable-overlap-detection"},description="Disable Paired-end read overlap detection")
	private boolean disable_overlap_detection = false;
	@Parameter(names={"--skip-illegal-variant"},description="just skip illegal variant in the VCF (not diallelic SNPs)")
	private boolean just_skip_illegal_variant = false;
	@Parameter(names={"--mode"},description="how to scan the bam. Using random-access (ok if small number of variant) or streaming (no index required, large number of variants)")
	private Mode run_mode = Mode.random_access;
	@Parameter(names={"--no-id"},description="don't print ID column")
	private boolean without_ID = false;
	@Parameter(names={"--no-info"},description="don't print INFO column")
	private boolean without_INFO = false;

	
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

    
    private static class Call implements Cloneable
	    {
    	int tid;
    	int position;
    	byte ref;
    	byte alt;
    	int sample_index;
    	int n_REF=0;
		int n_ALT=0;
		int n_OTHER=0;
		
    	public int compare1(final Call o) {
    		int i = Integer.compare(this.tid, o.tid);
    		if(i!=0) return i;
			i = Integer.compare(this.position, o.position);
			if(i!=0) return i;
			i = Byte.compare(this.ref, o.ref);
			if(i!=0) return i;
			i = Byte.compare(this.alt, o.alt);
	    	return i;
    		}
    	public int compare2(final Call o) {
    		int i = compare1(o);
    		if(i!=0) return i;
			i = Integer.compare(this.sample_index, o.sample_index);
	    	return i;
    		}
    	public int getDepth() {
    		return n_REF+n_ALT+n_OTHER; 
    		}
		public int compareOnPosition(final Call other) {
			int i = Integer.compare(this.tid, other.tid);
			if(i!=0) return i;
			i = Integer.compare(this.position, other.position);
			if(i!=0) return i;
			i = Byte.compare(this.ref, other.ref);
			if(i!=0) return i;
			i = Byte.compare(this.alt, other.alt);
			return i;
			}

    	@Override
    	protected Call clone()   {
    		final Call c = new Call();
    		c.tid = this.tid;
    		c.position = this.position;
    		c.ref = this.ref;
    		c.alt = this.alt;
    		c.sample_index = this.sample_index;
    		c.n_REF = this.n_REF;
    		c.n_ALT = this.n_ALT;
    		c.n_OTHER = this.n_OTHER;
    		return c;
    		}
    	@Override
    	public String toString()
    		{
    		return "sample_index:"+sample_index+" position:"+position+" ref:"+ref+" alt:"+alt;
    		}
	    }
    
    private static class CallCodec extends AbstractDataCodec<Call> {
		@Override
		public void encode(DataOutputStream dos, final Call c) throws IOException
			{
			dos.writeInt(c.tid);
			dos.writeInt(c.position);
			dos.writeByte(c.ref);
			dos.writeByte(c.alt);
			dos.writeInt(c.sample_index);
			dos.writeInt(c.n_REF);
			dos.writeInt(c.n_ALT);
			dos.writeInt(c.n_OTHER);
			}
		@Override
		public Call decode(DataInputStream dis) throws IOException
			{
			final Call c = new Call();
			try {
				c.tid = dis.readInt();
				}
			catch(EOFException err) {
				return null;
				}
			c.position = dis.readInt();
			c.ref = dis.readByte();
			c.alt = dis.readByte();
			c.sample_index  = dis.readInt();
			c.n_REF = dis.readInt();
			c.n_ALT = dis.readInt();
			c.n_OTHER = dis.readInt();
			return c;
			}
		@Override
		public AbstractDataCodec<Call> clone()
			{
			return new CallCodec();
			}
	    }


    
    private boolean checkValidVariant(final VariantContext ctx) {
    	if(ctx.getNAlleles()!=2) {
    		final String msg ="only diallelic variants enabled: "+ctx;
    		if(just_skip_illegal_variant) {
    			LOG.warn(msg);
    			return false;
    			}
    				
    		throw new IllegalArgumentException(msg);
    		}
    	if(!ctx.isSNP()) {
    		final String msg = "only SNP variants enabled: "+ctx;
    		if(just_skip_illegal_variant) {
    			LOG.warn(msg);
    			return false;
    			}
    				
    		throw new IllegalArgumentException(msg);
    		}
    	for(Allele a: ctx.getAlleles()) {
    		if(a.length()!=1)  {
    			final String msg = "only Allele length=1 enabled: "+ctx;
        		if(just_skip_illegal_variant) {
        			LOG.warn(msg);
        			return false;
        			}
        				
        		throw new IllegalArgumentException(msg);
    			}
    		if(!AcidNucleics.isATGC(a)) {
    			final String msg = "only Allele [ATGC] enabled: "+ctx;
        		if(just_skip_illegal_variant) {
        			LOG.warn(msg);
        			return false;
        			}        				
        		throw new IllegalArgumentException(msg);
    			}
    		}
    	return true;
    	}
    /** get mate end position or mate-start */
    private static int getMateEnd(final SAMRecord rec0) {
        return rec0.hasAttribute(SAMTag.MC)?SAMUtils.getMateAlignmentEnd(rec0):rec0.getMateAlignmentStart();
        }
    
    private boolean scanRead(final SAMRecord rec, final Call call) {
    	boolean got_position = false;
		if(call.position< rec.getStart() ||  rec.getEnd()<call.position) {
			return false;
			}
		if(!disable_overlap_detection && 
				rec.getReadPairedFlag() && 
				!rec.getMateUnmappedFlag() && 
				rec.getReferenceName().equals(rec.getMateReferenceName()) && 
				rec.getAlignmentStart() <= rec.getMateAlignmentStart() && 
				rec.getMateAlignmentStart() <= call.position &&
				getMateEnd(rec) >= call.position
				) {
				return false;
				}
		
        final Cigar cigar= rec.getCigar();
        if(cigar==null) return false;
        final byte[] bases = rec.getReadBases();
        final byte[] quals = rec.getBaseQualities();
        if(bases ==  SAMRecord.NULL_SEQUENCE) {
        	return false;
        	}
        if(quals !=  SAMRecord.NULL_QUALS && bases.length != bases.length) {
        	return false;
        	}
        for(AlignmentBlock block: rec.getAlignmentBlocks()) {
        	if(block.getReferenceStart() > call.position) {
        		break;
        		}
        	if(block.getReferenceStart()+block.getLength()-1 < call.position) {
        		continue;
        		}
        	for(int i=0;i< block.getLength();i++) {
        		final int pos1 = block.getReferenceStart() +i;
				if(pos1 < call.position ) continue;
				if(pos1 > call.position ) break;
				final int readpos0 = block.getReadStart()-1+i;
				
				if(quals !=  SAMRecord.NULL_QUALS && quals[readpos0] < this.min_base_quality) {
					break;
					}
				got_position = true;
				final byte readBase = (byte)Character.toUpperCase(  bases[readpos0] );
				if( readBase==call.ref) {
					call.n_REF++;
					}
				else if(readBase==call.alt) {
					call.n_ALT++;
					}
				else
					{
					call.n_OTHER++;
					}
				break;		                            				
        		}
        	}
        return got_position;
    	}
    
    private void streaming(
	    final SamReader sr,
		final int sample_index,
		final SAMSequenceDictionary variantDict,
		final List<Call> variants0,
		SortingCollection<Call> sorter,
		final Consumer<String> logger
		) throws IOException  {
    	final long startMillisec = System.currentTimeMillis();
    	for(Call c: variants0) {
    		c.sample_index=sample_index;
    		c.n_ALT=0;
    		c.n_REF=0;
    		c.n_OTHER=0;
    		}
		int prev_tid=-1;
		// call for the current chromosome
		final List<Call> variants = new ArrayList<MiniGenotyper.Call>(variants0.size());
		final Function<String,String> contigConverter = ContigNameConverter.fromOneDictionary(variantDict);
		int call_index=-1;
		try(CloseableIterator<SAMRecord> iter2= sr.iterator()) {
			while(iter2.hasNext()) {
				final SAMRecord rec = iter2.next();
				if(!SAMRecordDefaultFilter.accept(rec, this.mapq)) {
					continue;
					}
				final String ctg= contigConverter.apply(rec.getContig());
				if(StringUtils.isBlank(ctg)) {
					continue;
					}
				final int tid = variantDict.getSequenceIndex(rec.getContig());
				if(tid<0) {
					continue;
					}
				if(prev_tid==-1 || prev_tid!=tid) {
					variants.clear();
					for(Call c: variants0) {
						if(c.tid == tid  && c.position>=rec.getAlignmentStart()) {
							variants.add(c);
							}
						}
					Collections.sort(variants,(A,B)->Integer.compare(A.position,B.position));
					if(!variants.isEmpty()) {
						logger.accept(ctg+" "+StringUtils.niceDuration(System.currentTimeMillis()- startMillisec));
						}
					prev_tid=tid;
					call_index=0;
					}
				while(call_index < variants.size() && variants.get(call_index).position < rec.getAlignmentStart()) {
					call_index++;
					}
				for(int i=call_index ; i< variants.size() ;++i) {
					final Call c = variants.get(i);
					if(c.position> rec.getAlignmentEnd()) break;
					scanRead(rec, c);
					}
				}
			}
		for(Call c: variants0) {
			sorter.add(c.clone());
			}
		}
    
    private void randomAccess(
    		final SamReader sr,
    		final int sample_index,
    		final SAMSequenceDictionary variantDict,
    		final SortingCollection<Call> sorter,
    		final Consumer<Long> logger
    		) throws IOException {

		final SAMFileHeader samHeader  = sr.getFileHeader();
		final Function<String,String> contigConverter = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(samHeader));
		long n_processed=0;
		try(CloseableIterator<VariantContext> iter= new VCFIteratorBuilder().open(this.variantPath)) {
			while(iter.hasNext()) {
				final VariantContext ctx = iter.next();
				if(!checkValidVariant(ctx)) continue;
				final String contig = contigConverter.apply(ctx.getContig());
				final Call call = new Call();
				call.tid = variantDict.getSequenceIndex(ctx.getContig());
				if(call.tid<0) continue;
				call.ref  = (byte)Character.toUpperCase( ctx.getReference().getBases()[0]);
				call.alt =  (byte)Character.toUpperCase( ctx.getAlternateAllele(0).getBases()[0]);
				call.sample_index = sample_index;
				call.position = ctx.getStart();
				if(StringUtils.isBlank(contig)) {
					//nothing
					}
				else
					{
					try(CloseableIterator<SAMRecord> iter2= sr.query(contig, call.position, call.position ,false)) {
						while(iter2.hasNext()) {
							final SAMRecord rec = iter2.next();
							if(!SAMRecordDefaultFilter.accept(rec, this.mapq)) continue;
        					scanRead(rec,call);
        					} /* end loop SAM Record */
						} /* end contig found in SAM */
					}
				sorter.add(call);
				n_processed++;
				if(n_processed%1000L==0) {
					logger.accept(n_processed);
					}
				} /* while variant has Next */
			} /* open VCF */
    	}
    
    @Override
    public int doWork(final List<String> args) {

    	try {
    		final List<Call> variantList = this.run_mode.equals(Mode.streaming)?new ArrayList<>():null;	
			final SampleAndBamFactory sampleAndBamFacory = new SampleAndBamFactory();
			sampleAndBamFacory.setEnableNonIndexedBam(this.run_mode.equals(Mode.streaming));
			//sampleAndBamFacory.setEnableWithoutReference(true);
			
			for(final String fna: this.otherReferences) {
				final Path fa = Paths.get(fna);
				IOUtil.assertFileIsReadable(fa);
				sampleAndBamFacory.reference(fa);
				}
			
			final List<SampleAndBamFactory.SampleAndBamRecord> bams = sampleAndBamFacory.parse(IOUtils.unrollPaths(args));
			
    		if(bams.isEmpty()) {
    			LOG.error("Bam path(s) missing");
    			return -1;
    			}
    		
    		final SAMSequenceDictionary variantDict ;
    		int count_variant = 0; 
    		try(VCFIterator iter= new VCFIteratorBuilder().open(this.variantPath)) {
    			variantDict =  SequenceDictionaryUtils.extractRequired(iter.getHeader());
    			while(iter.hasNext()) {
    				final VariantContext ctx = iter.next();
    				if(!checkValidVariant(ctx)) continue;
    				if(variantList!=null) {
    					final Call call = new Call();
        				call.tid = variantDict.getSequenceIndex(ctx.getContig());
        				if(call.tid<0) continue;
        				call.ref  = (byte)Character.toUpperCase( ctx.getReference().getBases()[0]);
    					call.alt =  (byte)Character.toUpperCase( ctx.getAlternateAllele(0).getBases()[0]);
    					call.position = ctx.getStart();
    					variantList.add(call);
    					}
    				count_variant++;
    				}
    			}
    		if(count_variant==0) {
    			LOG.error("VCF is empty : "+this.variantPath);
    			return -1;
    			}
    		if(variantList!=null) {
    			Collections.sort(variantList, (A,B)->A.compareOnPosition(B));
    			}
    		
    		
    		final SortingCollection<Call> sorter = SortingCollection.newInstance(Call.class, new CallCodec(),
    				(A,B)->A.compare2(B),
    				writingSortingCollection.getMaxRecordsInRam(),
    				writingSortingCollection.getTmpPaths()
    				);	
    		sorter.setDestructiveIteration(true);
    		
	        for(int bam_idx=0; bam_idx < bams.size();++bam_idx) {
	        		final int final_bam_idx = bam_idx;
	        		final int final_count_variant = count_variant;
	        		final SampleAndBamFactory.SampleAndBamRecord bamRecord = bams.get(bam_idx);
	        		LOG.info("Scanning "+bamRecord.getPath()+ "("+(bam_idx+1)+"/"+bams.size()+")");
	        		final SamReaderFactory srf=createSamReaderFactory().referenceSequence(bamRecord.getReference().orElse(null));
	        		try(SamReader sr = srf.open(bamRecord.getPath())) {
		        		if(this.run_mode.equals(Mode.random_access)) {
			        			randomAccess(sr, bam_idx, variantDict, sorter, N->{
			        				LOG.info(""+N+"/"+final_count_variant+" ("+(N/(double)final_count_variant)+") bam:"+"("+(final_bam_idx+1)+"/"+bams.size()+")");
			        				});
		        			}
		        		else
		        			{
		        			streaming(sr, bam_idx, variantDict, variantList, sorter,CTG->{
		        				LOG.info(CTG+" bam:"+"("+(final_bam_idx+1)+"/"+bams.size()+")");
		        				});	
		        			}
	        			}
        			} /* end loop over each bam */
    			sorter.doneAdding();
        	    /*  VCF metadata */
        	    final Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();	        		
	        		
		        try(CloseableIterator<Call> iter0 = sorter.iterator()) {
		        	try(EqualIterator<Call> iter  = new EqualIterator<MiniGenotyper.Call>(iter0, (A,B)->A.compare1(B))) {
		            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
		            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
		            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
		            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
		            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY));
		            metaData.add(new VCFFormatHeaderLine("N", 1, VCFHeaderLineType.Integer, "Number of other bases"));

		            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));
		            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
		            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
		            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
		            
			        final VCFHeader header =new VCFHeader(
		                    metaData, 
		                    bams.stream().map(REC->REC.getSampleName()).collect(Collectors.toList())
		                    );
			        header.setSequenceDictionary(variantDict);
			        JVarkitVersion.getInstance().addMetaData(this, header);
			        
		            try(VariantContextWriter vcw = this.writingVariantsDelegate.dictionary(variantDict).open(this.outputFile)) {
		            	vcw.writeHeader(header);
		            	while(iter.hasNext()) {
		            			final List<Call> row = iter.next();
		            			final Call first = row.get(0);
		            			final Allele REF=  Allele.create(first.ref, true);
		            			final Allele ALT=  Allele.create(first.alt, false);
		            			final VariantContextBuilder vcb = new VariantContextBuilder("jvarkit", variantDict.getSequence(first.tid).getContig(), first.position, first.position, Arrays.asList(REF,ALT));
		            			final Map<String,Genotype> sample2genotypes = new HashMap<>(bams.size());
		            			int ac=0;
		            			int an=0;
		            			int dp= 0;
		            			for(Call c:row) {
		            				String filter=null;
		            				final String sample = bams.get(c.sample_index).getSampleName();
		            				dp+= c.n_ALT+c.n_REF;
		            				final GenotypeBuilder gb=new GenotypeBuilder(sample);
		            				gb.DP(c.getDepth());
		            				if((c.n_REF+c.n_ALT)< Math.max( this.min_genotype_depth ,0) ) {
		            					gb.alleles(Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
		            					gb.filter("LowQual");
		            					gb.GQ(0);
		            					}
		            				else
		            					{
		            					gb.AD(new int[] {c.n_REF,c.n_ALT});
		            					double qual = Math.min(c.getDepth(),30)/30.0;
		            					
		            					if(c.n_REF>0 && c.n_ALT==0) {
		            						gb.alleles(Arrays.asList(REF,REF));
		            						qual *= (c.n_REF/(double)c.getDepth());
		            						an+=2;
		            						if(c.n_OTHER >= c.n_REF*this.bad_ad_ratio) filter="LowQual";
		            						}
		            					else if(c.n_REF==0 && c.n_ALT>=0) {
		            						gb.alleles(Arrays.asList(ALT,ALT));
		            						qual *= (c.n_ALT/(double)c.getDepth());
		            						if(c.n_OTHER >= c.n_ALT*this.bad_ad_ratio) filter="LowQual";
		            						an+=2;
		            						ac+=2;
		            						}
		            					else
		            						{
		            						final double ratio = c.n_ALT/(double)(c.n_REF+c.n_ALT);
		            						if(ratio <= this.bad_ad_ratio) {
		            							gb.alleles(Arrays.asList(REF,REF));
		            							qual *= (c.n_REF/(double)c.getDepth());
		            							if(c.n_OTHER >= c.n_REF*this.bad_ad_ratio ) filter="LowQual";
		            							an+=2;
		            							}
		            						else if(ratio >= (1.0-this.bad_ad_ratio)) {
		            							gb.alleles(Arrays.asList(ALT,ALT));
		            							qual *= (c.n_ALT/(double)c.getDepth());
		            							an+=2;
		            							ac+=2;
		            							if(c.n_OTHER >= c.n_ALT*this.bad_ad_ratio ) filter="LowQual";
		            							}
		            						else
		            							{
		            							gb.alleles(Arrays.asList(REF,ALT));
		            							qual *= Math.abs(0.5-((c.n_ALT+c.n_REF)/(double)c.getDepth()))/0.5;
		            							if( c.n_OTHER>=c.n_ALT || c.n_OTHER>=c.n_REF) filter="LowQual";
		            							an+=2;
		            							ac++;
		            							}
		            						
		            						}
		            					gb.attribute("N", c.n_OTHER);
	            						gb.GQ((int)(qual*99.0));
		            					if(filter!=null) gb.filter(filter);
		            					
		            					}
		            				sample2genotypes.put(sample,gb.make());
		            				}
		            			if(!without_ID) {
			            			vcb.id(variantDict.getSequence(first.tid).getContig()+":"+first.position+":"+(char)first.ref+":"+(char)first.alt);
			            			}
		            			vcb.genotypes(new ArrayList<>(sample2genotypes.values()));
		            			
		            			if(!without_INFO) {
			            			vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY, an);
			            			vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, ac);
			            			vcb.attribute(VCFConstants.DEPTH_KEY, dp);
			            			if(an>0) vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, ac/(double)an);
			            			}
		            			
		            			vcw.add(vcb.make());
		            			}
		            		} /* end writer */
		        		} /* end equal iterator */
		        	} /* end sort iterator */
		        sorter.cleanup();
		        return 0;
	    	} catch(final Throwable err) {
	    		LOG.error(err);
	    		return -1;
	    	}
	    }

    public static void main(String[] args)
        {
        new MiniGenotyper().instanceMainWithExit(args);
        }

    }
