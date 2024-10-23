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
package com.github.lindenb.jvarkit.tools.calling;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFSampleHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/*
BEGIN_DOC


## Cited-In

  * "Direct Head-to-Head Evaluation of Recombinant Adeno-associated Viral Vectors Manufactured in Human versus Insect Cells". Kondratov & al. Molecular Therapy. [https://doi.org/10.1016/j.ymthe.2017.08.003](https://doi.org/10.1016/j.ymthe.2017.08.003).
  *  METHODS OF ENHANCING BIOLOGICAL POTENCY OF BACULOVIRUS SYSTEM-PRODUCED RECOMBINANT ADENO-ASSOCIATED VIRUS   United States Patent Application 20200123572 . United States Patent Application 20200123572

## Example

```bash
$  java -jar dist/minicaller.jar -R ref.fa  bam.list > out.vcf
```

```
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Depth ReforAlt|Strand : RF,RR,AF,AR">
##FORMAT=<ID=DPG,Number=G,Type=Integer,Description="Depth for each allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Variant is indel">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	4	.	T	A	.	.	DP=65	GT:DP:DP4:DPG	0/1:22:20,0,2,0:20,2	./.	./.	./.
rotavirus	5	.	T	A	.	.	DP=83	GT:DP:DP4:DPG	0:27:27,0,0,0:27	./.	0/1:25:20,0,5,0:20,5	./.
rotavirus	6	.	T	A	.	.	DP=97	GT:DP:DP4:DPG	0:30:30,0,0,0:30	0/1:21:20,0,1,0:20,1	0/1:33:31,0,2,0:31,2	./.
rotavirus	7	.	T	A	.	.	DP=112	GT:DP:DP4:DPG	0/1:38:36,0,2,0:36,2	0/1:23:21,0,2,0:21,2	0:37:37,0,0,0:37	./.
rotavirus	8	.	A	C	.	.	DP=122	GT:DP:DP4:DPG	0/1:41:38,0,3,0:38,3	0/1:26:25,0,1,0:25,1	0/1:40:38,0,2,0:38,2	./.
rotavirus	9	.	A	C	.	.	DP=139	GT:DP:DP4:DPG	0/1:46:44,0,2,0:44,2	0/1:29:27,0,2,0:27,2	0/1:48:44,0,4,0:44,4	./.
```

END_DOC
 */
@Program(name="minicaller",
	description="Simple and Stupid Variant Caller designed for @AdrienLeger2",
	keywords={"bam","sam","calling","vcf"},
	modificationDate="20220705",
	creationDate="201500306",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class MiniCaller extends Launcher   {
	private static final Logger LOG = Logger.build(MiniCaller.class).make();
	private static final String CLIP_FLAG = "!";
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-d","--mindepth"},description="Min depth")
	private int min_depth = 20 ;
	@Parameter(names={"-R","--reference"},description="Main fasta reference. "+INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxPath = null;
	@Parameter(names={"--other-reference"},description="Other fasta references if you mix bam mapped on different fasta (will try to convert chromosomes names). " + INDEXED_FASTA_REFERENCE_DESCRIPTION,required=false)
	private List<String> otherReferences = new ArrayList<>();
	@Parameter(names={"-r","--region"},description=IntervalParser.OPT_DESC,required=true)
	private String rgnStr  = null;
	@Parameter(names={"-mapq","--mapq"},description="min mapping quality")
	private int mapq  = 1;
	@Parameter(names={"--min-base-quality"},description="min base quality")
	private int min_base_quality = 1;
	@Parameter(names={"--min-gt-depth"},description="min genotype DP")
	private int min_genotype_depth = 1;
	@Parameter(names={"--min-gt-allele-depth"},description="min genotype allele DP")
	private int min_genotype_allele_depth = 1;
	@Parameter(names={"--gt-fraction"},description="ignore genotype ALT/(REF+ALT) < x")
	private double min_genotype_fraction = 1.0/20.0;
	@Parameter(names={"--bad-ad-ratio"},description="Filter Genotype if x< ALT/(REF+ALT) < (1-x).")
	private double bad_ad_ratio = 0.2;
	@DynamicParameter(names = "-D", description = "extra parameters. Undocumented.",hidden=true)
	private Map<String, String> __dynaParams = new HashMap<>();


	private static final int PLOIDY = 2;
	private final AttributeMap attMap = AttributeMap.verbose(AttributeMap.wrap(__dynaParams),(K)->{
		LOG.info("Using default parameter for '"+K+"'.");
		});
	
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();


    private static class AlleleInfo {
    	final Allele allele;
    	int count_F=0;
    	int count_R=0;
    	long mapq_sum = 0;
    	boolean masked = false;
    	AlleleInfo(final Allele allele) {
    		this.allele = allele;
    		}
    	int count() {
    		return count_F+count_R;
    		}
    	}
    
    private static class SampleInfo {
    	final String sn;
    	final Allele refAllele;
    	final Map<Allele,AlleleInfo> allele2info = new TreeMap<>();
    	String gtFilter = null;
    	int clip_sum = 0;
    	SampleInfo(String sn,final Allele refAllele) {
    		this.sn = sn;
    		this.refAllele = refAllele;
    		}
    	void visit(final Call c) {
    		if(c.alt.equals(CLIP_FLAG)) {
    			clip_sum++;
    			return;
    			}
    		final Allele alt = createAllele(c.alt);
    		AlleleInfo ai = allele2info.get(alt);
    		if(ai==null) {
    			ai = new AlleleInfo(alt);
    			allele2info.put(alt, ai);
    			}
    		if(c.negativeStrand) {
    			ai.count_R++;
    			}
    		else
    			{
    			ai.count_F++;
    			}
    		ai.mapq_sum  += c.mq;
    		}
    	private Allele createAllele(final String a) {
    		final Allele alt = Allele.create(a,false);
    		if(alt.equals(this.refAllele, true)) return this.refAllele;
    		return alt;
    		}
    	
    	}
    
    private static class Call
	    {
    	int sample_index;
    	int position;
    	String ref;
    	String alt;
    	short mq;
    	boolean negativeStrand;
    	public int compare1(final Call o) {
			int i = Integer.compare(this.position, o.position);
			if(i!=0) return i;
			i = Integer.compare(this.sample_index, o.sample_index);
			if(i!=0) return i;
			i = ref.compareTo(o.ref);
	    	return i;
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
			dos.writeInt(c.sample_index);
			dos.writeInt(c.position);
			dos.writeUTF(c.ref);
			dos.writeUTF(c.alt);
			dos.writeShort(c.mq);
			dos.writeBoolean(c.negativeStrand);
			}
		@Override
		public Call decode(DataInputStream dis) throws IOException
			{
			final Call c = new Call();
			try {
				c.sample_index = dis.readInt();
				}
			catch(EOFException err) {
				return null;
				}
			c.position = dis.readInt();
			c.ref = dis.readUTF();
			c.alt = dis.readUTF(); if(c.alt.equals(CLIP_FLAG)) c.alt= CLIP_FLAG;
			c.mq = dis.readShort();
			c.negativeStrand = dis.readBoolean();
			return c;
			}
		@Override
		public AbstractDataCodec<Call> clone()
			{
			return new CallCodec();
			}
	    }


    
    private static class SAMSequenceDictionaryPath {
    	final Path path;
    	final SAMSequenceDictionary dict;
    	SAMSequenceDictionaryPath(final Path path,final SAMSequenceDictionary dict) {
    		this.path = path;
    		this.dict = dict;
    		}
    	}
    
    @Override
    public int doWork(final List<String> args) {

    	try {
			final List<Path> bamPaths = IOUtils.unrollPaths(args);
			final List<String> samplesList = new ArrayList<>();
    		if(bamPaths.isEmpty()) {
    			LOG.error("Bam path(s) missing");
    			return -1;
    			}
    		
    		final List<SAMSequenceDictionaryPath> all_dictionaries = new ArrayList<>();
    		for(Path otherFaidx: IOUtils.unrollPaths(this.otherReferences)) {
    			all_dictionaries.add(new SAMSequenceDictionaryPath(otherFaidx,SequenceDictionaryUtils.extractRequired(otherFaidx)));
    			}
    	    /*  VCF metadata */
    	    final Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();

        	final SAMSequenceDictionary dictionary = SequenceDictionaryUtils.extractRequired(this.faidxPath);
			all_dictionaries.add(new SAMSequenceDictionaryPath(this.faidxPath,dictionary));

			final SimpleInterval interval = new IntervalParser(dictionary).
					apply(this.rgnStr).
					orElseThrow(IntervalParser.exception("Cannot parse region "+this.rgnStr));
    			
	    		final SortingCollection<Call> sorter = SortingCollection.newInstance(Call.class, new CallCodec(),
	    				(A,B)->A.compare1(B),
	    				writingSortingCollection.getMaxRecordsInRam(),
	    				writingSortingCollection.getTmpPaths()
	    				);	
	    		sorter.setDestructiveIteration(true);
	        	for(final Path bamPath: bamPaths) {
	        		LOG.info("Scanning "+bamPath);
	        		final SAMSequenceDictionary bamDict = SequenceDictionaryUtils.extractRequired(bamPath);
		        	final ContigNameConverter bamCtgConverter = ContigNameConverter.fromOneDictionary(bamDict);
	
	        		final Path bamFastaPath = all_dictionaries.stream().
	        				filter(SDP->SequenceUtil.areSequenceDictionariesEqual(bamDict, SDP.dict)).
	        				map(SDP->SDP.path).
	        				findFirst().orElseThrow(()->new FileNotFoundException("Cannot find fasta Reference for "+bamPath));
	            	try(ReferenceSequenceFile bamReferenceSequenceFile= ReferenceSequenceFileFactory.getReferenceSequenceFile(bamFastaPath)) {
	            		final SamReaderFactory srf=createSamReaderFactory().referenceSequence(bamFastaPath);
		        		try(SamReader sr = srf.open(bamPath)) {
		        			final SAMFileHeader header = sr.getFileHeader();
		        			String sn = header.getReadGroups().stream().
		        				map(RG->RG.getSample()).
		        				filter(SN->!StringUtils.isBlank(SN)).
		        				findFirst().
		        				orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamPath));
		        			String sn2 = sn;
		        			int tmp_sn_index =-1;
		        			int try_count=1;
		        			for(;;) {
		        				tmp_sn_index = samplesList.indexOf(sn2);
		        				if( tmp_sn_index < 0 ) {
		        					tmp_sn_index = samplesList.size();
		        					samplesList.add(sn2);
		        					break;
		        					}
		        				try_count++;
		        				sn2 = sn+"_"+try_count;
		        				}
		        			final int sample_index = tmp_sn_index;
		        			metaData.add(new VCFSampleHeaderLine("<ID="+sn2+",path="+bamPath+">", VCFHeaderVersion.VCF4_3));
		        			
		        			final String bamCtg = bamCtgConverter.apply(interval.getContig());
		        			if(StringUtils.isBlank(bamCtg)) {
		        				LOG.warn("Cannot map contig "+interval+" for reference "+bamFastaPath);
		        				}
		        			else
		        				{
		        				final GenomicSequence genomicSequence = new GenomicSequence(bamReferenceSequenceFile, bamCtg);
		        				final int minSeqLength = this.attMap.getIntAttribute("min.read.length").orElse(1);
		        				final SimpleInterval bamInterval = interval.renameContig(bamCtg);
			        			try(CloseableIterator<SAMRecord> iter = sr.query(bamInterval.getContig(),bamInterval.getStart(),bamInterval.getEnd(), false)) {
			        				while(iter.hasNext()) {
			        					final SAMRecord rec = iter.next();
			        					if(!SAMRecordDefaultFilter.accept(rec, this.mapq)) continue;
			        					if(rec.getReadLength() < minSeqLength) continue;
			                            final Cigar cigar= rec.getCigar();
			                            if(cigar==null) continue;
			                            final byte[] bases = rec.getReadBases();
			                            final byte[] quals = rec.getBaseQualities();
			                            if(bases ==  SAMRecord.NULL_SEQUENCE) continue;
			                            if(quals ==  SAMRecord.NULL_QUALS) continue;
			                            if(bases.length != bases.length) continue;
			                            int refpos1 = rec.getUnclippedStart();
			                            int readpos0 = 0;
			                            for(CigarElement ce : cigar) {
			                            	final CigarOperator op = ce.getOperator();
			                            	final int len = ce.getLength();

			                            	try {
			                            	
				                            	switch(op) {
				                            		case H: case S:
				                            			{
				                            			for(int i=0;i<len;i++) {
				                            				int pos = refpos1 +i;
				                            				if(pos< interval.getStart()) continue;
				                            				if(pos> interval.getEnd()) break;
				                            				if(pos <1) continue;
				                            				if(pos-1 >= genomicSequence.length() ) break;
	
				                            				final Call call = new Call();
				                            				call.sample_index = sample_index;
				                            				call.position = pos;
				                            				call.ref = String.valueOf(Character.toUpperCase(genomicSequence.charAt(pos-1)));
				                            				call.alt = CLIP_FLAG;
				                            				call.mq = (short)rec.getMappingQuality();
				                            				call.negativeStrand = rec.getReadNegativeStrandFlag();
				                            				sorter.add(call);
				                            				}
				                            			if(op.equals(CigarOperator.S)) {
				                            				readpos0+=len;
				                            				}
				                            			refpos1+=len;
				                            			break;
				                            			}
				                            		case P: break;
				                            		case D: case N: {
				                            			if(CoordMath.overlaps(interval.getStart(), interval.getEnd(), refpos1, refpos1)) {
						                            		final Call call = new Call();
				                            				call.sample_index = sample_index;
				                            				call.position = refpos1-1;// one base before
				                            				final String prevBase = String.valueOf(genomicSequence.charAt(call.position-1)).toUpperCase();
				                            				call.ref = genomicSequence.subSequence(call.position-1,call.position-1+len).toString().toUpperCase();
				                            				call.alt = prevBase;
				                            				call.mq = (short)rec.getMappingQuality();
				                            				call.negativeStrand = rec.getReadNegativeStrandFlag();
				                            				sorter.add(call);
					                            			}
				                            			refpos1 += len;
				                            			break;
				                            			}
				                            		case I: {
				                            			if(CoordMath.overlaps(interval.getStart(), interval.getEnd(), refpos1, refpos1+len-1)) {
						                            		final Call call = new Call();
				                            				call.sample_index = sample_index;
				                            				call.position = refpos1-1;// one base before
				                            				final String prevBase = String.valueOf(genomicSequence.charAt(call.position-1)).toUpperCase();
				                            				call.ref = prevBase;
				                            				call.alt = prevBase + new String(bases,readpos0,len);
				                            				call.mq = (short)rec.getMappingQuality();
				                            				call.negativeStrand = rec.getReadNegativeStrandFlag();
				                            				sorter.add(call);
					                            			}
				                            			readpos0+= len;
				                            			break;
				                            			}
				                            		case X: case M: case EQ: {
				                            			for(int i=0;i<len;i++) {
				                            				int pos = refpos1 +i;
				                            				if(pos< interval.getStart()) continue;
				                            				if(pos> interval.getEnd()) break;
				                            				if(pos-1 >= genomicSequence.length() ) break;
				                            				if(quals[readpos0] < this.min_base_quality) continue;
				                            				char refB = Character.toUpperCase(genomicSequence.charAt(pos-1));
				                            				byte readB = bases[readpos0+i];
				                            				
				                            				final Call call = new Call();
				                            				call.sample_index = sample_index;
				                            				call.position = pos;
				                            				call.ref = String.valueOf(refB);
				                            				call.alt = String.valueOf((char)readB);
				                            				call.mq = (short)rec.getMappingQuality();
				                            				call.negativeStrand = rec.getReadNegativeStrandFlag();
				                            				sorter.add(call);
				                            				}
			                            				refpos1 += len;
			                            				readpos0 += len;
				                            			break;
				                            			}
				                            		default:throw new IllegalStateException();
				                            		}
			                            		}
			                            	catch(final Throwable err) {
			                            		err.printStackTrace();
			                            		throw new RuntimeException("Error with "+rec.getSAMString()+ " "+op.name()+len,err);
			                            		}
			                            	}
			                            
			        					}//end while iter.hasNext
			        				}/* end of iterator */
			        			} // end if contig is not blank
		        			} /* end of open current SAM */
		            	}/* end of REF sequence */
	        		}/* end of loop over BAM */
	        	
	        sorter.doneAdding();
	        try(CloseableIterator<Call> iter = sorter.iterator()) {
	            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
	            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY));
	            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
	            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
	            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
	            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));
	            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
	            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
	            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));

		        final VCFHeader header =new VCFHeader(
	                    metaData , samplesList
	                    );
		        header.setSequenceDictionary(dictionary);
	            try(VariantContextWriter vcw = this.writingVariantsDelegate.dictionary(dictionary).open(this.outputFile)) {
	            	vcw.writeHeader(header);
		        	int prev_pos = -1;
		        	String prev_ref = null;
		        	final Map<Integer,SampleInfo> sampleid2info = new HashMap<>(samplesList.size());
	        		for(;;) {
	        			final Call c = iter.hasNext()?iter.next():null;
	        			
	        			if(c==null || !(prev_pos==c.position && c.ref.equals(prev_ref))) {
		        			if(!sampleid2info.isEmpty()) {
		        				for(SampleInfo si: sampleid2info.values()) {
		        					final List<AlleleInfo> ailist = si.allele2info.values().
	        						stream().
		        		    			sorted((A,B)->Integer.compare(B.count(), A.count())).
		        		    			limit(PLOIDY).
		        		    			sorted((A,B)->A.allele.compareTo(B.allele)).
		        		    			collect(Collectors.toCollection(ArrayList::new));
	        		    		
		        					if(ailist.stream().mapToInt(AI->AI.count()).sum() <this.min_genotype_depth) {
			        					for(AlleleInfo ai: ailist) {
			        						ai.masked = true;
			        						}
			        					si.gtFilter = "LowDepth";
		        						}
		        					for(AlleleInfo ai: ailist) {
		        						if(ai.count() < this.min_genotype_allele_depth) {
			        						ai.masked = true;
			        						}
		        						}
	        						
		        					
		        		    		if(ailist.size()==2 ) {
		        		    			if(ailist.get(1).allele.isReference()) throw new IllegalArgumentException();
		        						double dp = ailist.get(0).count() + ailist.get(1).count();
		        						double fract = ailist.get(1).count()/dp;
		        						if(fract< this.min_genotype_fraction) {
		        							ailist.get(1).masked=true;
		        							}
		        						}
		        		    		si.allele2info.clear();
		        		    		for(AlleleInfo ai:ailist) {
		        		    			si.allele2info.put(ai.allele, ai);
		        		    			}
		        					}
		        				
		        				final Allele refAllele = Allele.create(prev_ref,true);
	        					final SortedSet<Allele> all_alleles_set = new TreeSet<>();
	        					all_alleles_set.add(refAllele);

		        				for(SampleInfo si: sampleid2info.values()) {
		        					for(AlleleInfo ai:si.allele2info.values()) {
		        						if(ai.masked) continue;
		        						all_alleles_set.add(ai.allele);
		        						}
		        					}
		        				all_alleles_set.remove(Allele.NO_CALL);
		        				if(all_alleles_set.size()>1) {
		        					int an = 0;
		        					int sum_gq = 0;
		        					final Counter<Allele> ac_counter = new Counter<>();
		        					final List<Allele> all_alleles_list = new ArrayList<>(all_alleles_set);
		        					if(!all_alleles_list.get(0).isReference()) {
		        						System.err.println(all_alleles_list);
		        						throw new IllegalStateException();
		        						}
		        					int ctx_dp  = 0;
			        				final VariantContextBuilder vcb = new VariantContextBuilder(null,
			        						interval.getContig(),
			        						prev_pos,
			        						prev_pos + prev_ref.length()-1,
			        						all_alleles_list
			        						);
			        				final List<Genotype> genotypes = new ArrayList<>(sampleid2info.size());

			        				/* create genotypes */
		        					for(SampleInfo si: sampleid2info.values()) {
		        						
		        						final List<AlleleInfo> ailist = si.allele2info.values().stream().
		        								filter(AI->!AI.masked).
		        								collect(Collectors.toCollection(ArrayList::new));
		        						final List<Allele> gt_alleles = new ArrayList<>(PLOIDY);
		        						
			        		    		if(ailist.isEmpty()) {
				        		    		gt_alleles.add(Allele.NO_CALL);
				        		    		gt_alleles.add(Allele.NO_CALL);
			        		    			}
			        		    		else if(ailist.size()==1) {
			        		    			final Allele alt= ailist.get(0).allele;
			        		    			if(alt.isReference()) {
			        		    				gt_alleles.add(refAllele);
			        		    				gt_alleles.add(refAllele);
			        		    				}
			        		    			else
			        		    				{
			        		    				gt_alleles.add(alt);
			        		    				gt_alleles.add(alt);
			        		    				}
			        		    			}
			        		    		else//  if(ailist.size()==2)
			        		    			{
			        		    			gt_alleles.add(ailist.get(0).allele);
			        		    			gt_alleles.add(ailist.get(1).allele);
			        		    			float dp = ailist.get(0).count()+ailist.get(1).count();
			        		    			float f = ailist.get(1).count()/dp;
			        		    			if(f<this.bad_ad_ratio || f>(1.0f-bad_ad_ratio)) {
			        		    				si.gtFilter="LowQual";
			        		    				}
			        		    			}
			        		    		for(Allele a: gt_alleles) {
			        		    			if(a.isNoCall()) continue;
			        		    			an++;
			        		    			ac_counter.incr(a);
			        		    			}
		        						final int[] ad  = new int[all_alleles_list.size()];
		        						int dp = 0;
		        						double gq = 0;
		        						Arrays.fill(ad, 0);

			        		    		for(int idx=0;idx< all_alleles_list.size();++idx) {
			        		    			final Allele a = all_alleles_list.get(idx);
			        		    			final AlleleInfo ai = si.allele2info.get(a);
			        		    			if(ai==null) continue;
			        		    			ad[idx] = ai.count();
			        		    			if(ai.masked) continue;
			        		    			dp += ai.count();
			        		    			ctx_dp += dp;
			        		    			gq += ai.mapq_sum;
			        		    			
			        		    			}
			        		    		final GenotypeBuilder gb = new GenotypeBuilder(si.sn, gt_alleles);
			        		    		if(!StringUtils.isBlank(si.gtFilter)) {
			        		    			gb.filter(si.gtFilter);
			        		    			}
			        		    		else if(((double)si.clip_sum/(dp+si.clip_sum)) >= 0.1) {
			        		    			gb.filter("HighClip");
			        		    			}
			        		    		gb.AD(ad);
			        					gb.DP(dp);
			        					final int gq2 = (int)(gq/(double)dp);
			        					sum_gq += gq2;
			        					gb.GQ(gq2);
			        		    		genotypes.add(gb.make());
			        					}
		        					vcb.attribute(VCFConstants.DEPTH_KEY, ctx_dp);
		        					vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,an);
		        					int[] ac  = new int[all_alleles_list.size()-1];
		        					double[] af  = new double[all_alleles_list.size()-1];
		        					Arrays.fill(ac, 0);
		        					Arrays.fill(af, 0.0);
		        					for(int idx=1;idx< all_alleles_list.size();++idx) {
		        						ac[idx-1] = (int)ac_counter.count(all_alleles_list.get(idx));
		        						if(an>0) {
		        							af[idx-1] = ac[idx-1]/(double)an;
		        							}
		        						}
		        					vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,ac);
		        					if(an>0) {
		        						vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,af);
		        						}
		        					if(sum_gq>0) vcb.log10PError(sum_gq/-10.0);
		        					vcb.genotypes(genotypes);
		        					vcw.add(vcb.make());
			        				}
	        					}
	        				if(c==null) break;
	        				prev_pos = c.position;
	        				prev_ref = c.ref;
	        				sampleid2info.clear();
	        				}
	        			SampleInfo si = sampleid2info.get(c.sample_index);
	        			if(si==null) {
	        				si = new SampleInfo(samplesList.get(c.sample_index),Allele.create(c.ref, true));
	        				sampleid2info.put(c.sample_index, si);
	        				}
	        			si.visit(c);
	        			}
		            }//end write vcf
	        	}
        	
        	
        	
        	return 0;
	    	} catch(Throwable err) {
	    		LOG.error(err);
	    		return -1;
	    	}
	    }

    public static void main(String[] args)
        {
        new MiniCaller().instanceMainWithExit(args);
        }

    }
