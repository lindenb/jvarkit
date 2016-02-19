/**
 * Author:
 * 	Pierre Lindenbaum PhD
 * Date:
 * 	Fev-2014
 * Contact:
 * 	plindenbaum@yahoo.fr
 * Motivation:
 * 	Idea from Solena: successive synonymous mutations are a stop codong
 */
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;



/**
 * VCFStopCodon
 * @SolenaLS 's idea: variant in the same codon give a new Amino acid undetected by annotaion tools.
 *
 */
public class VCFCombineTwoSnvs extends AbstractVCFCombineTwoSnvs
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VCFCombineTwoSnvs.class);
	/** known Gene collection */
	private final IntervalTreeMap<KnownGene> knownGenes=new IntervalTreeMap<KnownGene>();
	/** reference genome */
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	/** current genomic sequence */
	private GenomicSequence genomicSequence=null;
	/** all variants */
	private SortingCollection<Variant> variants= null;
	/** genetic Code used */
	private static final GeneticCode GENETIC_CODE = GeneticCode.getStandard();
	/** when INFO field for alleles are copied, we use that prefix */
	private static final String INFO_ALLELE_PREFIX="INFO_";
	/** FILTER header */
	private final VCFFilterHeaderLine vcfFilterHeaderLine = new VCFFilterHeaderLine(
			"TwoStrands",
			"(number of reads carrying both mutation) < (reads carrying variant 1 + reads carrying variant 2) "
			);
	/** when reading SAM file to get phasing */
	private SamReader samReader =null;

	
	
	/** mutated cdna */
	private static class MutedSequence extends DelegateCharSequence
		{
		private int begin=-1;
		private int end=-1;
		private String newseq=null;
	
		MutedSequence(CharSequence wild)
			{
			super(wild);
			}
		
		void setMutation(int begin,int end,final String newseq)
			{
			if(this.newseq!=null) throw new IllegalStateException();
			this.newseq = newseq;
			this.begin=begin;
			this.end=end;
			if(this.begin>this.end) throw new IllegalArgumentException();
			if(this.end> getDelegate().length()) throw new IndexOutOfBoundsException();
			}
		@Override
		public int length() {
			int L = getDelegate().length();
			if(this.newseq!=null) {
				L-=(this.end-this.begin);
				L+=this.newseq.length();
				}
			return L;
			}
		@Override
		public char charAt(int index)
			{
			if(this.newseq==null || index < this.begin ) return getDelegate().charAt(index);
			int idx2= index-this.begin;
			if(idx2 < this.newseq.length())
				{
				return this.newseq.charAt(idx2);
				}
			idx2-=this.newseq.length();
			return getDelegate().charAt(this.end+idx2);
			}
		
		}
	
	
	private static class ProteinCharSequence extends DelegateCharSequence
		{
		ProteinCharSequence(final CharSequence cDNA)
			{
			super(cDNA);
			}
		
		@Override
		public char charAt(int i)
			{
			return GENETIC_CODE.translate(
				getDelegate().charAt(i*3+0),
				getDelegate().charAt(i*3+1),
				getDelegate().charAt(i*3+2));
			}	
		
		@Override
		public int length()
			{
			return getDelegate().length()/3;
			}
		
		}
	
	/** load KnownGenes */
	private void loadKnownGenesFromUri() throws IOException
		{
		BufferedReader in = null;
		try {
			final SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
	        if(dict==null) throw new IOException("dictionary missing");

			LOG.info("loading genes from "+this.kgURI);
			in =IOUtils.openURIForBufferedReading(this.kgURI);
			final Pattern tab=Pattern.compile("[\t]");
			String line = null;
			while((line=in.readLine())!=null)
				{
				line= line.trim();
				if(line.isEmpty()) continue;
				String tokens[]=tab.split(line);
				final KnownGene g=new KnownGene(tokens);
				if(g.isNonCoding()) continue;
				if(dict.getSequence(g.getContig())==null)
					{
					LOG.warn("Unknown chromosome "+g.getContig()+" in dictionary");
					continue;
					}
				//use 1 based interval
				final Interval interval=new Interval(
						g.getContig(),
						g.getTxStart()+1,
						g.getTxEnd()
						);
				this.knownGenes.put(interval, g);
				}
			CloserUtil.close(in);in=null;
			LOG.info("genes:"+knownGenes.size());
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	
	
	private static int ID_GENERATOR=0;
	
	static private class AbstractContext {
		String contig;
		int genomicPosition1=0;
		String id=VCFConstants.EMPTY_ID_FIELD;
		Allele refAllele;
		Allele altAllele;
		int sorting_id;
		}
	
	static private class Variant extends AbstractContext
		{
		String transcriptName;
		byte strand;
		int position_in_cdna=-1;
		String wildCodon=null;
		String mutCodon=null;
		/** data from INFO column for this allele */
		Map<String,String> alleleInfo=null;
		
		Variant()
			{
			
			}
		
		Variant(final VariantContext ctx,final Allele allele,final KnownGene gene) {
			this.contig = ctx.getContig();
			this.genomicPosition1=ctx.getStart();
			this.id = (ctx.hasID()?ctx.getID():VCFConstants.EMPTY_ID_FIELD);
			this.transcriptName = gene.getName();
			this.strand = (byte)(gene.isNegativeStrand()?'-':'+');
			this.refAllele = ctx.getReference();
			this.altAllele = allele;
			}
		int positionInCodon() { return  position_in_cdna%3;}
		int codonStart() { return this.position_in_cdna - this.positionInCodon();}
		
		/** get specific info for this variant */
		private void _getInfo(final  Map<String,Object> info,int suffix) {
			//data about this
			info.put("POS"+suffix,String.valueOf(this.genomicPosition1));
			info.put("ID"+suffix,this.id);
			info.put("PosInCodon"+suffix,String.valueOf(this.positionInCodon()+1));
			info.put("Alt"+suffix,this.altAllele.getBaseString());
			info.put("Codon"+suffix,this.mutCodon);
			info.put("AA"+suffix,new ProteinCharSequence(this.mutCodon).getString());
				
			for(final String k:this.alleleInfo.keySet())
				{
				info.put(INFO_ALLELE_PREFIX+k+"_"+suffix,
						this.alleleInfo.get(k));
				}
			}
		
		public Map<String,Object> getInfo(final Variant other,final Variant third)
			{
			if(third!=null && third.position_in_cdna<other.position_in_cdna)
				{
				return getInfo(third,other);
				}
			final Map<String,Object> info = new LinkedHashMap<>();
			info.put("CHROM",this.contig);
			info.put("REF",this.refAllele.getBaseString());
			info.put("TRANSCRIPT",this.transcriptName);
			info.put("STRAND",this.strand==(byte)'+'?"plus":"minus");
			info.put("cDdnaPos",String.valueOf(this.position_in_cdna+1));
			info.put("CodonPos",String.valueOf(this.codonStart()+1));
			info.put("CodonWild",this.wildCodon);
			info.put("AAPos",String.valueOf(1+(this.codonStart()/3)));
			info.put("AAWild",new ProteinCharSequence(this.wildCodon).getString() );
			_getInfo(info, 1);
			other._getInfo(info, 2);
			
			if(third!=null) {
				third._getInfo(info, 3);
				}
			
			return info;
			}
		
		public Map<String,Object> getInfo(final Variant other)
			{
			return getInfo(other,null);
			}
			
		
		@Override
		public String toString() {
			return contig+"\t"+genomicPosition1+"\t"+refAllele.getBaseString()+"\t"+
					altAllele.getBaseString()+"\t"+
					transcriptName+"\t"+
					position_in_cdna+"\t"+
					codonStart()+"\t"+
					positionInCodon()+"\t"+
					wildCodon+"\t"+
					mutCodon
					;
			}
		
		}
	static private class VariantCodec extends AbstractDataCodec<Variant>
		{
		@Override
		public Variant decode(DataInputStream dis) throws IOException {
			String contig;
			try {
				contig = dis.readUTF();
			} catch (Exception e) {
				return null;
				}
			final Variant variant = new Variant();
			variant.contig = contig;
			variant.genomicPosition1 = dis.readInt();
			variant.id = dis.readUTF();
			variant.transcriptName = dis.readUTF();
			variant.strand = dis.readByte();
			variant.refAllele = Allele.create(dis.readUTF(), true);
			variant.altAllele = Allele.create(dis.readUTF(), false);
			variant.position_in_cdna = dis.readInt();
			variant.wildCodon = dis.readUTF();
			variant.mutCodon = dis.readUTF();
			variant.sorting_id = dis.readInt();
			
			int n = dis.readInt();
			variant.alleleInfo =new HashMap<>(n);
			for(int i=0;i< n;++i) {
				final String k = dis.readUTF();
				final String v = dis.readUTF();
				variant.alleleInfo.put(k, v);
			}
			
			return variant;
		}

		@Override
		public void encode(DataOutputStream dos, Variant v) throws IOException {
			dos.writeUTF(v.contig);
			dos.writeInt(v.genomicPosition1);
			dos.writeUTF(v.id);
			dos.writeUTF(v.transcriptName);
			dos.writeByte(v.strand);
			dos.writeUTF(v.refAllele.getBaseString());
			dos.writeUTF(v.altAllele.getBaseString());
			dos.writeInt(v.position_in_cdna);
			dos.writeUTF(v.wildCodon);
			dos.writeUTF(v.mutCodon);
			dos.writeInt(v.sorting_id);
			dos.writeInt(v.alleleInfo.size());
			for(final String k:v.alleleInfo.keySet()) {
				dos.writeUTF(k);
				dos.writeUTF(v.alleleInfo.get(k));
			}
		}

		@Override
		public VariantCodec clone() {
			return new VariantCodec();
		}
		
		}
	static private class VariantComparator implements Comparator<Variant>
		{
		SAMSequenceDictionary dict;
		VariantComparator(SAMSequenceDictionary dict) {
			this.dict = dict;
		}
		int contig(final Variant v) { return dict.getSequenceIndex(v.contig);}
		@Override
		public int compare(Variant o1, Variant o2) {
			int i= contig(o1) - contig(o2);
			if(i!=0) return i;
			i= o1.transcriptName.compareTo(o2.transcriptName);
			if(i!=0) return i;
			i= o1.position_in_cdna-o2.position_in_cdna;
			if(i!=0) return i;
			return o1.sorting_id - o2.sorting_id;
			}
		}
	
	
	private static class CombinedMutation extends AbstractContext
		{
		String info=null;
		String filter = VCFConstants.UNFILTERED;
		int depth =0;
		}

	static private class MutationCodec extends AbstractDataCodec<CombinedMutation>
		{
		@Override
		public CombinedMutation decode(final DataInputStream dis) throws IOException {
			String contig;
			try {
				contig = dis.readUTF();
			} catch (Exception e) {
				return null;
				}
			final CombinedMutation mutation = new CombinedMutation();
			mutation.contig = contig;
			mutation.genomicPosition1 = dis.readInt();
			mutation.id = dis.readUTF();
			mutation.refAllele = Allele.create(dis.readUTF(), true);
			mutation.altAllele = Allele.create(dis.readUTF(), false);
			mutation.info = dis.readUTF();
			mutation.sorting_id = dis.readInt();
			mutation.filter = dis.readUTF();
			mutation.depth = dis.readInt();
			return mutation;
		}
	
		@Override
		public void encode(final DataOutputStream dos, final CombinedMutation v) throws IOException {
			dos.writeUTF(v.contig);
			dos.writeInt(v.genomicPosition1);
			dos.writeUTF(v.id);
			dos.writeUTF(v.refAllele.getBaseString());
			dos.writeUTF(v.altAllele.getBaseString());
			dos.writeUTF(v.info);
			dos.writeInt(v.sorting_id);
			dos.writeUTF(v.filter);
			dos.writeInt(v.depth);
		}
	
		@Override
		public MutationCodec clone() {
			return new MutationCodec();
		}
		
		}
	
static private class MutationComparator implements Comparator<CombinedMutation>
	{
	final SAMSequenceDictionary dict;
	MutationComparator(SAMSequenceDictionary dict) {
		this.dict = dict;
	}
	int contig(final CombinedMutation v) { return dict.getSequenceIndex(v.contig);}
	@Override
	public int compare(CombinedMutation o1, CombinedMutation o2) {
		int i= contig(o1) - contig(o2);
		if(i!=0) return i;
		i= o1.genomicPosition1-o2.genomicPosition1;
		if(i!=0) return i;
		i =  o1.refAllele.compareTo(o2.refAllele);
		if(i!=0) return i;
		return o1.sorting_id - o2.sorting_id;
		}
	}

	/* used to compare a base in a BAM and an VCF allele */
	private static boolean same(byte baseInRead,final Allele a)
		{
		return Character.toUpperCase(a.getBases()[0])==Character.toUpperCase(baseInRead);
		}
	/* used to serialize a map to string in INFO column */
	private static String mapToString(final Map<String,Object> map)
		{
		final StringBuilder sb = new StringBuilder();
		for(final String key: map.keySet()) {
			if(sb.length()>0) sb.append("|");
			sb.append(key);
			sb.append("|");
			sb.append(String.valueOf(map.get(key)));
		}
		return sb.toString();
		}

	@Override
	protected Collection<Throwable> doVcfToVcf(final String inputName,final VcfIterator r,final VariantContextWriter w) throws IOException {
			SortingCollection<CombinedMutation> mutations = null;
			CloseableIterator<Variant> varIter = null;
			CloseableIterator<CombinedMutation> mutIter = null;
			try {
			LOG.info("opening REF:"+referenceFile);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(this.referenceFile);
	        final SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
	        if(dict==null) throw new IOException("dictionary missing");
	        
	        
	        final VCFInfoHeaderLine depthInfoHeaderLine = (super.bamIn == null ?
	        		null:
	        		new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth")
	        		);
	        
	        if(super.bamIn!=null)
	        	{
	        	LOG.info("opening BAM :"+super.bamIn);
	        	this.samReader = SamReaderFactory.makeDefault().
	        			referenceSequence(this.referenceFile).
	        			validationStringency(ValidationStringency.LENIENT).
	        			open(super.bamIn)
	        			;
	        	if(!this.samReader.hasIndex())
	        		{
	        		 throw new IOException("Sam file is NOT indexed.");
	        		}
	        	if(this.samReader.getFileHeader().getSequenceDictionary()==null || 
	        		!SequenceUtil.areSequenceDictionariesEqual(dict, samReader.getFileHeader().getSequenceDictionary())) {
	        		 throw new IOException("SamReader and REF don't have the same Sequence Dictionary.");
	        		}
	        	}
	        
	        
			loadKnownGenesFromUri();
			
			this.variants = SortingCollection.newInstance(Variant.class,
					new VariantCodec(),
					new VariantComparator(dict),
					super.getMaxRecordsInRam(),
					super.getTmpDirectories()
					);
			this.variants.setDestructiveIteration(true);
			
			/* get VCF header */
			final VCFHeader header=(VCFHeader)r.getHeader();
			/* collect VCF INFO describing ALT alleles */
			final List<VCFInfoHeaderLine> alleleInfoHeaderLines = 
					super.grabInfoA ?
					header.getInfoHeaderLines().stream().
					filter(INFO->INFO.getCountType()==VCFHeaderLineCount.A).
					collect(Collectors.toList()):
					Collections.emptyList()
					;
			
			
	       SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			while(r.hasNext())
				{
				final VariantContext ctx= progress.watch(r.next());
				/* discard non SNV variant */
				if(!ctx.isVariant() || ctx.isIndel())
					{
					continue;
					}
				
				/* find the overlapping genes : extend the interval of the variant to include the stop codon */
				final Collection<KnownGene> genes= this.knownGenes.getOverlapping(
						new Interval(ctx.getContig(),
						Math.max(1,ctx.getStart()-3),
						ctx.getEnd()+3
						));
				
				final List<Allele> alternateAlleles =  ctx.getAlternateAlleles();
				
				/* loop over overlapping genes */
				for(final KnownGene kg:genes) {
					/* loop over available alleles */
					for(int allele_idx=0;allele_idx< alternateAlleles.size();++allele_idx) {
						final Allele alt = alternateAlleles.get(allele_idx);
						
						/* get specific INFO for this allele */
						final Map<String,String> alleleInfo = new HashMap<>(alleleInfoHeaderLines.size());
						for(final VCFInfoHeaderLine vihl:alleleInfoHeaderLines)
							{
							if(!ctx.hasAttribute(vihl.getID())) continue;
							final List<Object> atts = ctx.getAttributeAsList(vihl.getID());
							if(allele_idx< atts.size()) 
								{
								alleleInfo.put(vihl.getID(), String.valueOf(atts.get(allele_idx)));
								}
							}
						
						challenge(ctx,alt,kg,alleleInfo);
						}
					}
				}
			progress.finish();
			this.variants.doneAdding();
			
			mutations = SortingCollection.newInstance(CombinedMutation.class,
					new MutationCodec(),
					new MutationComparator(dict),
					super.getMaxRecordsInRam(),
					super.getTmpDirectories()
					);
			mutations.setDestructiveIteration(true);
			
			
			
			varIter = this.variants.iterator();
			progress=new SAMSequenceDictionaryProgress(header);
			final ArrayList<Variant> buffer= new ArrayList<>();
			for(;;)
				{
				Variant variant = null;
				if(varIter.hasNext())
					{
					variant = varIter.next();
					progress.watch(variant.contig, variant.genomicPosition1);
					}
				if(variant==null || !(!buffer.isEmpty() &&
							buffer.get(0).contig.equals(variant.contig) && 
							buffer.get(0).transcriptName.equals(variant.transcriptName)))
					{
					searchSnvInSameCodon(mutations,buffer);	
					buffer.clear();
					if(variant==null) break;
					}
				buffer.add(variant);
				}
			progress.finish();
			mutations.doneAdding();
			varIter.close();varIter=null;
			variants.cleanup();variants=null;
			final ArrayList<CombinedMutation> mBuffer= new ArrayList<>();
			
			final VCFHeader header2 = new VCFHeader();
			header2.setSequenceDictionary(header.getSequenceDictionary());
			
			final StringBuilder infoDesc =new StringBuilder("Variant affected by two distinct mutation. Format is defined in the INFO column. ");
			
			/* copy original 'A' info */
			for(final VCFInfoHeaderLine ihl: alleleInfoHeaderLines)
				{
				infoDesc.append(INFO_ALLELE_PREFIX+ihl.getID()+":"+ihl.getDescription()+".");
				}
			
			final VCFInfoHeaderLine infoHeaderLine = new VCFInfoHeaderLine(
					"CodonVariant",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					infoDesc.toString()
					);
			
			
			
			super.addMetaData(header2);
			header2.addMetaDataLine(infoHeaderLine);
			if(depthInfoHeaderLine!=null)
				{
				header2.addMetaDataLine(depthInfoHeaderLine);
				}
			if(samReader!=null)
				{
				header2.addMetaDataLine(this.vcfFilterHeaderLine);
				}
			w.writeHeader(header2);
			
			progress=new SAMSequenceDictionaryProgress(header);
			mutIter = mutations.iterator();
			for(;;)
				{
				CombinedMutation mutation = null;
				if(mutIter.hasNext())
					{
					mutation = mutIter.next();
					progress.watch(mutation.contig, mutation.genomicPosition1);
					}
				if(mutation==null || !(!mBuffer.isEmpty() &&
						mBuffer.get(0).contig.equals(mutation.contig) &&  
						mBuffer.get(0).genomicPosition1 == mutation.genomicPosition1 &&  
						mBuffer.get(0).refAllele.equals(mutation.refAllele)))
					{
					if(!mBuffer.isEmpty())
						{
						//default filter fails
						String filter= this.vcfFilterHeaderLine.getID();
						final CombinedMutation first  = mBuffer.get(0);
						final Set<Allele> alleles = new HashSet<>();
						final Set<String> info = new HashSet<>();
						alleles.add(first.refAllele);
						final VariantContextBuilder vcb=new VariantContextBuilder();
						vcb.chr(first.contig);
						vcb.start(first.genomicPosition1);
						vcb.stop(first.genomicPosition1 + first.refAllele.length()-1);
						if( !first.id.equals(VCFConstants.EMPTY_ID_FIELD)) vcb.id(first.id);
						for(final CombinedMutation m:mBuffer){
							alleles.add(m.altAllele);
							info.add(m.info);
							if(VCFConstants.UNFILTERED.equals(m.filter)) {
								filter = null; //at least one SNP is ok one this line
							}
						}
						if(filter!=null && samReader!=null) vcb.filter(filter);
						vcb.attribute(infoHeaderLine.getID(), new ArrayList<String>(info));
						if(depthInfoHeaderLine!=null)
							{
							vcb.attribute(depthInfoHeaderLine.getID(),first.depth);
							}
						
						vcb.alleles(alleles);
						w.add(vcb.make());
						}
					mBuffer.clear();
					if(mutation==null) break;
					}
				mBuffer.add(mutation);
				}
			progress.finish();
			mutIter.close();
			mutations.cleanup();mutations=null;
			
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(mutIter);
			CloserUtil.close(varIter);
			
			if(this.variants!=null) this.variants.cleanup();
			if(mutations!=null) mutations.cleanup();
			this.variants=null;
			CloserUtil.close(this.samReader);
			this.samReader=null;
			}
		}
	
	private void searchSnvInSameCodon(final SortingCollection<CombinedMutation> mutations,final List<Variant> buffer)
		{
		if(buffer.size()<2) return;
		for(int i=0;i< buffer.size();++i)
			{
			final Variant v1  = buffer.get(i);
			for(int j=i+1;j< buffer.size();++j)
				{
				final Variant v2  = buffer.get(j);
				if(v1.codonStart() != v2.codonStart()) continue;
				if(v1.positionInCodon() == v2.positionInCodon()) continue;
				if(!v1.wildCodon.equals(v2.wildCodon))
					{
					throw new IllegalStateException();
					}
				if(foundVariantsInTheSameCodon(mutations,new Variant[]{v1,v2})) {
					//we already found something with those two variants
					continue;
				}
				
				//optional 3rd base in the same codon
				for(int k=j+1;k< buffer.size();++k)
					{
					final Variant v3  = buffer.get(k);
					if(v1.codonStart() != v3.codonStart()) continue;
					if(v1.positionInCodon() == v3.positionInCodon()) continue;
					if(v2.positionInCodon() == v3.positionInCodon()) continue;
					if(!v1.wildCodon.equals(v3.wildCodon))
						{
						throw new IllegalStateException();
						}
					foundVariantsInTheSameCodon(mutations,new Variant[]{v1,v2,v3});
					}
				}
			}
		}
	
	private boolean foundVariantsInTheSameCodon(
			SortingCollection<CombinedMutation> mutations,
			final Variant variants[]) 	{
		final Variant firstVariant = variants[0];
		int chromStart = firstVariant.genomicPosition1; 
		int chromEnd = firstVariant.genomicPosition1; 
		
				
		final StringBuilder combinedCodon = new StringBuilder(firstVariant.wildCodon);
		final String mutedAminoAcids[]= new String[variants.length];
		final String pwild = new ProteinCharSequence(firstVariant.wildCodon).getString();

		boolean mutations_is_synonymous=true;
		boolean mutations_same_as_combined=true;
		
		for(int i=0;i< variants.length;++i) {
			final Variant v=variants[i];
			combinedCodon.setCharAt(v.positionInCodon(), v.mutCodon.charAt(v.positionInCodon()));
			mutedAminoAcids[i] = new ProteinCharSequence(v.mutCodon).getString();
			if(!pwild.equals(mutedAminoAcids[i] )) {
				mutations_is_synonymous = false;
				}
			chromStart = Math.min(chromStart,v.genomicPosition1);
			chromEnd = Math.max(chromEnd,v.genomicPosition1);
			}
		
		final String pCombined = new ProteinCharSequence(combinedCodon).getString();
		for(int i=0;i< variants.length;++i)
			{
			if(!pCombined.equals(mutedAminoAcids[i])) {
				mutations_same_as_combined = false;
				}
			}
		
		final String combinedSO;
		final String combinedType;
		
		/* both AA are synonymous, while combined is not */
		if(!pCombined.equals(pwild) &&  mutations_is_synonymous ) {
			combinedType = "combined_is_nonsynonymous";
			if(pCombined.equals("*"))
				{
				/* http://www.sequenceontology.org/browser/current_svn/term/SO:0001587 */
				combinedSO="stop_gained";
				}
			else if(pwild.equals("*"))
				{
				/* http://www.sequenceontology.org/browser/current_svn/term/SO:0002012 */
				combinedSO="stop_lost";
				}
			else
				{
				/* http://www.sequenceontology.org/miso/current_svn/term/SO:0001992 */
				combinedSO="nonsynonymous_variant";
				}
			}
		else if(!mutations_same_as_combined && !pCombined.equals(pwild))
			{
			combinedType = "combined_is_new";
			if(pCombined.equals("*"))
				{
				/* http://www.sequenceontology.org/browser/current_svn/term/SO:0001587 */
				combinedSO="stop_gained";
				}
			else
				{
				/* http://www.sequenceontology.org/miso/current_svn/term/SO:0001992 */
				combinedSO="nonsynonymous_variant";
				}
			}
		else
			{
			//System.err.println("ignored "+pCombined+" vs "+pwild+" "+Arrays.asList(mutedAminoAcids));
			return false;
			}
			/** ok, there is something interesting here ,
			 * create two new Mutations carrying the
			 * two variants 
			 */
			final List<Map<String,Object>> infos=new ArrayList<>(variants.length);
		
			if(variants.length==2)
				{
				infos.add( variants[0].getInfo(variants[1]));
				infos.add( variants[1].getInfo(variants[0]));
				}
			else if(variants.length==3) {
				infos.add( variants[0].getInfo(variants[1],variants[2]));
				infos.add( variants[1].getInfo(variants[0],variants[2]));
				infos.add( variants[2].getInfo(variants[0],variants[1]));
				}
			else
				{
				throw new RuntimeException();
				}
			
			
			/** info that will be displayed in the vcf */
			//filter for this combined: default it fails the filter
			String filter = this.vcfFilterHeaderLine.getID();
			
			final Map<String,Object> combinedMap = new LinkedHashMap<>();
			combinedMap.put("CombinedCodon",combinedCodon);
			combinedMap.put("CombinedAA",pCombined);
			combinedMap.put("CombinedSO",combinedSO);
			combinedMap.put("CombinedType",combinedType);
			
			for(final Map<String,Object> info:infos)
				{
				info.putAll(combinedMap);
				}			
			/** depth at position */
			final int depthList[]= new int[variants.length];
			Arrays.fill(depthList,0);
			
			/* get phasing info */
			if(this.samReader!=null) {
				int count_reads_having_all_variants = 0;
				int count_reads_having_no_variants = 0;
				int count_reads_having_ONLY_variant[] = new int[variants.length];
				Arrays.fill(count_reads_having_ONLY_variant,0);


				final SAMRecordIterator iter = samReader.query(
						firstVariant.contig,
						chromStart,
						chromEnd,
						false
						);
				try {
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(rec.isSecondaryOrSupplementary()) continue;
						if(rec.getDuplicateReadFlag()) continue;
						if(rec.getReadFailsVendorQualityCheckFlag()) continue;
						
						// get DEPTh for variant x
						for(int x=0;x < variants.length;++x)
							{
							if(rec.getAlignmentStart()<= variants[x].genomicPosition1 && variants[x].genomicPosition1<=rec.getAlignmentEnd()) {
								depthList[x]++;
								}
							}
						
						if(rec.getAlignmentEnd()<chromEnd) continue;
						if(rec.getAlignmentStart()>chromStart) continue;
						final Cigar cigar  =  rec.getCigar();
						if(cigar==null) continue;
						final byte bases[] = rec.getReadBases();
						if(bases==null) continue;
						int refpos1=rec.getAlignmentStart();
						int readpos = 0;
						final boolean found_variant_on_this_read[] = new boolean[variants.length];
						Arrays.fill(found_variant_on_this_read, false);
						
						/** loop over cigar */
						for(final CigarElement ce:cigar.getCigarElements()) {
							final CigarOperator op = ce.getOperator();
							switch(op)
								{
								case P: continue;
								case S: case I: readpos+=ce.getLength();break;
								case D: case N: refpos1+=ce.getLength(); break;
								case H: continue;
								case EQ:case M:case X:
									for(int x=0;x< ce.getLength();++x)
										{
										for(int vidx=0;vidx< variants.length;++vidx)
											{
											final Variant v=variants[vidx];
											if(refpos1 == v.genomicPosition1 && same(bases[readpos],v.altAllele))
												{
												found_variant_on_this_read[vidx] = true;
												}
											}
										refpos1++;
										readpos++;
										}
									break;
								default: throw new IllegalStateException(op.name());
								}
							/* skip remaining bases after last variant */
							if(refpos1>chromEnd) break;
						}/* end of loop over cigar */
						
					/* sum-up what we found */
					int count_having_found_variant_on_this_read=0; 
					for(boolean b:found_variant_on_this_read) count_having_found_variant_on_this_read+=(b?1:0);
					
						
					if( count_having_found_variant_on_this_read == variants.length) {
						count_reads_having_all_variants++; }
					else if(count_having_found_variant_on_this_read==0) {
						count_reads_having_no_variants++; }
					else 
						{
						for(int idx1=0;idx1 < found_variant_on_this_read.length;++idx1) {
							boolean only_for_this_variant=found_variant_on_this_read[idx1];
							for(int idx2=0;idx2 < found_variant_on_this_read.length;++idx2)
								{
								if(idx1==idx2) continue;
								if(found_variant_on_this_read[idx2]) {
									only_for_this_variant=false;
								}
								}
							if(only_for_this_variant) count_reads_having_ONLY_variant[idx1]++;
							}
						}					
					}/* end of loop over reads */
					
				} finally {
					iter.close();
				}
				
				for(final Map<String,Object> info:infos)
					{
					info.put("N_READS_ALL_VARIANTS", count_reads_having_all_variants);
					info.put("N_READS_NO_VARIANTS", count_reads_having_no_variants);
					info.put("N_READS_NO_VARIANTS", count_reads_having_no_variants);
					info.put("N_READS_TOTAL",
							count_reads_having_all_variants +
							count_reads_having_no_variants +
							IntStream.of(count_reads_having_ONLY_variant).sum()
							);
					}
				
				for(int idx1=0;idx1 < variants.length;++idx1)
					{
					infos.get(idx1).put("N_READS_ONLY_1", count_reads_having_ONLY_variant[idx1]);
					if(variants.length==2)
						{
						infos.get(idx1).put("N_READS_ONLY_2", count_reads_having_ONLY_variant[idx1==0?1:0]);
						}
					else if(variants.length==3)
						{
						List<Integer> others= new ArrayList<>(2);
						for(int idx2=0;idx2< variants.length;++idx2) if(idx2!=idx1) others.add(idx2);
						if(variants[others.get(0)].position_in_cdna < variants[others.get(1)].position_in_cdna)
							{
							infos.get(idx1).put("N_READS_ONLY_2", count_reads_having_ONLY_variant[others.get(0)]);
							infos.get(idx1).put("N_READS_ONLY_3", count_reads_having_ONLY_variant[others.get(1)]);
							}
						else
							{
							infos.get(idx1).put("N_READS_ONLY_2", count_reads_having_ONLY_variant[others.get(1)]);
							infos.get(idx1).put("N_READS_ONLY_3", count_reads_having_ONLY_variant[others.get(0)]);
							}
						}
						
					}
				
				
				/* number of reads with both variant is greater than
				 * reads carrying only one variant: reset the filter 
				 */
				if(count_reads_having_all_variants> IntStream.of(count_reads_having_ONLY_variant).sum() ) {
					filter = VCFConstants.UNFILTERED;
					}
				}
			
			for(int vidx=0;vidx< variants.length;++vidx)
				{
				final Variant v1 = variants[vidx];
				final CombinedMutation m1 = new CombinedMutation();
				m1.contig = v1.contig;
				m1.genomicPosition1 = v1.genomicPosition1;
				m1.id = v1.id ;
				m1.refAllele = v1.refAllele;
				m1.altAllele = v1.altAllele;
				m1.info = mapToString(infos.get(vidx));
				m1.filter = filter;
				m1.depth = depthList[vidx];
	
				m1.sorting_id = ID_GENERATOR++;
				mutations.add(m1);
				}
		return true;
		}
		
	
	
	private void challenge(
			final VariantContext ctx,
			final Allele allele,
			final KnownGene gene,
			final Map<String,String> alleleInfo
			) throws IOException
		{
		if(allele.isSymbolic()) return;
		if(allele.isNoCall()) return;
		if(!allele.isCalled())return;
		if(allele.length()!=1) return;
		if(ctx.getReference().length()!=1) return;
		if(allele.equals(Allele.SPAN_DEL)) return;
		
		if(genomicSequence==null || !genomicSequence.getChrom().equals(ctx.getContig()))
			{
			LOG.info("getting genomic Sequence for "+gene.getContig());
			genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, gene.getContig());
			}
		final int positionContext0 = ctx.getStart() -1;
		
		
		Variant variant=null;
		final StringBuilder wildRNA=new StringBuilder(1000);
		final MutedSequence mutRNA=new MutedSequence(wildRNA);
    	
		if(gene.isPositiveStrand())
			{
	    	for(int exon_index=0;exon_index< gene.getExonCount();++exon_index)
	    		{
	    		final KnownGene.Exon exon= gene.getExon(exon_index);
	    		int genomicPosition = Math.max(gene.getCdsStart() , exon.getStart());
	    		for(;;)
					{
					// we need to consider the last stop codon
					if( genomicPosition>=exon.getEnd() ||
						genomicPosition>= gene.getCdsEnd())
						{
						break;
						}

					wildRNA.append(genomicSequence.charAt(genomicPosition));
					
					if(variant==null && positionContext0 ==genomicPosition)
						{
						variant = new Variant(ctx,allele,gene);
						variant.sorting_id = ID_GENERATOR++;
						variant.position_in_cdna=wildRNA.length()-1;
						char mutBase= allele.getBaseString().charAt(0);
						mutRNA.setMutation( wildRNA.length()-1, wildRNA.length(),""+mutBase 	);
						}
					++genomicPosition;
					}
	    		}
			}
		else
			{
			int exon_index = gene.getExonCount()-1;
			while(exon_index >=0)
				{
				final KnownGene.Exon exon= gene.getExon(exon_index);
				int genomicPosition = Math.min(gene.getCdsEnd()-1, exon.getEnd()-1);
				for(;;) {					
					if( genomicPosition< exon.getStart() ||
						genomicPosition<gene.getCdsStart())
						{
						break;
						}
	
					wildRNA.append(AcidNucleics.complement(genomicSequence.charAt(genomicPosition)));
					
					if(variant==null && positionContext0 == genomicPosition )
						{
						variant = new Variant(ctx,allele,gene);
						variant.sorting_id = ID_GENERATOR++;
						variant.position_in_cdna=wildRNA.length()-1;
						char mutBase=AcidNucleics.complement(allele.getBaseString().charAt(0));
						mutRNA.setMutation( wildRNA.length()-1, wildRNA.length(),""+mutBase 	);
						}
					--genomicPosition;
    				}
    			--exon_index;
    			}
			}
        	
    	if(variant!=null)
    		{
    		variant.wildCodon="";
    		variant.mutCodon="";
    		for(int i=0;i< 3;++i)
    			{
    			int pos = variant.codonStart()+i;
    			variant.wildCodon += (pos< wildRNA.length()?wildRNA.charAt(pos):'*');
    			variant.mutCodon +=  (pos< mutRNA.length()?mutRNA.charAt(pos):'*');
    			}
    		variant.wildCodon = variant.wildCodon.toUpperCase();
    		variant.mutCodon = variant.mutCodon.toUpperCase();
    		variant.alleleInfo = alleleInfo;
    		
    		if(variant.wildCodon.equals(variant.mutCodon)) {
    			LOG.info("Uh??????? "+allele+" "+ctx);
    			return;
    			}
    		this.variants.add(variant);
    		}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(this.referenceFile==null)
			{
			return wrapException("Undefined REFERENCE. option: -"+OPTION_REFERENCEFILE);
			}
		if(this.kgURI==null || this.kgURI.trim().isEmpty())
			{
			return wrapException("Undefined kgURI. option: -"+OPTION_KGURI);
			}
		return doVcfToVcf(inputName);
		}
	
	public static void main(String[] args)
		{
		new VCFCombineTwoSnvs().instanceMainWithExit(args);
		}
	}
