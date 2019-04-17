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
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.WeakHashMap;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.Paranoid;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.GranthamScore;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
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
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**

BEGIN_DOC

## Motivation

@SolenaLS 's idea: variant in the same codon give a new Amino acid undetected by annotaion tools.

## Example

```
##fileformat=VCFv4.2
##FILTER=<ID=TwoStrands,Description="(number of reads carrying both mutation) < (reads carrying variant 1 + reads carrying variant 2)">
##INFO=<ID=CodonVariant,Number=.,Type=String,Description="Variant affected by two distinct mutation. Format is defined in the INFO column. INFO_AC:Allele count in genotypes, for each ALT allele, in the same order as listed.INFO_AF:Allele Frequency, for each ALT allele, in the same order as listed.INFO_MLEAC:Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed.INFO_MLEAF:Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed.">
##VCFCombineTwoSnvsCmdLine=-k jeter.knownGene.txt -tmpdir tmp/ -R /commun/data/pubdb/broadinstitute.org/bundle/1.5/b37/human_g1k_v37.fasta -B /commun/data/projects/plateforme/NTS-017_HAL_Schott_mitral/20141106/align20141106/Samples/CD13314/BAM/Haloplex20141106_CD13314_final.bam
##VCFCombineTwoSnvsHtsJdkHome=/commun/data/packages/htsjdk/htsjdk-2.0.1
##VCFCombineTwoSnvsHtsJdkVersion=2.0.1
##VCFCombineTwoSnvsVersion=c5af7d1bd367562b3578d427d24ec62856835d38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	120612013	rs200646249	G	A	.	.	CodonVariant=CHROM|1|REF|G|TRANSCRIPT|uc001eil.3|cDdnaPos|8|CodonPos|7|CodonWild|GCC|AAPos|3|AAWild|A|POS1|120612013|ID1|rs200646249|PosInCodon1|2|Alt1|A|Codon1|GTC|AA1|V|INFO_MLEAC_1|1|INFO_AC_1|1|INFO_MLEAF_1|0.500|INFO_AF_1|0.500|POS2|120612014|ID2|.|PosInCodon2|1|Alt2|A|Codon2|TCC|AA2|S|INFO_MLEAC_2|1|INFO_AC_2|1|INFO_MLEAF_2|0.500|INFO_AF_2|0.500|CombinedCodon|TTC|CombinedAA|F|CombinedSO|nonsynonymous_variant|CombinedType|combined_is_new|N_READS_BOTH_VARIANTS|168|N_READS_NO_VARIANTS|1045|N_READS_TOTAL|1213|N_READS_ONLY_1|0|N_READS_ONLY_2|0,CHROM|1|REF|G|TRANSCRIPT|uc001eik.3|cDdnaPos|8|CodonPos|7|CodonWild|GCC|AAPos|3|AAWild|A|POS1|120612013|ID1|rs200646249|PosInCodon1|2|Alt1|A|Codon1|GTC|AA1|V|INFO_MLEAC_1|1|INFO_AC_1|1|INFO_MLEAF_1|0.500|INFO_AF_1|0.500|POS2|120612014|ID2|.|PosInCodon2|1|Alt2|A|Codon2|TCC|AA2|S|INFO_MLEAC_2|1|INFO_AC_2|1|INFO_MLEAF_2|0.500|INFO_AF_2|0.500|CombinedCodon|TTC|CombinedAA|F|CombinedSO|nonsynonymous_variant|CombinedType|combined_is_new|N_READS_BOTH_VARIANTS|168|N_READS_NO_VARIANTS|1045|N_READS_TOTAL|1213|N_READS_ONLY_1|0|N_READS_ONLY_2|0;EXAC03_AC_NFE=641;EXAC03_AN_NFE=48948
1	120612014	.	C	A	.	.	CodonVariant=CHROM|1|REF|C|TRANSCRIPT|uc001eik.3|cDdnaPos|7|CodonPos|7|CodonWild|GCC|AAPos|3|AAWild|A|POS1|120612014|ID1|.|PosInCodon1|1|Alt1|A|Codon1|TCC|AA1|S|INFO_MLEAC_1|1|INFO_AC_1|1|INFO_MLEAF_1|0.500|INFO_AF_1|0.500|POS2|120612013|ID2|rs200646249|PosInCodon2|2|Alt2|A|Codon2|GTC|AA2|V|INFO_MLEAC_2|1|INFO_AC_2|1|INFO_MLEAF_2|0.500|INFO_AF_2|0.500|CombinedCodon|TTC|CombinedAA|F|CombinedSO|nonsynonymous_variant|CombinedType|combined_is_new|N_READS_BOTH_VARIANTS|168|N_READS_NO_VARIANTS|1045|N_READS_TOTAL|1213|N_READS_ONLY_1|0|N_READS_ONLY_2|0,CHROM|1|REF|C|TRANSCRIPT|uc001eil.3|cDdnaPos|7|CodonPos|7|CodonWild|GCC|AAPos|3|AAWild|A|POS1|120612014|ID1|.|PosInCodon1|1|Alt1|A|Codon1|TCC|AA1|S|INFO_MLEAC_1|1|INFO_AC_1|1|INFO_MLEAF_1|0.500|INFO_AF_1|0.500|POS2|120612013|ID2|rs200646249|PosInCodon2|2|Alt2|A|Codon2|GTC|AA2|V|INFO_MLEAC_2|1|INFO_AC_2|1|INFO_MLEAF_2|0.500|INFO_AF_2|0.500|CombinedCodon|TTC|CombinedAA|F|CombinedSO|nonsynonymous_variant|CombinedType|combined_is_new|N_READS_BOTH_VARIANTS|168|N_READS_NO_VARIANTS|1045|N_READS_TOTAL|1213|N_READS_ONLY_1|0|N_READS_ONLY_2|0;EXAC03_AC_NFE=640;EXAC03_AN_NFE=48228

```

#### Fields

```
KEYEXAMPLEDescription
CHROM1Chromosome for current variant.
REFCReference Allele for current variant
TRANSCRIPTuc001eik.3UCSC knownGene Transcript
cDdnaPos7+1 based position in cDNA
CodonPos7+1 based position of the codon in cNA
CodonWildGCCWild codon
AAPos3+1 based position of amino acid
AAWildAWild amino acid
POS1120612014+1 based position of variant 1
ID1.RS ID of variant 1
PosInCodon11Position in codon (1,2,3) of variant 1
Alt1AAlternate allele of variant 1
Codon1TCC Codon with variant 1 alone
AA1SAmino acid prediction for variant 1
INFO_*_11Data about alternate allele 1 taken out of original VCF
POS2120612013+1 based position of variant 1
ID2rs200646249RS ID of variant 1
PosInCodon22Position in codon (1,2,3) of variant 2
Alt2AAlternate allele of variant 2
Codon2GTC Codon with variant 2alone
AA2VAmino acid prediction for variant 2
INFO_*_21Data about alternate allele 2 taken out of original VCF
CombinedCodonTTCCombined codon with ALT1 and ALT2
CombinedAAFCombined amino acid with ALT1 and ALT2
CombinedSOnonsynonymous_variantSequence Ontology term
CombinedTypecombined_is_newtype of new mutation
N_READS_BOTH_VARIANTS168Number of reads carrying both variants
N_READS_NO_VARIANTS1045Number of reads carrying no variants
N_READS_TOTAL1213Total Number of reads
N_READS_ONLY_10Number of reads carrying onlt variant 1
N_READS_ONLY_20Number of reads carrying onlt variant 2
```

### See also

 * http://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-5-615
 * https://www.biorxiv.org/content/10.1101/573378v2

END_DOC
*/


@Program(name="vcfcombinetwosnvs",
	description="Detect Mutations than are the consequences of two distinct variants. This kind of variant might be ignored/skipped from classical variant consequence predictor. Idea from @SolenaLS and then @AntoineRimbert",
	keywords={"vcf","annotation","prediction","protein","mnv"},
	modificationDate="20190322"
	)
public class VCFCombineTwoSnvs extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFCombineTwoSnvs.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-k","--knownGene"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String kgURI  = KnownGene.getDefaultUri();
	@Parameter(names={"-B","--bam"},description="Optional indexed BAM file used to get phasing information. This can be a list of bam if the filename ends with '.list'")
	private Path bamIn = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true )
	private Path referencePath;
	@Parameter(names={"-P","--bedpe"},description="save optional report as bedpe")
	private Path bedPePath = null;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	
	/** known Gene collection */
	private final IntervalTreeMap<List<KnownGene>> knownGenes=new IntervalTreeMap<>();
	/** reference genome */
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	/** current genomic sequence */
	private GenomicSequence genomicSequence=null;
	/** all variants */
	private SortingCollection<Variant> variants= null;
	/** genetic Code used */
	private static final GeneticCode GENETIC_CODE = GeneticCode.getStandard();
	/**  map transcript to their sequence */
	private final WeakHashMap<String, CDNASequence> kgId2transcriptCache = new WeakHashMap<>();
	/** because 'I'm paranoid */
	private final Paranoid paranoiac = Paranoid.createThrowingInstance();
	
	private class CDNASequence extends AbstractCharSequence
		{
		final KnownGene knownGene;
		final StringBuilder sequence = new StringBuilder(1_000);
		final Map<Integer, Integer> genomic2cnda;
		CDNASequence(final KnownGene knownGene) {
			this.knownGene=knownGene;
			final List<Integer> genomicPositions = new ArrayList<>();
			
	    	for(int exon_index=0;exon_index< this.knownGene.getExonCount();++exon_index)
	    		{
	    		int pos0 = Math.max(this.knownGene.getCdsStart(),this.knownGene.getExonStart(exon_index));
	    		final int end0 =  Math.min(this.knownGene.getExonEnd(exon_index),this.knownGene.getCdsEnd());
	    		while(pos0<end0) {
	    			genomicPositions.add(pos0);
	    			char base = Character.toUpperCase(genomicSequence.charAt(pos0));
	    			this.sequence.append(base);
	    			pos0++;
	    			}
	    		}
	    	
	    	if(this.knownGene.isPositiveStrand()) {
	    		// nothing
	    		}
	    	else
	    		{
	    			
	    		Collections.reverse(genomicPositions);
	    		this.sequence.reverse();
	    		
	    		for(int i=0;i< this.sequence.length();i++) {
	    			this.sequence.setCharAt(i, AcidNucleics.complement(this.sequence.charAt(i)));
    				}
	    		}
	    	
	    	paranoiac.assertGt(this.sequence.length(), 3);
	    	
	    	this.genomic2cnda=new HashMap<>(genomicPositions.size());
	    	for(int i=0;i< genomicPositions.size();i++) {
	    		this.genomic2cnda.put(genomicPositions.get(i),i);
	    		}
	    	/*
	    	paranoiac.assertTrue(GENETIC_CODE.isStopCodon(
	    			this.sequence.charAt(this.sequence.length()-3),
	    			this.sequence.charAt(this.sequence.length()-2),
	    			this.sequence.charAt(this.sequence.length()-1)
	    			),knownGene.getName()+" "+knownGene.getStrand()+" "+knownGene.getContig()+" "+knownGene.getCdsStart()+" "+knownGene.getCdsEnd());
	    	*/
			}
		@Override
		public char charAt(int index) {
			return this.sequence.charAt(index);
			}
		
		@Override
		public int length() {
			return  this.sequence.length();
			}
		}
	
	/** mutated cdna */
	private static class MutedSequence extends DelegateCharSequence
		{
		private int pos0;
		private char newbase;
	
		MutedSequence(final CharSequence wild,int pos0,final char newbase)
			{
			super(wild);
			this.pos0 = pos0;
			this.newbase=newbase;
			if(this.pos0 > getDelegate().length()) throw new IndexOutOfBoundsException();
			}
		
		@Override
		public char charAt(int index)
			{
			return index==this.pos0?newbase:getDelegate().charAt(index);
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
			final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			final ContigNameConverter  ctgNameConverter = ContigNameConverter.fromOneDictionary(dict);
			LOG.info("loading genes from "+this.kgURI);
			in =IOUtils.openURIForBufferedReading(this.kgURI);
			final CharSplitter tab=CharSplitter.TAB;
			String line = null;
			int count_non_coding = 0;
			final Set<String> ignored=new HashSet<>();
			while((line=in.readLine())!=null)
				{
				line= line.trim();
				if(line.isEmpty()) continue;
				String tokens[]=tab.split(line);
				final KnownGene g=new KnownGene(tokens);
				if(g.isNonCoding()) {
					++count_non_coding;
					continue;
				}
				final String ctg = ctgNameConverter.apply(g.getContig());
				
				if(StringUtils.isBlank(ctg)) {
					if(ignored.add(g.getContig()))
						{
						LOG.warn("Unknown chromosome "+g.getContig()+" in dictionary");
						}
					continue;
					}
				g.setChrom(ctg);
				//use 1 based interval
				final Interval interval=new Interval( ctg, g.getTxStart()+1, g.getTxEnd() );
				List<KnownGene> lkg= this.knownGenes.get(interval);
				if(lkg==null) {
					lkg=new ArrayList<>(2);
					this.knownGenes.put(interval, lkg);
					}
				lkg.add(g);
				}
			CloserUtil.close(in);in=null;
			LOG.info("transcripts:"+knownGenes.size()+" discarded:" +count_non_coding +" non-coding transcript(s).");
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	
	
	private static int ID_GENERATOR=0;
	

	private static class CoverageInfo
		{
		int depth1=0;
		int depth2=0;
		int count_reads_having_both_variants = 0;
		int count_reads_having_no_variants = 0;
		int count_reads_having_variant1 = 0;
		int count_reads_having_variant2 = 0;
		}
	
	static private class AbstractContext implements Locatable{
		String contig;
		int genomicPosition1=0;
		String id=VCFConstants.EMPTY_ID_FIELD;
		Allele refAllele;
		Allele altAllele;
		int sorting_id;
		/* original vcf line */
		String vcfLine=null;
		
		@Override
		public String getContig() {
			return this.contig;
			}
		@Override
		public int getStart() {
			return genomicPosition1;
			}
		@Override
		public int getEnd() {
			return genomicPosition1;
			}
		}
	
	/** A mutation in a Transcript */
	static private class Variant extends AbstractContext
		{
		String transcriptName;
		byte strand;
		int position_in_cdna=-1;
		String wildCodon=null;
		String mutCodon=null;
		// genotypes carrying ALT may be empty but not null. 0 or 1
		byte genotypes[];
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
			
			this.genotypes = new byte[ctx.getNSamples()];
			for(int i=0;i< ctx.getNSamples();i++) {
				final Genotype gt = ctx.getGenotype(i);
				this.genotypes[i] = (byte)(gt.getAlleles().contains(allele)?1:0);
			}
			
			}
		/** position in codin 0 or 1 or 2 */
		int positionInCodon() { return  position_in_cdna%3;}
		/** start position in cdna of the codon overlapped by the mutation */
		int codonStart() { return this.position_in_cdna - this.positionInCodon();}
		
		
		Set<Integer> getSharedSampleIndexes(final Variant other) {
			if(this.genotypes.length==0) return Collections.emptySet();
			final Set<Integer> set = new TreeSet<>();
			for(int i=0;i< this.genotypes.length;i++) {
				if(this.genotypes[i]!=(byte)1) continue;
				if(other.genotypes[i]!=(byte)1) continue;
				set.add(i);
				}
			return set;
		}
		
		/** get specific info for this variant */
		private void _getInfo(final  Map<String,Object> info,int suffix) {
			//data about this
			info.put("POS"+suffix,String.valueOf(this.genomicPosition1));
			info.put("ID"+suffix,this.id);
			info.put("PosInCodon"+suffix,String.valueOf(this.positionInCodon()+1));
			info.put("Alt"+suffix,this.altAllele.getBaseString());
			info.put("Codon"+suffix,this.mutCodon);
			info.put("AA"+suffix,new ProteinCharSequence(this.mutCodon).getString());
			}
		
		public Map<String,Object> getInfo(final Variant other) {
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
			return info;
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
		VariantCodec() {
		}
		
		@Override
		public Variant decode(final DataInputStream dis) throws IOException {
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
			
			
			variant.vcfLine = readString(dis);
			
			final int n_genotypes = dis.readInt();
			variant.genotypes = new byte[n_genotypes];
			dis.readFully(variant.genotypes);
			
			return variant;
		}

		@Override
		public void encode(final DataOutputStream dos, final Variant v) throws IOException {
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
			
			writeString(dos, v.vcfLine);
			
			dos.writeInt(v.genotypes.length);
			dos.write(v.genotypes);
		}

		@Override
		public VariantCodec clone() {
			return new VariantCodec();
		}
		
		}
	/** equal range comparator */
	static private class VariantComparatorOne implements Comparator<Variant>
		{
		final ContigDictComparator contigDictComparator;
		VariantComparatorOne(final SAMSequenceDictionary dict) {
			this.contigDictComparator = new ContigDictComparator(dict);
			}
		@Override
		public int compare(final Variant o1,final  Variant o2) {
			int i= this.contigDictComparator.compare(o1.getContig(),o2.getContig());
			if(i!=0) return i;
			i= o1.transcriptName.compareTo(o2.transcriptName);
			return i;
			}
		}
	
	/** sorting collection comparator */
	static private class VariantComparatorTwo extends VariantComparatorOne
		{
		VariantComparatorTwo(final SAMSequenceDictionary dict) {
			super(dict);
			}
		@Override
		public int compare(final Variant o1,final  Variant o2) {
			int i= super.compare(o1, o2);
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
		int grantham_score = -1;
		// genotypes carrying ALT may be empty but not null. 0 or 1
		final Set<Integer> sampleIndexes=new HashSet<>();
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
			mutation.info = readString(dis);
			mutation.sorting_id = dis.readInt();
			mutation.filter = dis.readUTF();
			mutation.grantham_score = dis.readInt();
			mutation.vcfLine = readString(dis);
			
			int n_samples = dis.readInt();
			while(n_samples>0) {
				mutation.sampleIndexes.add(dis.readInt());
				n_samples--;
			}
			
			return mutation;
		}
	
		@Override
		public void encode(final DataOutputStream dos, final CombinedMutation v) throws IOException {
			dos.writeUTF(v.contig);
			dos.writeInt(v.genomicPosition1);
			dos.writeUTF(v.id);
			dos.writeUTF(v.refAllele.getBaseString());
			dos.writeUTF(v.altAllele.getBaseString());
			writeString(dos,v.info);
			dos.writeInt(v.sorting_id);
			dos.writeUTF(v.filter);
			dos.writeInt(v.grantham_score);
			writeString(dos, v.vcfLine);
			dos.writeInt(v.sampleIndexes.size());
			for(final int sampleId: v.sampleIndexes) {
				dos.writeInt(sampleId);
			}
		}
	
		@Override
		public MutationCodec clone() {
			return new MutationCodec();
		}
		
		}

static private class MutationComparatorOne implements Comparator<CombinedMutation>
	{
	final ContigDictComparator contigDictComparator;
	MutationComparatorOne(final SAMSequenceDictionary dict) {
		this.contigDictComparator = new ContigDictComparator(dict);
	}
	@Override
	public int compare(final CombinedMutation o1, final CombinedMutation o2) {
		int i= contigDictComparator.compare(o1.getContig(), o2.getContig());
		if(i!=0) return i;
		i= o1.genomicPosition1-o2.genomicPosition1;
		if(i!=0) return i;
		i =  o1.refAllele.compareTo(o2.refAllele);
		return i;
		}
	}

	
static private class MutationComparatorTwo extends MutationComparatorOne
	{
	MutationComparatorTwo(final SAMSequenceDictionary dict) {
		super(dict);
		}
	
	@Override
	public int compare(final CombinedMutation o1, final CombinedMutation o2) {
		int i= super.compare(o1, o2);
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
		return map.entrySet().stream().map(KV->KV.getKey()+"|"+KV.getValue()).collect(Collectors.joining("|"));
		}

	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		throw new IllegalStateException("should be never called");
		}
	
	
	
	@Override
	protected int doVcfToVcf(final String inputName,File saveAs)
		{
		BufferedReader bufferedReader = null;
		htsjdk.variant.variantcontext.writer.VariantContextWriter w=null;
		SortingCollection<CombinedMutation> mutations = null;
		CloseableIterator<Variant> varIter = null;
		CloseableIterator<CombinedMutation> mutIter = null;
		final Map<String,SamReader> sample2samReader = new HashMap<>();
		PrintWriter bedPeReport = null;
		try {
			bufferedReader = inputName==null?
						IOUtils.openStreamForBufferedReader(stdin()):
						IOUtils.openURIForBufferedReading(inputName);
			final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(bufferedReader);
			
			/* get VCF header */
			final VCFHeader header= cah.header;
			final List<String> sampleList = header.getSampleNamesInOrder();

			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(this.referencePath);
	        final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
	        
	        if(this.bamIn!=null)
	        	{
				final Set<String> sampleSet = new HashSet<>(sampleList);

	        	/** unroll and open bam file */
	        	for(final Path bamFile : IOUtils.unrollPaths(Collections.singletonList(this.bamIn.toString())))
		        	{
		        	LOG.info("opening BAM :"+this.bamIn);
		        	final SamReader samReader = SamReaderFactory.makeDefault().
		        			referenceSequence(this.referencePath).
		        			validationStringency(ValidationStringency.LENIENT).
		        			open(this.bamIn)
		        			;
		        	if(!samReader.hasIndex())
		        		{
		        		samReader.close();
		        		throw new IOException("Sam file is NOT indexed: "+bamFile);
		        		}
		        	final SAMFileHeader samHeader = samReader.getFileHeader();
		        	if(samHeader.getSequenceDictionary()==null || 
		        		!SequenceUtil.areSequenceDictionariesEqual(dict, samReader.getFileHeader().getSequenceDictionary())) {
		        		samReader.close(); 
		        		throw new JvarkitException.DictionariesAreNotTheSame(dict,samReader.getFileHeader().getSequenceDictionary());
		        		}
		        	/* get sample name */
		        	String sampleName=null;
		        	for(final SAMReadGroupRecord rg:samHeader.getReadGroups()) {
		        		if(rg.getSample()==null) continue;
		        		if(sampleName!=null && !sampleName.equals(rg.getSample())) {
		        			samReader.close();
			        		 throw new IOException(bamFile+" Contains two samples "+sampleName+" "+rg.getSample());
		        			}
		        		sampleName= rg.getSample();
		        		}
		        	if(sampleName==null) {
		        		samReader.close();
		        		LOG.warn("no sample in "+ bamFile);
		        		continue;
		        	}
		        	if(!sampleSet.contains(sampleName)) {
		        		samReader.close();
		        		LOG.warn("no sample "+sampleName+" in vcf. Ignoring "+bamFile);
		        		continue;
		        	}
		        	sample2samReader.put(sampleName, samReader);
		        	}
	        	}
	        
			loadKnownGenesFromUri();

			this.variants = SortingCollection.newInstance(Variant.class,
					new VariantCodec(),
					new VariantComparatorTwo(dict),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			this.variants.setDestructiveIteration(true);
			
			
			
	       ProgressFactory.Watcher<VariantContext> progress1= ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
		   String vcfLine=null;
	       while((vcfLine=bufferedReader.readLine())!=null)
				{
				final VariantContext ctx= progress1.apply(cah.codec.decode(vcfLine));
								
				/* discard non SNV variant */
				if(!ctx.isVariant() || ctx.isIndel())
					{
					continue;
					}
				
				/* find the overlapping genes : extend the interval of the variant to include the stop codon */
				final Collection<KnownGene> genes=  this.knownGenes.getOverlapping(
					new Interval(ctx.getContig(),
					Math.max(1,ctx.getStart()-3),
					ctx.getEnd()+3
					)).stream().flatMap(L->L.stream()).
					collect(Collectors.toList());
					
				
				final List<Allele> alternateAlleles =  ctx.getAlternateAlleles();
				
				/* loop over overlapping genes */
				for(final KnownGene kg:genes) {
					/* loop over available alleles */
					for(int allele_idx=0;allele_idx< alternateAlleles.size();++allele_idx) {
						final Allele alt = alternateAlleles.get(allele_idx);
						challenge(ctx,alt,kg,vcfLine);
						}
					}
				}
	       	progress1.close();
			this.variants.doneAdding();
			
			bedPeReport = this.bedPePath==null?
						new PrintWriter(new NullOuputStream()):
						IOUtils.openPathForPrintWriter(this.bedPePath)
						;
			
			mutations = SortingCollection.newInstance(CombinedMutation.class,
					new MutationCodec(),
					new MutationComparatorTwo(dict),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			mutations.setDestructiveIteration(true);
			
			
			final VCFFilterHeaderLine vcfFilterHeaderLine = new VCFFilterHeaderLine(
				"TwoHaplotypes",
				"(number of reads carrying both mutation) < (reads carrying variant 1 + reads carrying variant 2) "
				);
			
			varIter = this.variants.iterator();
			@SuppressWarnings("resource")
			EqualRangeIterator<Variant> eqVarIter=new EqualRangeIterator<>(varIter, new VariantComparatorOne(dict));
			
		    ProgressFactory.Watcher<Variant> progress2= ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			
			while(eqVarIter.hasNext())
				{
				final List<Variant> buffer= eqVarIter.next();
				
				if(buffer.size()<2) continue;
				
				for(int i=0;i+1< buffer.size();++i)
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
						// no sample share the two variants
						final Set<Integer> sharedSamplesIdx = v1.getSharedSampleIndexes(v2);
						if(sharedSamplesIdx.isEmpty() && !sampleList.isEmpty()) continue;
						
						final StringBuilder combinedCodon = new StringBuilder(v1.wildCodon);
						combinedCodon.setCharAt(v1.positionInCodon(), v1.mutCodon.charAt(v1.positionInCodon()));
						combinedCodon.setCharAt(v2.positionInCodon(), v2.mutCodon.charAt(v2.positionInCodon()));
						
						final String pwild = new ProteinCharSequence(v1.wildCodon).getString();
						final String p1 = new ProteinCharSequence(v1.mutCodon).getString();
						final String p2 = new ProteinCharSequence(v2.mutCodon).getString();
						final String pCombined = new ProteinCharSequence(combinedCodon).getString();
						final String combinedSO;
						final String combinedType;
						
						/* both AA are synonymous, while combined is not */
						if(!pCombined.equals(pwild) && p1.equals(pwild) && p2.equals(pwild)) {
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
						else if(!pCombined.equals(p1) && !pCombined.equals(p2) && !pCombined.equals(pwild))
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
							combinedType = null;
							combinedSO = null;
							}
						/** ok, there is something interesting here ,
						 * create two new Mutations carrying the
						 * two variants 
						 */
						if( combinedSO != null)
							{
							/** grantham score is max found combined vs (p1/p2/wild)*/
							int grantham_score= GranthamScore.score(pCombined.charAt(0), pwild.charAt(0));
							grantham_score = Math.max(grantham_score,GranthamScore.score(pCombined.charAt(0), p1.charAt(0)));
							grantham_score = Math.max(grantham_score,GranthamScore.score(pCombined.charAt(0), p2.charAt(0)));
							
							
							/** info that will be displayed in the vcf */
							final Map<String,Object> info1 = v1.getInfo(v2);
							final Map<String,Object> info2 = v2.getInfo(v1);
							//filter for this combined: default it fails the filter
							String filter = vcfFilterHeaderLine.getID();
							
							final Map<String,Object> combinedMap = new LinkedHashMap<>();
							combinedMap.put("CombinedCodon",combinedCodon);
							combinedMap.put("CombinedAA",pCombined);
							combinedMap.put("CombinedSO",combinedSO);
							combinedMap.put("CombinedType",combinedType);
							combinedMap.put("GranthamScore", grantham_score);
							
							info1.putAll(combinedMap);
							info2.putAll(combinedMap);
							
							final Map<String,CoverageInfo> sample2coverageInfo = new HashMap<>(sample2samReader.size());
							final int chromStart = Math.min(v1.genomicPosition1,v2.genomicPosition1);
							final int chromEnd = Math.max(v1.genomicPosition1,v2.genomicPosition1);

							
							
							/* get phasing info for each sample*/
							for(final String sampleName : sample2samReader.keySet() ) {
								final SamReader samReader = sample2samReader.get(sampleName);
								final CoverageInfo covInfo = new CoverageInfo();
								sample2coverageInfo.put(sampleName, covInfo);

								SAMRecordIterator  iter = null;
								try {
									iter = samReader.query(v1.contig,
											chromStart,
											chromEnd,
											false
											);
									while(iter.hasNext()) {
										final SAMRecord rec = iter.next();
										if(rec.getReadUnmappedFlag()) continue;
										if(rec.isSecondaryOrSupplementary()) continue;
										if(rec.getDuplicateReadFlag()) continue;
										if(rec.getReadFailsVendorQualityCheckFlag()) continue;
										
										// get DEPTh for variant 1
										if(rec.getAlignmentStart()<= v1.genomicPosition1 && v1.genomicPosition1<=rec.getAlignmentEnd()) {
											covInfo.depth1++;
											}
										// get DEPTh for variant 2
										if(rec.getAlignmentStart()<= v2.genomicPosition1 && v2.genomicPosition1<=rec.getAlignmentEnd()) {
											covInfo.depth2++;
											}
										
										if(rec.getAlignmentEnd()<chromEnd) continue;
										if(rec.getAlignmentStart()>chromStart) continue;
										final Cigar cigar  =  rec.getCigar();
										if(cigar==null) continue;
										final byte bases[] = rec.getReadBases();
										if(bases==null) continue;
										int refpos1=rec.getAlignmentStart();
										int readpos = 0;
										boolean found_variant1_on_this_read=false;
										boolean found_variant2_on_this_read=false;
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
														if(refpos1 == v1.genomicPosition1 && same(bases[readpos],v1.altAllele))
															{
															found_variant1_on_this_read = true;
															}
														else if(refpos1 == v2.genomicPosition1  && same(bases[readpos],v2.altAllele))
															{
															found_variant2_on_this_read = true;
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
									if( found_variant1_on_this_read && found_variant2_on_this_read) {
										covInfo.count_reads_having_both_variants++; }
									else if( !found_variant1_on_this_read && !found_variant2_on_this_read) {
										covInfo.count_reads_having_no_variants++; }
									else if( found_variant1_on_this_read) {
										covInfo.count_reads_having_variant1++;
									}else if( found_variant2_on_this_read) {
										covInfo.count_reads_having_variant2++;
									}
									
									
									}/* end of loop over reads */
									
								} finally {
									iter.close();
									iter=null;
									
								}
								
								info1.put("N_READS_BOTH_VARIANTS_"+sampleName, covInfo.count_reads_having_both_variants);
								info2.put("N_READS_BOTH_VARIANTS_"+sampleName, covInfo.count_reads_having_both_variants);
								info1.put("N_READS_NO_VARIANTS_"+sampleName, covInfo.count_reads_having_no_variants);
								info2.put("N_READS_NO_VARIANTS_"+sampleName, covInfo.count_reads_having_no_variants);
								info1.put("N_READS_TOTAL_"+sampleName,
										covInfo.count_reads_having_both_variants +
										covInfo.count_reads_having_no_variants +
										covInfo.count_reads_having_variant1+
										covInfo.count_reads_having_variant2
										);
								info2.put("N_READS_TOTAL_"+sampleName,
										covInfo.count_reads_having_both_variants +
										covInfo.count_reads_having_no_variants +
										covInfo.count_reads_having_variant1+
										covInfo.count_reads_having_variant2
										);
								//count for variant 1
								info1.put("N_READS_ONLY_1_"+sampleName, covInfo.count_reads_having_variant1);
								info1.put("N_READS_ONLY_2_"+sampleName, covInfo.count_reads_having_variant2);
								info1.put("DEPTH_1_"+sampleName, covInfo.depth1);
								//inverse previous count
								info2.put("N_READS_ONLY_1_"+sampleName, covInfo.count_reads_having_variant2);
								info2.put("N_READS_ONLY_2_"+sampleName, covInfo.count_reads_having_variant1);
								info2.put("DEPTH_2_"+sampleName, covInfo.depth2);
								
								/* number of reads with both variant is greater than
								 * reads carrying only one variant: reset the filter 
								 */
								if(2*covInfo.count_reads_having_both_variants>(covInfo.count_reads_having_variant1+covInfo.count_reads_having_variant2)) {
									/* reset filter */
									filter = VCFConstants.UNFILTERED;
									info1.put("FILTER_1_"+sampleName,".");
									info2.put("FILTER_2_"+sampleName,".");
									}
								else
									{
									info1.put("FILTER_1_"+sampleName,vcfFilterHeaderLine.getID());
									info2.put("FILTER_2_"+sampleName,vcfFilterHeaderLine.getID());
									}
							}/* end of loop over bams */
							
							
							final CombinedMutation m1 = new CombinedMutation();
							m1.contig = v1.contig;
							m1.genomicPosition1 = v1.genomicPosition1;
							m1.id = v1.id ;
							m1.refAllele = v1.refAllele;
							m1.altAllele = v1.altAllele;
							m1.vcfLine = v1.vcfLine;
							m1.info = mapToString(info1);
							m1.filter = filter;
							m1.grantham_score = grantham_score;
							m1.sampleIndexes.addAll(sharedSamplesIdx);
							m1.sorting_id = ID_GENERATOR++;
							mutations.add(m1);
							
							final CombinedMutation m2 = new CombinedMutation();
							m2.contig = v2.contig;
							m2.genomicPosition1 = v2.genomicPosition1;
							m2.id = v2.id ;
							m2.refAllele = v2.refAllele;
							m2.altAllele = v2.altAllele;
							m2.vcfLine = v2.vcfLine;
							m2.info = mapToString(info2);
							m2.filter = filter;
							m2.grantham_score = grantham_score;
							m2.sampleIndexes.addAll(sharedSamplesIdx);
							m2.sorting_id = ID_GENERATOR++;
							mutations.add(m2);
							
							bedPeReport.print(m1.contig);
							bedPeReport.print('\t');
							bedPeReport.print(m1.genomicPosition1-1);
							bedPeReport.print('\t');
							bedPeReport.print(m1.genomicPosition1);
							bedPeReport.print('\t');
							bedPeReport.print(m2.contig);
							bedPeReport.print('\t');
							bedPeReport.print(m2.genomicPosition1-1);
							bedPeReport.print('\t');
							bedPeReport.print(m2.genomicPosition1);
							bedPeReport.print('\t');
							bedPeReport.print(v1.transcriptName);//name
							bedPeReport.print('\t');
							bedPeReport.print(grantham_score==GranthamScore.getDefaultScore()?0:(int)((grantham_score/255.0)*1000.0));//score
							bedPeReport.print('\t');
							final KnownGene kg= this.knownGenes.getOverlapping(new Interval(v1.getContig(),v1.genomicPosition1-1,v1.genomicPosition1+1)).
									stream().
									flatMap(L->L.stream()).
									filter(P->P.getContig().equals(v1.contig) && P.getName().equals(v1.transcriptName)).
									findFirst().orElseThrow(IllegalStateException::new);
							
							bedPeReport.print(kg.isNegativeStrand()?"-":"+");//strand1
							bedPeReport.print('\t');
							bedPeReport.print(kg.isNegativeStrand()?"-":"+");//strand2
							bedPeReport.print('\t');
							if(sharedSamplesIdx.isEmpty()){
								bedPeReport.print('.');
								}
							else
								{
								bedPeReport.print(sharedSamplesIdx.stream().map(I->sampleList.get(I.intValue())).collect(Collectors.joining(";")));
								}
							bedPeReport.print('\t');
							bedPeReport.print(combinedSO);
							bedPeReport.print('\t');
							bedPeReport.print(String.join(":",pwild,p1,p2,pCombined));
							bedPeReport.println();
							}
						}
					}					
				}
			progress2.close();
			mutations.doneAdding();
			eqVarIter.close();eqVarIter=null;
			varIter.close();varIter=null;
			variants.cleanup();variants=null;
			
			bedPeReport.flush();
			bedPeReport.close();
			bedPeReport=null;
			
			
			final VCFHeader header2 = new VCFHeader(header);
			
			header2.addMetaDataLine(new VCFHeaderLine(getProgramName()+"AboutQUAL", "QUAL is filled with Grantham Score  http://www.ncbi.nlm.nih.gov/pubmed/4843792"));
			
			final StringBuilder infoDesc =new StringBuilder("Variant affected by two distinct mutation. Format is defined in the INFO column. ");
			
			
			final VCFInfoHeaderLine CodonVariantHeader = new VCFInfoHeaderLine(
					"CodonVariant",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					infoDesc.toString()
					);
			header2.addMetaDataLine(CodonVariantHeader);
			final VCFInfoHeaderLine CodonSampleHeader = new VCFInfoHeaderLine(
					"Samples",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					"Samples that could be affected"
					);
			header2.addMetaDataLine(CodonSampleHeader);
			
			
			JVarkitVersion.getInstance().addMetaData(this, header2);
			
			
			
			if(!sample2samReader.isEmpty())
				{
				header2.addMetaDataLine(vcfFilterHeaderLine);
				}
			
			w = super.openVariantContextWriter(saveAs);
			w.writeHeader(header2);
			
		    ProgressFactory.Watcher<CombinedMutation> progress3= ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			mutIter = mutations.iterator();
			EqualRangeIterator<CombinedMutation> eqRangeMutIter = new EqualRangeIterator<>(mutIter, new MutationComparatorOne(dict));
			
			while(eqRangeMutIter.hasNext())
				{
				final List<CombinedMutation> mBuffer = eqRangeMutIter.next();
				if(mBuffer.isEmpty()) break;				
				progress3.apply(mBuffer.get(0));
			
				//default grantham score used in QUAL
				int grantham_score = -1;
				//default filter fails
				String filter=vcfFilterHeaderLine.getID();
				final CombinedMutation first  = mBuffer.get(0);
				final Set<String> info = new HashSet<>();
				final VariantContext ctx  = cah.codec.decode(first.vcfLine);
				
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.chr(first.contig);
				vcb.start(first.genomicPosition1);
				vcb.stop(first.genomicPosition1 + first.refAllele.length()-1);
				if( !first.id.equals(VCFConstants.EMPTY_ID_FIELD)) vcb.id(first.id);
				
				for(final CombinedMutation m:mBuffer){	
					
					info.add(m.info);
					grantham_score=Math.max(grantham_score, m.grantham_score);
					if(VCFConstants.UNFILTERED.equals(m.filter)) {
						filter = null; //at least one SNP is ok one this line
					}
				}
				
				if(!sampleList.isEmpty())
					{
					vcb.attribute(CodonSampleHeader.getID(),
							new ArrayList<>(
								mBuffer.
								stream().
								flatMap(S->S.sampleIndexes.stream()).
								map(IDX->sampleList.get(IDX)).
								collect(Collectors.toSet())
								)
							);
					}
				vcb.unfiltered();
				if(filter!=null && !sample2samReader.isEmpty())
					{
					vcb.filter(filter);
					}
				else
					{
					vcb.passFilters();
					}
				
				vcb.attribute(CodonVariantHeader.getID(), new ArrayList<String>(info));
				
				if(grantham_score>0) {
					vcb.log10PError(grantham_score/-10.0);
					}
				else {
					vcb.log10PError(VariantContext.NO_LOG10_PERROR);
					}
				
				w.add(vcb.make());	
				}
			progress3.close();
			eqRangeMutIter.close();
			mutIter.close();
			mutations.cleanup();mutations=null;
			
			
			return RETURN_OK;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(mutIter);
			CloserUtil.close(varIter);
			CloserUtil.close(bedPeReport);
			if(this.variants!=null) this.variants.cleanup();
			if(mutations!=null) mutations.cleanup();
			this.variants=null;
			for(SamReader r: sample2samReader.values()) CloserUtil.close(r);
			CloserUtil.close(w);
			CloserUtil.close(bufferedReader);
			}
		
		
		}
	
	
	
	private CDNASequence getTranscript(final KnownGene kg) {
		CDNASequence mrna= this.kgId2transcriptCache.get(kg.getName());
		if(mrna!=null && mrna.knownGene.getContig().equals(kg.getContig())) return mrna;
		
		if(genomicSequence==null || !genomicSequence.getChrom().equals(kg.getContig()))
			{
			LOG.info("getting genomic Sequence for "+kg.getContig());
			genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, kg.getContig());
			}

		
		mrna = new CDNASequence(kg);
		this.kgId2transcriptCache.put(kg.getName(), mrna);
		return mrna;
	}
	
	private void challenge(
			final VariantContext ctx,
			final Allele allele,
			final KnownGene gene,
			final String vcfLine
			) throws IOException
		{
		if(allele.isSymbolic()) return;
		if(allele.isNoCall()) return;
		if(!allele.isCalled())return;
		if(allele.length()!=1) return;
		if(ctx.getReference().length()!=1) return;
		if(allele.equals(Allele.SPAN_DEL)) return;
		
		final CDNASequence wildRNA = this.getTranscript(gene);
		final int positionContext0 = ctx.getStart() -1;
		final Integer position_in_cnda = wildRNA.genomic2cnda.get(positionContext0);
		if(position_in_cnda==null) return;
				
		final Variant variant = new Variant(ctx,allele,gene);
		variant.sorting_id = ID_GENERATOR++;
		variant.position_in_cdna= position_in_cnda.intValue();
		char mutBase= Character.toUpperCase(allele.getDisplayString().charAt(0));
		if(gene.isNegativeStrand()) mutBase= AcidNucleics.complement(mutBase);
		
		
		final MutedSequence mutRNA = new MutedSequence(
				wildRNA,
				variant.position_in_cdna,
				mutBase
				);
		
		
		variant.wildCodon="";
		variant.mutCodon="";
		for(int i=0;i< 3;++i)
			{
			final int pos = variant.codonStart()+i;
			variant.wildCodon += (pos< wildRNA.length()?wildRNA.charAt(pos):'*');
			variant.mutCodon +=  (pos< mutRNA.length()?mutRNA.charAt(pos):'*');
			}
		variant.wildCodon = variant.wildCodon.toUpperCase();
		variant.mutCodon = variant.mutCodon.toUpperCase();
		variant.vcfLine  = vcfLine;
		
		
		
		if(variant.wildCodon.equals(variant.mutCodon)) {
			LOG.info("Uh??????? "+allele+" "+ctx.getContig()+":"+ctx.getStart()+" "+position_in_cnda+" "+wildRNA.length()+" "+variant.wildCodon+" "+variant.mutCodon);
			return;
			}
		
		this.variants.add(variant);
		
		}
	@Override
	public int doWork(final List<String> args) {
		if(StringUtils.isBlank(this.kgURI))
			{
			LOG.error("Undefined kgURI.");
			return -1;
			}
		return doVcfToVcf(args,this.outputFile);
		}
	
	public static void main(final String[] args)
		{
		new VCFCombineTwoSnvs().instanceMainWithExit(args);
		}
	}

