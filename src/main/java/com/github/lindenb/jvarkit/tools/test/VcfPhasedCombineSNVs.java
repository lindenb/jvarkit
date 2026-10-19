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
package com.github.lindenb.jvarkit.tools.test;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.AminoAcids;
import com.github.lindenb.jvarkit.bio.AminoAcids.AminoAcid;
import com.github.lindenb.jvarkit.bio.GeneticCode;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.grantham.GranthamScore;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.RNASequence;
import com.github.lindenb.jvarkit.util.bio.structure.RNASequenceFactory;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

@Program(name="vcfphasedcombinesnvs",
description="Detect Mutations than are the consequences of two distinct variants using phased vcf. Idea from @SolenaLS and then @AntoineRimbert",
keywords={"vcf","annotation","prediction","protein","mnv"},
modificationDate="20160215",
creationDate="20160219",
jvarkit_amalgamion = true,
menu="VCF Manipulation"
)
public class VcfPhasedCombineSNVs extends Launcher {
	private static final Logger LOG = Logger.of(VcfPhasedCombineSNVs.class);

	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true )
	private Path referencePath;
	@Parameter(names={"--gtf"},description=GtfReader.OPT_DESC,required=true )
	private Path gtfPath;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@Parameter(names={"--debug"},description="debug mode")
	protected boolean debug_flag = false;
	@Parameter(names={"--bed","--regions-file"},description="limit to the genes overlapping that bed file")
	protected Path regions_file_bed = null;
	@Parameter(names={"--phase","--phased"},description="Only consider phased variants.")
	protected boolean phased_only=false;

	/** reference genome */
	private ReferenceSequenceFile indexedFastaSequenceFile=null;

	
	private static class Codon {
		int query_pos0_in_cdna;
		int shifted_position0;
		String mrna_codon;
		AminoAcid aminoAcid;
		@Override
		public String toString() {
			return aminoAcid.getName()+"("+mrna_codon+") pos0:"+query_pos0_in_cdna+"/shifted0:"+shifted_position0;
			}
		}
	
	/** extended sequence where I can */
	private static class TranslatableRNASequence extends DelegateCharSequence
		{
		TranslatableRNASequence(final CharSequence wild)
			{
			super(wild);
			}
		
		int inFramePosition(int index0) {
			return index0 - (index0%3);
			}
		
		Codon getCodonAt(int index0) {
			final Codon codon=new Codon();
			codon.query_pos0_in_cdna  = index0;
			int shifted0 = inFramePosition(index0);
			codon.shifted_position0 = shifted0;
			codon.mrna_codon = new StringBuilder(3)
					.append(charAt(shifted0+0))
					.append(charAt(shifted0+1))
					.append(charAt(shifted0+2))
					.toString();
			final GeneticCode gencode = GeneticCode.getStandard();
			final char c =  gencode.translate(
					codon.mrna_codon.charAt(0),
					codon.mrna_codon.charAt(1),
					codon.mrna_codon.charAt(2)
					);
			if(gencode.isStop(c)) {
				codon.aminoAcid = AminoAcids.Stop;
				}
			else
				{
				codon.aminoAcid = AminoAcids.getAminoAcidFromOneLetterCode(c);
				}
			if(codon.aminoAcid==null) throw new IllegalStateException("cannot find AA for "+c+" "+codon.mrna_codon);
			return codon;
			}
		}

	
	/** mutated cdna */
	private static class MutedRNASequence extends TranslatableRNASequence
		{
		private final int pos0;
		private final char newbase;
	
		MutedRNASequence(final CharSequence wild,int pos0,final char newbase)
			{
			super(wild);
			this.pos0 = pos0;
			this.newbase=newbase;
			if(this.pos0 > getDelegate().length()) throw new IndexOutOfBoundsException();
			}
		
		@Override
		public char charAt(final int index)
			{
			return index==this.pos0 ? this.newbase : getDelegate().charAt(index);
			}
	
		}
	
	private static boolean isAltAllele(Allele alt) {
		if(alt.isReference() || alt.isSymbolic() || alt.equals(Allele.SPAN_DEL) || alt.isNoCall()) return false;
		if(!AcidNucleics.isATGC(alt)) return false;
		return true;
	}
	
	private static String toString(VariantContext ctx) {
		return  new StringBuilder(ctx.getContig())
				.append(":")
				.append(ctx.getStart())
				.append("/")
				.append(ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/")))
				.toString();
		}
	
	private static String toString(VariantContext ctx,Allele alt) {
		return  new StringBuilder(ctx.getContig())
				.append(":")
				.append(ctx.getStart())
				.append(":")
				.append(ctx.getReference().getDisplayString())
				.append("/")
				.append(alt.getDisplayString())
				.toString();
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Path vcfPath = Paths.get(oneAndOnlyOneFile(args));
			this.indexedFastaSequenceFile= ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referencePath);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			final String buildName = SequenceDictionaryUtils.getBuildName(dict).orElse(null);
			final Hyperlink hyperlink = Hyperlink.compile(dict);
			ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(dict);
			final List<Gene> genes;
			final Predicate<Locatable> acceptInterval;
			if(regions_file_bed!=null) {
				try(BedLineReader blr=new BedLineReader(this.regions_file_bed)) {
					blr.setContigNameConverter(ctgConverter);
					final IntervalTreeMap<Boolean> treeMap= blr.toIntervalTreeMap(B->true);
					acceptInterval = L->treeMap.containsOverlapping(L);
					}
				}
			else
				{
				acceptInterval = L->true;
				}
			
			try (GtfReader gtfReader=new GtfReader(this.gtfPath)) {
				gtfReader.setContigNameConverter(ctgConverter);
				genes = gtfReader.getAllGenes()
					.stream()
					.filter(G->acceptInterval.test(G))
					.filter(G->G.getTranscripts().stream() .anyMatch(T->T.isCoding() && T.hasCodonStartDefined() && T.hasCodonStopDefined()))
					.sorted(new ContigDictComparator(dict).createLocatableComparator())
					.collect(Collectors.toList());
				}
			if(this.debug_flag) {
				LOG.error("number of genes: "+ genes.size());
				LOG.error("number of transcripts: "+ genes.stream().mapToLong(G->G.getTranscripts().stream().count()).sum());
				}
			
			if(genes.isEmpty()) {
				LOG.error("no gene was found");
				return -1;
				}
			
			final RNASequenceFactory rnaSequenceFactory = new RNASequenceFactory(this.indexedFastaSequenceFile);
			
			try(VCFFileReader vcfReader0 = new VCFFileReader(vcfPath,true)) {
				final VCFHeader header = vcfReader0.getFileHeader();
				if(!header.hasGenotypingData()) {
					LOG.error("no genotype in "+vcfPath);
					return -1;
					}
				final int buffer_size= 100_000;
				try(BufferedVCFReader vcfReader=new BufferedVCFReader(vcfReader0, buffer_size)) {
					try(JsonWriter js = new JsonWriter(super.openPathOrStdoutAsPrintWriter(this.outputFile))) {
						js.beginArray();
						for(Gene gene: genes) {
							if(this.debug_flag) {
								LOG.debug(gene.getGeneName()+" "+gene.getContig()+":"+gene.getStart()+"-"+gene.getEnd());
								}
							for(Transcript transcript: gene.getTranscripts()) {
								if(!transcript.isCoding()) continue;
								if(!transcript.hasCodonStartDefined()) continue;
								if(!transcript.hasCodonStopDefined()) continue;
								if(this.debug_flag) {
									LOG.debug(transcript.getId()+" "+transcript.getContig()+":"+transcript.getStart()+"-"+transcript.getEnd());
									}
								RNASequence codingRNA = null;
								for(String sn: header.getGenotypeSamples()) {
									if(this.debug_flag) {
										LOG.debug("sample "+sn);
										}
									final List<VariantContext> variants = new ArrayList<>(transcript.getLengthOnReference());
									boolean non_diploid_flag=false;
									try(CloseableIterator<VariantContext> iter =  vcfReader.query(transcript)) {
										while(iter.hasNext()) {
											final VariantContext ctx = iter.next();
											
											final Genotype gt = ctx.getGenotype(sn);
											if(gt.getPloidy()!=2) {
												if(this.debug_flag) {
													LOG.debug("found non-diploid");
													}
												non_diploid_flag =  true;
												break;
												}
											
											variants.add(new VariantContextBuilder(ctx)
												.attributes(Collections.emptyMap())
												.genotypes(Collections.singletonList(gt))
												.make());
											}
										}
									if(non_diploid_flag) {
										if(this.debug_flag) {
											LOG.debug("skipping "+sn+" "+transcript.getId()+" because found non-diploid");
											}
										continue;
										}
									// at least two genotypes with ALT allele
									if(variants.stream()
											.flatMap(G->G.getGenotypes().stream())
											.filter(G->G.getAlleles().stream().anyMatch(A->isAltAllele(A)))
											.count()<2L) {
										if(this.debug_flag) {
											LOG.info("not enough variants in the interval for "+sn);
											}
										continue;
										}
									if(this.debug_flag) {
										LOG.info("number of variants : "+variants.size()+" for "+sn);
										}
									if(codingRNA==null) {
										codingRNA = rnaSequenceFactory.getCodingRNA(transcript);
										}
									/** loop over each phase for variant 1*/
									for(int phase1_idx=0;phase1_idx<2;++phase1_idx) {
										if(this.debug_flag) {
											LOG.info("phase1=" + phase1_idx +" for "+sn);
											}
										/** loop over each variant in that transcript */
										for(int x1=0;x1+1< variants.size();++x1) {
											final VariantContext vc1= variants.get(x1);
											if(vc1.getReference().length()!=1) {
												if(this.debug_flag) {
													LOG.info("vc1 bad ref  " + toString(vc1));
													}
												continue;
												}
											if(this.debug_flag) {
												LOG.info("vc1 = " + toString(vc1));
												}
											final Genotype g1 = vc1.getGenotype(sn);
											
											
											
											if(this.phased_only && !g1.isPhased()) {
												if(this.debug_flag) {
													LOG.info("phase1=" + phase1_idx +" for "+sn);
													}
												continue;
												}
											
											if(phase1_idx>=g1.getPloidy()) {
												if(this.debug_flag) {
													LOG.info("type=" +g1.getType().name() +" for "+sn);
													}
												continue;
												}
											
											final Allele alt1 = g1.getAllele(phase1_idx);
											if(!isAltAllele(alt1)) {
												if(this.debug_flag) {
													LOG.info("bad alt allele1 "+alt1+" vs "+vc1.getReference()+" "+alt1.isNonReference()+" "+alt1.isReference()+" "+alt1.isNonRefAllele());
													}
												continue;
												}
											if(alt1.length()!=1) {
												if(this.debug_flag) {
													LOG.info("bad length(alt) allele1 "+alt1+" "+alt1.length());
													}
												continue;
												}
											/* get position in mRNA*/
											final OptionalInt pos1 =  codingRNA.convertGenomic0ToRnaIndex0(vc1.getStart()-1);
											if(!pos1.isPresent()) {
												if(this.debug_flag) {
													LOG.debug("pos1 not in cDNA "+vc1.getStart());
													}
												continue;
												}
											if(this.debug_flag) {
												LOG.debug("cdnapos1 :  "+pos1);
												}
											char base1 = alt1.getBaseString().charAt(0);
											if(transcript.isNegativeStrand()) base1 = AcidNucleics.complement(base1);
											final MutedRNASequence rna1 = new MutedRNASequence(codingRNA, pos1.getAsInt(), base1);
											// position shifted to the first base of the codong
											final Codon wildCodon = new TranslatableRNASequence(codingRNA).getCodonAt(pos1.getAsInt()-1);
											final Codon codon1 = rna1.getCodonAt(pos1.getAsInt()-1);
											if(this.debug_flag) {
												LOG.debug("codon1 :  "+codon1);
												}
											
											if(!wildCodon.aminoAcid.equals(codon1.aminoAcid)) {
												if(this.debug_flag) {
													LOG.debug("non synonymous codon 1 "+codon1+" vs wild "+wildCodon);
													}
												continue;//we want synonymous
												}
											
											if(this.debug_flag) {
												LOG.debug("ok: synonymous codon 1 "+codon1+" vs wild "+wildCodon);
												}
											/** loop over each phase of second variant */
											for(int phase2_idx=0;phase2_idx<2;++phase2_idx) {
												/* both phases must be the same if this.phased_only*/
												if(this.phased_only && phase1_idx!=phase2_idx) {
													continue;
													}
												/** loop over each variant in that transcript */
												for(int x2=x1+1; x2< variants.size();++x2) {
													final VariantContext vc2= variants.get(x2);
													
													if(this.debug_flag) {
														LOG.info("vc2 = " + toString(vc2));
														}
													
													if(vc2.getReference().length()!=1) {
														if(this.debug_flag) {
															LOG.info("vc2 bad REF length : " + toString(vc2));
															}
														continue;
														}
													
													final Genotype g2 = vc2.getGenotype(sn);
													if(this.phased_only && !g2.isPhased()) {
														/* end of PHASE, let's stop here */
														if(this.debug_flag) {
															LOG.info("gt2 is not phased");
															}
														break;
														}
													

													
													if(phase2_idx >= g2.getPloidy()) {
														if(this.debug_flag) {
															LOG.info("gt2 has bad ploidy "+g2.getType().name());
															}
														continue;
														}
													
													if(vc1.getStart()>= vc2.getStart()) {
														if(this.debug_flag) {
															LOG.info("Hum.. strange order ?"+vc1+" / "+vc2);
															}
														continue;
														}
													final Allele alt2 = g2.getAllele(phase2_idx);
													if(!isAltAllele(alt2)) {
														if(this.debug_flag) {
															LOG.info("bad alt allele2 "+alt2);
															}
														continue;
														}
													if(alt2.length()!=1) {
														if(this.debug_flag) {
															LOG.debug("bad length for allele2 "+alt2);
															}
														continue;
														}
													/* get position in mRNA*/
													final OptionalInt pos2 =  codingRNA.convertGenomic0ToRnaIndex0(vc2.getStart()-1);
													if(!pos2.isPresent()) {
														if(this.debug_flag) {
															LOG.debug("pos2 is not in cdna");
															}
														continue;
														}
													
													if(this.debug_flag) {
														LOG.debug("cdnapos2 :  "+pos2);
														}
													
													char base2 = alt2.getBaseString().charAt(0);
													if(transcript.isNegativeStrand()) base2 = AcidNucleics.complement(base2);
													
													final MutedRNASequence rna2 = new MutedRNASequence(codingRNA, pos2.getAsInt(), base2);
													/* pos1 and pos2 must be in the same codon ! */
													
													
													final Codon codon2 = rna2.getCodonAt(pos2.getAsInt()-1);
													
													/* both mutations must fall in the same codon */
													if( codon1.shifted_position0!=codon2.shifted_position0 ) {
														if(this.debug_flag) {
															LOG.debug(" shifted_pos1!=shifted_pos2 "+codon1.shifted_position0+" "+codon2.shifted_position0);
															}
														continue;
														}
													
													if(!wildCodon.aminoAcid.equals(codon2.aminoAcid)) {
														if(this.debug_flag) {
															LOG.debug("non synonymous codon2  "+codon2+" vs wild "+wildCodon);
															}
														continue;//we want synonymous
														}
													if(this.debug_flag) {
														LOG.debug("ok : synonymous codon2  "+codon2+" vs wild "+wildCodon);
														}
		
													final MutedRNASequence rna_combined = new MutedRNASequence(rna1, pos2.getAsInt(), base2);
													final Codon codonCombined = rna_combined.getCodonAt(pos2.getAsInt()-1);
													
													/* even combined, the AA is the same e.g: leucine */
													if(wildCodon.aminoAcid.equals(codonCombined.aminoAcid)) {
														if(this.debug_flag) {
															LOG.debug("no change, even combined is same combine:" +
																		codonCombined+" vs wild "+
																	wildCodon+" codon1:"+codon1+" codon2:"+codon2+
																	" gt1:"+toString(vc1,alt1)+ " gt2:"+toString(vc2,alt2));
															}
														continue;//we don't want synonymous
														}
													
													
													js.beginObject();
													js.name("contig"); js.value(transcript.getContig());
													js.name("transcript_start"); js.value(transcript.getStart());
													js.name("transcript_end"); js.value(transcript.getEnd());
													js.name("strand"); js.value(transcript.getStrand());
													js.name("transcript_id"); js.value(transcript.getId());
													js.name("gene_id"); js.value(transcript.getGene().getId());
													js.name("gene_name"); js.value(transcript.getGene().getGeneName());
													js.name("pos_cdna"); js.value(pos2.getAsInt()+1);
													js.name("wild_codon"); js.value(wildCodon.mrna_codon);
													js.name("wild_aa"); js.value(wildCodon.aminoAcid.getThreeLettersCode());
													js.name("sample"); js.value(sn);
													if(this.phased_only) {
														js.name("phase"); js.value(phase1_idx+1);
														}
													else
														{
														// label if all variants are phased in the interval
														boolean all_phased=true;
														for(int x3=x1;x3<=x2;++x3) {
															final Genotype gt3 = variants.get(x3).getGenotype(sn);
															if(!gt3.isPhased()) {
																all_phased = false;
																break;
																}
															}
														js.name("continuous_phase"); js.value(all_phased);
														}
													js.name("variants");
													js.beginArray();
													for(int side=0;side< 2;++side) {
														final VariantContext vc =(side==0?vc1:vc2);
														js.beginObject();
														if(buildName!=null) {
															js.name("build"); js.value(buildName);
															}
														js.name("contig"); js.value(vc.getContig());
														js.name("pos"); js.value(vc.getStart());
														js.name("ref"); js.value(vc.getReference().getBaseString());
														if(vc.hasID()) {
															js.name("id"); js.value(vc.getID());
															}
														if(vc.isFiltered()) {
															js.name("filter"); 
															js.beginArray();
															for(String fmt: vc.getFilters()) {
																js.value(fmt);
																}
															js.endArray();
															}
														final Genotype gt =(side==0?g1:g2);
														if(gt.hasGQ()) {
															js.name("gq"); js.value(gt.getGQ());
															}
														if(gt.hasDP()) {
															js.name("dp"); js.value(gt.getDP());
															}
														if(gt.hasAD()) {
															js.name("ad");
															js.beginArray();
															for(int p: gt.getAD())
																{
																js.value(p);
																}
															js.endArray();
															}
														if(gt.hasPL()) {
															js.name("pl");
															js.beginArray();
															for(int p: gt.getPL())
																{
																js.value(p);
																}
															js.endArray();
															}
														if(gt.isFiltered()) {
															js.name("filter"); 
															js.value(gt.getFilters());
															}
														js.name("type"); js.value(gt.getType().name());
	
														js.name("alleles");
														js.beginArray();
														for(Allele a:gt.getAlleles()) {
															js.value(a.getBaseString());
															}
														js.endArray();
														
														js.name("alt");
														js.value((side==0?alt1:alt2).getBaseString());
														
														js.name("phase"); js.value(1+(side==0?phase1_idx:phase2_idx));
														
														js.name("pos_cdna");
														js.value((side==0?codon1:codon2).query_pos0_in_cdna+1);
														js.name("pos_codon");
														js.value((side==0?codon1:codon2).shifted_position0);
														js.name("mrna_codon");
														js.value((side==0?codon1:codon2).mrna_codon);
														js.name("amino_acid");
														js.value((side==0?codon1:codon2).aminoAcid.getThreeLettersCode());
														
														final Optional<String> url= hyperlink.apply(vc);
														if(url.isPresent()) {
															js.name("url");
															js.value(url.get());
															}
														
														js.endObject();
														}
													js.endArray();
													js.name("distance");
													js.value(vc2.getStart()-vc1.getEnd());
													js.name("combined_mrna_codon");
													js.value(codonCombined.mrna_codon);
													js.name("combined_amino_acid");
													js.value(codonCombined.aminoAcid.getName());
													js.name("combined_grantham_score");
													js.value(GranthamScore.score(
															wildCodon.aminoAcid.getOneLetterCode(),
															codonCombined.aminoAcid.getOneLetterCode())
															);
	
													
													js.endObject();
													} // end variant2
												}//end phase2
											}
										}
									}
								}
							}
						js.endArray();
						js.flush();
						}// end of output writer
					}
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			if(indexedFastaSequenceFile!=null) {
				try {
					indexedFastaSequenceFile.close();
					}
				catch(IOException err) {}
				}
			}
		}
	
	public static void main(final String[] args) {
		new VcfPhasedCombineSNVs().instanceMainWithExit(args);
	}

}
