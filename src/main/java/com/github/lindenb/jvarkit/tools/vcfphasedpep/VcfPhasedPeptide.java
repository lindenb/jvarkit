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
package com.github.lindenb.jvarkit.tools.vcfphasedpep;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.GeneticCode;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.PeptideSequence;
import com.github.lindenb.jvarkit.util.bio.structure.RNASequence;
import com.github.lindenb.jvarkit.util.bio.structure.RNASequenceFactory;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
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
/**
BEGIN_DOC

# Input/Output

Input is an indexed VCF file

Output is a JSON file containing the predictions for each transcript

# Example

```
$ java -jar dist/jvarkit.jar vcfphasedpep -R ref.fa --gtf input.gtf.gz input.vcf.gz > out.json

(...)
```
END_DOC

 */

@Program(name="vcfphasedpep",
description="Find altered peptides in phased VCF",
keywords= {"protein","phased","vcf"},
creationDate="20261026",
modificationDate="20261027",
jvarkit_amalgamion = true,
menu="VCF Manipulation"
)
public class VcfPhasedPeptide extends Launcher {
	private static final Logger LOG = Logger.of(VcfPhasedPeptide.class);
	
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
	@Parameter(names={"--regex-samples"},description="only process samples matching that regex.")
	protected String regexSamples=null;
	@Parameter(names={"--both-phases"},description="both phases must be altered")
	protected boolean both_phases_must_be_altered=false;
	@Parameter(names={"--accept-filtered"},description="Accept Filtered Genotype (default is to skip FILTERed genotypes")
	protected boolean accept_filtered_genotype_flag =false;

	/** reference genome */
	private ReferenceSequenceFile indexedFastaSequenceFile=null;

	
	private static class Base {
		char base;
		}
	
	private static class BaseIsRef extends Base {
		int genomic0;
		
		}
	private static class BaseIsAlt extends Base {
		@SuppressWarnings("unused")
		VariantContext variant;
		}        

	
	private static class SequenceOfBases extends AbstractCharSequence {
		private final List<Base> bases;
		
		SequenceOfBases(final List<Base> bases) {
			this.bases = bases;
			}
		@Override
		public int length() {
			return this.bases.size();
			}
		@Override
		public char charAt(int index) {
			return bases.get(index).base;
			}
		}

	
	private static boolean isAltAllele(final Allele alt) {
		if(alt.isReference() || alt.isSymbolic() || alt.equals(Allele.SPAN_DEL) || alt.isNoCall()) return false;
		if(!AcidNucleics.isATGC(alt)) return false;
		return true;
	}
	
	private static String toString(final VariantContext ctx) {
		return  new StringBuilder(ctx.getContig())
				.append(":")
				.append(ctx.getStart())
				.append("/")
				.append(ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/")))
				.toString();
		}
	
	
	
	
	private String peptideSeqUntilStop(final PeptideSequence<CharSequence>  pep) {
		final StringBuilder sb=new StringBuilder(pep.length());
		final GeneticCode code = GeneticCode.getStandard();
		for(int i=0;i< pep.length();i++ ) {
			char c= pep.charAt(i);
			if(code.isStop(c)) break;
			sb.append(c);
			}
		return sb.toString();
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Path vcfPath = Paths.get(oneAndOnlyOneFile(args));
			this.indexedFastaSequenceFile= ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referencePath);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			//final String buildName = SequenceDictionaryUtils.getBuildName(dict).orElse(null);
			//final Hyperlink hyperlink = Hyperlink.compile(dict);
			ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(dict);
			GenomicSequence genomicSequence = null;
			final List<Gene> genes;
			final Predicate<Locatable> acceptInterval;
			if(regions_file_bed!=null) {
				try(BedLineReader blr=BedLineReader.openPath(this.regions_file_bed)) {
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
							/* loop over transcripts */
							for(Transcript transcript: gene.getTranscripts()) {
								if(!transcript.isCoding()) continue;
								if(!transcript.hasCodonStartDefined()) continue;
								if(!transcript.hasCodonStopDefined()) continue;
								if(this.debug_flag) {
									LOG.debug(transcript.getId()+" "+transcript.getContig()+":"+transcript.getStart()+"-"+transcript.getEnd());
									}
								
								//position used to detect if a variant is in UTR ?
								final OptionalInt leftmost_cds = transcript.getAllCds().stream().mapToInt(X->X.getStart()).min();
								final OptionalInt righttmost_cds = transcript.getAllCds().stream().mapToInt(X->X.getEnd()).max();

								
								RNASequence codingRNA = null;
								PeptideSequence<CharSequence> wildProt = null;
								String wildPepUntilStop = null;
								
								/* loop over samples */
								for(String sn: header.getGenotypeSamples()) {
									if(this.debug_flag) {
										LOG.debug("sample "+sn);
										}
									/* skip that sample ?*/
									if(this.regexSamples!=null && !sn.matches(regexSamples)) {
										if(this.debug_flag) {
											LOG.debug("skipping sample "+sn);
											}
										continue;
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
									// at least one genotype with ALT allele
									if(variants.stream()
											.flatMap(G->G.getGenotypes().stream())
											.flatMap(L->L.getAlleles().stream())
											.noneMatch(A->isAltAllele(A))
											) {
										if(this.debug_flag) {
											LOG.info("not ALT variants in the interval for "+sn);
											}
										continue;
										}
									if(this.debug_flag) {
										LOG.info("number of variants : "+variants.size()+" for "+sn);
										}
									if(codingRNA==null) {
										codingRNA = rnaSequenceFactory.getCodingRNA(transcript);
										wildProt = PeptideSequence.of(codingRNA);
										wildPepUntilStop = peptideSeqUntilStop(wildProt);
										}
									
									
									final CharSequence[] mutated_peptides_each_phase=new CharSequence[] {null,null};
									/** loop over each phase for variant 1*/
									for(int phase1_idx=0;phase1_idx<2;++phase1_idx) {
										if(this.debug_flag) {
											LOG.info("phase1=" + phase1_idx +" for "+sn);
											}
										int prev_variant_index=0;
										boolean phase_interrupted=false;
										final List<Base> bases= new ArrayList<Base>();
										int transcript_index1=transcript.getStart();
										
										while(transcript_index1<=transcript.getEnd() && !phase_interrupted) {
											
											
											boolean got_base_flag=false;
											for(int i=prev_variant_index;i< variants.size() && !got_base_flag;i++) {
												final VariantContext ctx= variants.get(i);
												if(!accept_filtered_genotype_flag && ctx.isFiltered()) {
													continue;
													}
												final Genotype gt = ctx.getGenotype(sn);
												if(!gt.isPhased() && transcript.isPositiveStrand()) {
													//could be a variant in UTR ?
													if(leftmost_cds.getAsInt()<= transcript_index1) {
														phase_interrupted=true;
														if(this.debug_flag) {
															LOG.info("phase is interrupted for "+sn+" at "+toString(ctx)+" for "+transcript.getId());
															}
														break;
														}
													}
												if(ctx.getStart() > transcript_index1) {
													break;
													}
												if(ctx.getStart()==transcript_index1) {
													
													// unphased downstream ?
													if(transcript.isNegativeStrand()) {
														//test variants after that one
														for(int j=i+1;j< variants.size();j++) {
															final VariantContext ctx_downstream= variants.get(j);
															if(ctx_downstream.getEnd()>=righttmost_cds.getAsInt()) break;
															if(!this.accept_filtered_genotype_flag &&  ctx_downstream.isFiltered()) {
																continue;
																}
															final Genotype gt_downstream = ctx_downstream.getGenotype(sn);
															if(!gt_downstream.isPhased() && ctx_downstream.getEnd() <= righttmost_cds.getAsInt()) {
																phase_interrupted=true;
																if(this.debug_flag) {
																	LOG.info("phase is interrupted downstream for "+sn+" at "+toString(ctx_downstream)+" for "+transcript.getId());
																	}
																break;
																}
															}
														}
													
													
													prev_variant_index = i;
													final Allele alt_allele = gt.getAllele(phase1_idx);
													// test wether position falls in transcript
													final OptionalInt rnaIndex0 =  codingRNA.convertGenomic0ToRnaIndex0(ctx.getStart()-1);
													if(isAltAllele(alt_allele) && rnaIndex0.isPresent()) {
														final Allele refAllele = ctx.getReference();
														if(refAllele.length()<0) throw new IllegalArgumentException();
														final String dna = alt_allele.getBaseString();
														
														if(this.debug_flag) {
															LOG.info("got variant in cdna for sample "+sn+" at position "+transcript_index1);
															}
														
														for(int z=0;z< dna.length();z++) {
															final BaseIsAlt altbase=new BaseIsAlt();
															altbase.variant = ctx;
															altbase.base = dna.charAt(z);
															
															bases.add(altbase);
															got_base_flag=true;
															}
														transcript_index1+=refAllele.length();
														}
													break;
													}
												}
											//no variant found overlaping position, use reference sequence
											if(!got_base_flag) {
												if(codingRNA.convertGenomic0ToRnaIndex0( transcript_index1-1).isPresent()) {
													if(genomicSequence==null || !genomicSequence.hasName(transcript.getContig())) {
														genomicSequence = new GenomicSequence(indexedFastaSequenceFile, transcript.getContig());
														}
													final BaseIsRef refbase = new BaseIsRef();
													refbase.genomic0 = transcript_index1-1;
													refbase.base = genomicSequence.charAt(refbase.genomic0);
													bases.add(refbase);
													}
												++transcript_index1;
												}
											}
										
										if(phase_interrupted) {
											if(this.debug_flag) {
												LOG.info("phase is interrupted for "+sn);
												}
											continue;
											}
										
										if(transcript.isNegativeStrand()) {
											if(this.debug_flag) {
												LOG.info("reverse complementing "+transcript.getId()+" "+transcript.getStrand());
												}
											Collections.reverse(bases);
											for(int i=0;i< bases.size();i++) {
												final Base b = bases.get(i);
												b.base = AcidNucleics.complement(b.base);
												}
											}
										
										
										
										final SequenceOfBases mutatedRNA = new SequenceOfBases(bases);
										final PeptideSequence<CharSequence> mutatedProt = PeptideSequence.of(mutatedRNA);
										
										final String pepMut = peptideSeqUntilStop(mutatedProt);
										if(wildPepUntilStop.length()!= pepMut.length()) {
											if(this.debug_flag) {
												LOG.info("got a bad guy ! "+transcript.getId()+" "+transcript.getStrand()+" "+
														wildPepUntilStop.length()+" vs "+ pepMut.length()
														);
												}
											mutated_peptides_each_phase[phase1_idx] = pepMut;
											}
										else
											{
											if(this.both_phases_must_be_altered) {
												break;
												}
											}
										} // phase
									if(mutated_peptides_each_phase[0]==null && mutated_peptides_each_phase[1]==null) {

										if(this.debug_flag) {
											LOG.info("both sequences are ok ! "+transcript.getId()+" "+transcript.getStrand());
											}
										continue;
										}
									
									if(this.both_phases_must_be_altered && (mutated_peptides_each_phase[0]==null || mutated_peptides_each_phase[1]==null)) {
										if(this.debug_flag) {
											LOG.info("both sequences altered=false false "+transcript.getId()+" "+transcript.getStrand());
											}
										continue;
										}
										
									
									
									if(this.debug_flag) {
										LOG.info("Pushing JSON");
										}
									js.beginObject();
									js.name("contig"); js.value(transcript.getContig());
									js.name("transcript_start"); js.value(transcript.getStart());
									js.name("transcript_end"); js.value(transcript.getEnd());
									js.name("strand"); js.value(String.valueOf(transcript.getStrand()));
									js.name("transcript_id"); js.value(transcript.getId());
									js.name("gene_id"); js.value(transcript.getGene().getId());
									js.name("gene_name"); js.value(transcript.getGene().getGeneName());
									js.name("wild_peptide"); js.value(wildPepUntilStop);
									js.name("wild_peptide_length"); js.value(wildPepUntilStop.length());
									js.name("sample"); js.value(sn);
									js.name("both_phases_are_mutated"); js.value(mutated_peptides_each_phase[0]!=null && mutated_peptides_each_phase[1]!=null);
									js.name("phases"); 
									js.beginArray();
									for(int phase1_idx=0;phase1_idx<2;++phase1_idx) {
										if(mutated_peptides_each_phase[phase1_idx]==null) continue;
										js.beginObject();
										js.name("phase"); js.value(phase1_idx+1);
										js.name("mutated_peptide"); js.value(mutated_peptides_each_phase[phase1_idx].toString());
										js.name("mutated_peptide_length"); js.value(mutated_peptides_each_phase[phase1_idx].length());
										js.endObject();
										}
									
									js.endArray();
									js.endObject();
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
	new VcfPhasedPeptide().instanceMain(args);
	}
}