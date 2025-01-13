/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.fastqsw;

import java.io.File;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.Predicate;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.template.SequenceView;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.jcommander.OnePassFastqLauncher;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.LogicalOp;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.fastq.FastqRecord;


/**
BEGIN_DOC

## Example:

```
$ java -jar dist/fastqsw.jar -s 'TCGTGACCTCTCAGGTTAGACCTACAAAATCTTCGATCAAT' \
	--min-identity 0.9 --save-align jeter.txt --min-length 40 --pair-logical OR \
	src/test/resources/S1.R1.fq.gz src/test/resources/S1.R2.fq.gz


(...)

@RF11_137_644_0:0:0_2:0:0_23/1
GCTCCCTCGTGACCTCTCAGGTTAGACCTACAAATCTTCGATCAATAGCATTGCGACTTGTTTCATTCTC
+
2222222222222222222222222222222222222222222222222222222222222222222222
@RF11_137_644_0:0:0_2:0:0_23/2
GTACATTTCACCAGATGCAGAAGCATTCAGTAAATACATGCTGTCAAAGTCTCCAGAAGATATTGGACCA
+
2222222222222222222222222222222222222222222222222222222222222222222222

$ tail -6 jeter.txt 
#RF11_137_644_0:0:0_2:0:0_23/1 distance:0.46 score:189.0 similarity:0.54 length:41 identicals:40 similars:40 identity:0.975609756097561
    1                                      41
  7 TCGTGACCTCTCAGGTTAGACCTAC-AAATCTTCGATCAAT 46
    ||||||||||||||||||||||||| |||||||||||||||
  1 TCGTGACCTCTCAGGTTAGACCTACAAAATCTTCGATCAAT 41
```



END_DOC
*/
//https://github.com/broadinstitute/Drop-seq/blob/72920f2b3824c03ccfd0abf6fb57a179e01962a8/src/java/org/broadinstitute/dropseqrna/sbarro/utils/ConsensusSequenceFactory.java

@Program(name="fastqsw",
	keywords={"fastq","align","sw"},
	description="align fasta sequences vs fastq",
	modificationDate="20220207",
	creationDate="20220207",
	biostars=9534472
	)
public class FastqSW extends OnePassFastqLauncher {
	private static final Logger LOG = Logger.build(FastqSW.class).make();

	
	@Parameter(names={"--pairwise-aligner-type"},description="One item from org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType")
	private PairwiseSequenceAlignerType pairwiseSequenceAlignerType = PairwiseSequenceAlignerType.LOCAL;
	@Parameter(names={"--disable-reverse-complement"},description="disable search in reverse-complement.")
	boolean disable_reverse_complement = false;
	@Parameter(names={"-f","--queries"},description="Fasta file containing queries")
	private File fastaQueries = null;
	@Parameter(names={"-s","--string-queries"},description="Query as a DNA string")
	private List<String> user_string_queries = new ArrayList<>();

	@Parameter(names={"-open-penalty"},description="Gap open penalty.")
	private int open_penalty = 10;
	@Parameter(names={"-extension-penalty"},description="Gap extension penalty.")
	private int extension_penalty = 1;
	@Parameter(names={"--min-distance"},description="Minimum distance")
	private double min_distance = 0.0;
	@Parameter(names={"--min-similarity"},description="Minimum similarity. "+ FractionConverter.OPT_DESC,splitter=NoSplitter.class,converter=FractionConverter.class)
	private double min_similarity = 0.0;
	@Parameter(names={"--min-length"},description="Minimum length")
	private int min_length = 0;
	@Parameter(names={"--max-gap"},description="Maximum number of gaps. -1 to ignore")
	private int max_gap = -1;

	
	@Parameter(names={"--min-identity"},description="Minimum identity. "+ FractionConverter.OPT_DESC,splitter=NoSplitter.class,converter=FractionConverter.class)
	private double min_identity = 9.0;
	@Parameter(names={"--disable-count-gap"},description="Do not count gaps when calculating percentage of identity.")
	private boolean disable_count_gap = false;
	@Parameter(names={"--save-align"},description="Save alignments in this file (mostly for debugging)")
	private Path saveAlignFile= null;
	@Parameter(names={"--pair-logical"},description="Logical choice for pair of reads.")
	private LogicalOp pairLogicalOp = LogicalOp.OR;

	@Parameter(names={"--inverse"},description="Inverse logic, discard matching reads.")
	private boolean inverse_logic = false;
	
	private static final DNACompoundSet DNA_COMPOUND_SET = AmbiguityDNACompoundSet.getDNACompoundSet();
	private final GapPenalty gapPenalty = new SimpleGapPenalty();
	private PrintStream saveAlignStream = null;
	
	
	private final List<DNASequence> queries = new ArrayList<>();
	/* Returns Nuc 4.4 matrix by Lowe Both of the nucleotide sequences to align can contain ambiguous nucleotides */
	private final SubstitutionMatrix<NucleotideCompound> substitutionMatrix =  SubstitutionMatrixHelper.getNuc4_4();

	
	private static class MyAlign {
		PairwiseSequenceAligner<Sequence<NucleotideCompound>, NucleotideCompound> sw;
		SequencePair<Sequence<NucleotideCompound>,NucleotideCompound> pair;
	}
	
	
	private void saveAlign(final FastqRecord rec, MyAlign a) {
		if(this.saveAlignStream==null) return;
		this.saveAlignStream.println("#"+rec.getReadName()+
				" distance:"+ a.sw.getDistance()+ 
				" score:"+ a.sw.getScore()+ 
				" similarity:"+ a.sw.getSimilarity()+ 
				" length:"+a.pair.getLength()+
				" identicals:"+a.pair.getNumIdenticals()+
				" similars:"+a.pair.getNumSimilars()+
				" identity:"+a.pair.getPercentageOfIdentity(!this.disable_count_gap)
				);
		this.saveAlignStream.println(a.pair.toString(150));
		}
	
	private Optional<MyAlign> 
		align(final FastqRecord rec) {
	final boolean count_gaps = !this.disable_count_gap;
	final DNASequence recseq;
	try {
		recseq = new DNASequence(rec.getReadString(),DNA_COMPOUND_SET);
		}
	catch(final CompoundNotFoundException err) {
		throw new IllegalArgumentException("Cannot get reverse complement for "+rec.getReadName(),err);
		}
	final SequenceView<NucleotideCompound> rcseq = disable_reverse_complement?
			null:
			recseq.getReverseComplement();
	
	for(final DNASequence query : queries) {
		for(int side=0;side<2;++side) {
			if(side==1 && disable_reverse_complement) continue;
			final PairwiseSequenceAligner<Sequence<NucleotideCompound>, NucleotideCompound> sw  = Alignments.getPairwiseAligner(
					(side==0?recseq:rcseq),
					query,
					this.pairwiseSequenceAlignerType,
					this.gapPenalty,
					this.substitutionMatrix
					);
			final SequencePair<Sequence<NucleotideCompound>,NucleotideCompound> pair = sw.getPair();
			if(pair ==null) continue;
			if(pair.getPercentageOfIdentity(count_gaps) < this.min_identity) continue;
			if(pair.getLength() < this.min_length) continue;
			if(sw.getSimilarity() < this.min_similarity) continue;
			if(sw.getDistance() < this.min_distance) continue;
			if(this.max_gap>=0 && 
					pair.getAlignedSequences().stream().
						anyMatch(AS->AS.getNumGaps()>  max_gap)) {
				continue;
				}
			
			
			
			final MyAlign a  = new MyAlign();
			a.sw = sw;
			a.pair = pair;
			return Optional.of(a);
			}
		}
	return Optional.empty();
	}
	
	@Override
	protected Predicate<FastqRecord> createPredicateForFastqRecord() {
		return (REC)->{
			 final Optional<MyAlign> opt = align(REC);
			 if(opt.isPresent()) {
				saveAlign(REC,opt.get());
				return !inverse_logic;
			 	}
			 return inverse_logic;
			};
		}
	@Override
	protected Predicate<FastqRecordPair> createPredicateForFastqRecordPair() {
		return (RECS)->{
			 final Optional<MyAlign> optR1 = align(RECS.getFirstInPair());
			 final Optional<MyAlign> optR2;
			 //no need to align R2 if 'OR' mode and R1 is aligned
			 if(pairLogicalOp.equals(LogicalOp.OR) && optR1.isPresent()) {
				optR2 = Optional.empty(); 
			 	}
			//no need to align R2 if 'AND' mode and R1 is not aligned
			 else if(pairLogicalOp.equals(LogicalOp.AND) && !optR1.isPresent()) {
				return inverse_logic;
			 	}
			 else
			 	{
				optR2 =  align(RECS.getSecondInPair());
			 	}

			 if(pairLogicalOp.test(optR1.isPresent(),optR2.isPresent())) {
				if(optR1.isPresent())saveAlign(RECS.getFirstInPair(),optR1.get());
				if(optR2.isPresent())saveAlign(RECS.getSecondInPair(),optR2.get());
				return !inverse_logic;
			 	}
			 return inverse_logic;
			};
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeFastq() {
		try {
			if(fastaQueries!=null) {
				FastaReaderHelper.readFastaDNASequence(fastaQueries, false).values().
					stream().
					forEach(S->this.queries.add(S));
			}
			
			for(final String s:this.user_string_queries) {
				this.queries.add(new DNASequence(s,DNA_COMPOUND_SET));
			}
			
			if(this.queries.isEmpty()) {
				LOG.error("No query defined");
				return -1;
				}
			
			if(this.saveAlignFile!=null) {
				this.saveAlignStream = openPathOrStdoutAsPrintStream(this.saveAlignFile);
				}
			
			
			this.gapPenalty.setOpenPenalty(this.open_penalty);
			this.gapPenalty.setExtensionPenalty(this.extension_penalty);
			}
		catch(final Throwable err) {
			LOG.error(err);
			}
		return super.beforeFastq();
		}
	
@Override
protected void afterFastq() {
	if(this.saveAlignStream !=null) {
		this.saveAlignStream.flush();
		this.saveAlignStream.close();
		}
	super.afterFastq();
	}	
	
	public static void main(final String[] args) {
		new FastqSW().instanceMainWithExit(args);
		}
	}
