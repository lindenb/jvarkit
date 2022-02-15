package com.github.lindenb.jvarkit.tools.fastqsw;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
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
import com.github.lindenb.jvarkit.fastq.FastqPairedReaderFactory;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriterFactory;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassFastqLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.CloseableIterator;


/**
BEGIN_DOC

## Input


END_DOC
*/
//https://github.com/broadinstitute/Drop-seq/blob/72920f2b3824c03ccfd0abf6fb57a179e01962a8/src/java/org/broadinstitute/dropseqrna/sbarro/utils/ConsensusSequenceFactory.java

@Program(name="fastqsw",
	keywords={"fastq","align","sw"},
	description="align fasta sequences vs fastq",
	modificationDate="20220207",
	creationDate="20220207",
	generate_doc=false
	)
public class FastqSW extends OnePassFastqLauncher {
	private static final Logger LOG = Logger.build(FastqSW.class).make();

	
	@Parameter(names={"--pairwise-aligner-type"},description="One item from org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType")
	private PairwiseSequenceAlignerType pairwiseSequenceAlignerType = PairwiseSequenceAlignerType.LOCAL;
	@Parameter(names={"--disable-reverse-complement"},description="disable search in reverse-complement.")
	boolean disable_reverse_complement = false;
	@Parameter(names={"--queries"},description="Fasta file containing queries")
	private File fastaQueries = null;

	
	
	private static final DNACompoundSet DNA_COMPOUND_SET = AmbiguityDNACompoundSet.getDNACompoundSet();

	private final List<DNASequence> queries = new ArrayList<>();
	private GapPenalty gapPenalty = new GapPenalty();
	private SubstitutionMatrix<NucleotideCompound> substitutionMatrix;
	
	private boolean match(final FastqRecord rec) throws CompoundNotFoundException {
	final DNASequence recseq = new DNASequence(rec.getReadString(),DNA_COMPOUND_SET);
	final SequenceView<NucleotideCompound> rcseq = disable_reverse_complement?null:recseq.getReverseComplement();
	
	for(final DNASequence query : queries) {
		for(int side=0;side<2;++side) {
			if(side==1 && disable_reverse_complement) continue;
			final PairwiseSequenceAligner<Sequence<NucleotideCompound>, NucleotideCompound> sw  = Alignments.getPairwiseAligner(
					(side==0?recseq:rcseq),
					query,
					pairwiseSequenceAlignerType,
					gapPenalty,
					substitutionMatrix
					);
			if(sw.getSimilarity() < 0.0 ) continue;
			}
		}
	return false;
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
			
			if(this.queries.isEmpty()) {
				LOG.error("No query defined");
				return -1;
				}
			}
		catch(final Throwable err) {
			LOG.error(err);
			}
		return super.beforeFastq();
		}
	public static void main(String[] args) {
		new FastqSW().instanceMainWithExit(args);
		}
	}
