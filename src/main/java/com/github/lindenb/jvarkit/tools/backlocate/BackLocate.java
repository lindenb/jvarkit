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
package com.github.lindenb.jvarkit.tools.backlocate;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.AminoAcids;
import com.github.lindenb.jvarkit.util.bio.AminoAcids.AminoAcid;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Exon;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.PeptideSequence;
import com.github.lindenb.jvarkit.util.bio.structure.RNASequence;
import com.github.lindenb.jvarkit.util.bio.structure.RNASequenceFactory;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
/**
 BEGIN_DOC

## Example

mutation P->M at 1090 in NOTCH2

```
$  echo -e "NOTCH2\tP1090M" | java -jar dist/backlocate.jar -R hg19.fa --gtf ucsc.gtf
(...)
#User.Gene	AA1	petide.pos.1	AA2	knownGene.name	knownGene.strandknownGene.AA	index0.in.rna	codon	base.in.rna	chromosome	index0.in.genomic	exon
##uc001eik.3
NOTCH2	P	1090	M	uc001eik.3	NEGATIVE	P	3267	CCA	C	chr1	120480548	Exon 20
NOTCH2	P	1090	M	uc001eik.3	NEGATIVE	P	3268	CCA	C	chr1	120480547	Exon 20
NOTCH2	P	1090	M	uc001eik.3	NEGATIVE	P	3269	CCA	A	chr1	120480546	Exon 20
##uc001eil.3
NOTCH2	P	1090	M	uc001eil.3	NEGATIVE	P	3267	CCA	C	chr1	120480548	Exon 20
NOTCH2	P	1090	M	uc001eil.3	NEGATIVE	P	3268	CCA	C	chr1	120480547	Exon 20
NOTCH2	P	1090	M	uc001eil.3	NEGATIVE	P	3269	CCA	A	chr1	120480546	Exon 20
```

```
$ echo -e "NOTCH2\tPro1090M\tInteresting" | java -jar dist/backlocate.jar --gtf ucsc.gtf -R /path/to/human_g1k_v37.fasta | grep -v "##" | java -jar dist/prettytable.jar 

+------------+-----+--------------+-----+----------------+------------------+--------------+---------------+------------+----------------------+-------------+------------+-------------------+---------+-----------------+
| #User.Gene | AA1 | petide.pos.1 | AA2 | knownGene.name | knownGene.strand | knownGene.AA | index0.in.rna | wild.codon | potential.var.codons | base.in.rna | chromosome | index0.in.genomic | exon    | extra.user.data |
+------------+-----+--------------+-----+----------------+------------------+--------------+---------------+------------+----------------------+-------------+------------+-------------------+---------+-----------------+
| NOTCH2     | Pro | 1090         | Met | uc001eik.3     | -                | P            | 3267          | CCA        | .                    | C           | 1          | 120480548         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eik.3     | -                | P            | 3268          | CCA        | .                    | C           | 1          | 120480547         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eik.3     | -                | P            | 3269          | CCA        | .                    | A           | 1          | 120480546         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eil.3     | -                | P            | 3267          | CCA        | .                    | C           | 1          | 120480548         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eil.3     | -                | P            | 3268          | CCA        | .                    | C           | 1          | 120480547         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eil.3     | -                | P            | 3269          | CCA        | .                    | A           | 1          | 120480546         | Exon 20 | Interesting     |
+------------+-----+--------------+-----+----------------+------------------+--------------+---------------+------------+----------------------+-------------+------------+-------------------+---------+-----------------+
```


## See also

 * http://plindenbaum.blogspot.fr/2011/03/mapping-mutation-on-protein-to-genome.html
 * https://github.com/lindenb/jvarkit/issues/14
 * https://github.com/lindenb/jvarkit/issues/13


## History

 * 2019: move to GTF
 * 2019: add extra user data
 * 2017: Moved to jcommander
 * 2014: Moved to jvarkit
 * Nov 2014 : removed all the dependencies to SQL and DAS; use a local indexed genome
 * Aug 2015 : Added a new column "potention var codon" (as https://twitter.com/_ramrs/status/631123002005061633 ) , renamed "codon" to "wild codon"

## Cited in

backlocate was cited in:

 * CRISPR-STOP: gene silencing through base-editing-induced nonsense mutations. 2017 Nat Meth. [http://dx.doi.org/10.1038/nmeth.4327](http://dx.doi.org/10.1038/nmeth.4327).
 * "Differential 3' Processing of Specific Transcripts Expands Regulatory and Protein Diversity Across Neuronal Cell Types" Sasa Jereb, Hun-Way Hwang, Eric Van Otterloo, Eve-Ellen Govek, John J Fak, Yuan Yuan, Mary E Hatten, Robert B Darnell BioRxiv [https://www.biorxiv.org/content/biorxiv/early/2018/01/10/245886.full.pdf](https://www.biorxiv.org/content/biorxiv/early/2018/01/10/245886.full.pdf)
 
 END_DOC
 */
@Program(name="backlocate",
	description="Mapping a mutation on a protein back to the genome.",
	keywords={"vcf","annotation","prediction","protein"},
	biostars={15992,116366,425422},
	creationDate="20140619",
	modificationDate="20190820"
	)
public class BackLocate
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BackLocate.class).make();
	@Parameter(names={"-p","--printSeq"},description="print mRNA & protein sequences")
	private boolean printSequences = false;

	@Parameter(names={"-g","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfPath = null;

	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;

	

	private ReferenceSequenceFile referenceGenome = null;
	private GenomicSequence genomicContig = null;
	private final Map<String,List<Transcript>>  name2transcripts = new HashMap<>(100_000);
	/** get a genetic code from a chromosome name (either std or mitochondrial */
	private static GeneticCode getGeneticCodeByChromosome(final String chr)
		{
		if(chr.equalsIgnoreCase("chrM") || chr.equalsIgnoreCase("MT")) return GeneticCode.getMitochondrial();
		return GeneticCode.getStandard();
		}	


	private void backLocate(
		final PrintStream out,
		final Transcript transcript,
		final String geneName,
		final AminoAcid aa1,
		final AminoAcid aa2,
		int peptidePos1,
		final String extraUserData
		) throws IOException
		{
		final Set<String> messages = new LinkedHashSet<>();
		final GeneticCode geneticCode = getGeneticCodeByChromosome(transcript.getContig());
		final RNASequenceFactory rnaSequenceFactory=new RNASequenceFactory();
		
	        		
	        		
		if(this.genomicContig==null ||
		   !this.genomicContig.getChrom().equals(transcript.getContig())
	       )
	        	{
	        	final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.referenceGenome);
	        	final SAMSequenceRecord ssr =dict.getSequence(transcript.getContig());
	        	if(ssr==null) {
	        		LOG.warn(JvarkitException.ContigNotFoundInDictionary.getMessage(transcript.getContig(), dict));
	        		return;
	        		}
	        	this.genomicContig= new GenomicSequence(this.referenceGenome,transcript.getContig());
	        	}
    	rnaSequenceFactory.setContigToGenomicSequence(C->this.genomicContig);
        
		 final RNASequence wildRNA = rnaSequenceFactory.getCodingRNA(transcript);
		 final PeptideSequence<RNASequence> wildProt = PeptideSequence.of(wildRNA,geneticCode);
		 
	       		
	     final int peptideIndex0= peptidePos1-1;
         if(peptideIndex0 >=wildProt.length())
        	{
        	out.println("##index out of range for :"+transcript.getId()+" petide length="+wildProt.length());
        	return;
        	}
    
        if(wildProt.charAt(peptideIndex0)!= Character.toUpperCase(aa1.getOneLetterCode()))
        	{
        	messages.add("REF aminod acid ["+peptidePos1+"] is not the same ("+wildProt.charAt(peptideIndex0)+"/"+aa1+")");
        	}
       
        final int indexesInRNA[]=new int[]{
        	0+ peptideIndex0*3,
        	1+ peptideIndex0*3,
        	2+ peptideIndex0*3
        	};
        final String wildCodon=""
        		+ wildRNA.charAt(indexesInRNA[0])
        		+ wildRNA.charAt(indexesInRNA[1])
        		+ wildRNA.charAt(indexesInRNA[2])
        		;
        /* 2015 : adding possible mut codons */
        final Set<String> possibleAltCodons = new HashSet<>();
        final char bases[]=new char[]{'A','C','G','T'};
        for(int codon_pos=0;codon_pos<3;++codon_pos)
        	{
        	final StringBuilder sb=new StringBuilder(wildCodon);
        	for(char mutBase:bases)
        		{
        		sb.setCharAt(codon_pos, mutBase);
        		if(geneticCode.translate(sb.charAt(0), sb.charAt(1), sb.charAt(2))==Character.toUpperCase(aa2.getOneLetterCode()))
        			{
        			possibleAltCodons.add(sb.toString());
        			}
        		}
        	}
        
        for(int indexInRna: indexesInRNA)
        	{
        	out.print(geneName);
        	out.print('\t');
        	out.print(aa1.getThreeLettersCode());
        	out.print('\t');
        	out.print(peptidePos1);
        	out.print('\t');
        	out.print(aa2.getThreeLettersCode());
        	out.print('\t');
        	out.print(transcript.getGene().getGeneName());
        	out.print('\t');
        	out.print(transcript.getId());
        	out.print('\t');
        	out.print(transcript.getStrand());
        	out.print('\t');
        	out.print(wildProt.charAt(peptideIndex0));
        	out.print('\t');
        	out.print(indexInRna);
        	out.print('\t');
        	out.print(wildCodon);
        	out.print('\t');
        	if(possibleAltCodons.isEmpty())
        		{
        		out.print('.');
        		}
        	else
        		{
        		out.print(String.join("|", possibleAltCodons));
        		}
        	out.print('\t');
        	out.print(wildRNA.charAt(indexInRna));
        	out.print('\t');
        	out.print(transcript.getContig());
        	out.print('\t');
        	out.print(wildRNA.convertRnaIndex0ToGenomic0(indexInRna));
        	out.print('\t');
        	String exonName=null;
        	for(final Exon exon : transcript.getExons())
				{
				int genome=wildRNA.convertRnaIndex0ToGenomic0(indexInRna);
				if(exon.contains(genome))
					{
					exonName=exon.getName();
					break;
					}
				}
        	out.print(exonName);
        	if(this.printSequences)
        		{
        		String s=wildRNA.toString();
        		out.print('\t');
            	out.print(s.substring(0,indexInRna)+"["+s.charAt(indexInRna)+"]"+(indexInRna+1<s.length()?s.substring(indexInRna+1):""));
            	s=wildProt.toString();
            	out.print('\t');
            	out.print(s.substring(0,peptideIndex0)+"["+aa1+"/"+aa2+"/"+wildProt.charAt(peptideIndex0)+"]"+(peptideIndex0+1<s.length()?s.substring(peptideIndex0+1):""));
        		}
        	out.print('\t');
        	out.print(messages.isEmpty()?".":String.join(";",messages));
        	out.print('\t');
        	out.print(extraUserData);
        	out.println();
        	}
		}

	
	
	
	private void run(final PrintStream out,final BufferedReader in) throws IOException
		{
		final CharSplitter tab=CharSplitter.TAB;
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
			final String tokens[] = JvarkitException.TokenErrors.atLeast(2,tab.split(line,3));
			final String geneName=tokens[0].trim();
			if(StringUtils.isBlank(geneName)) throw new IOException("Bad line. No gene in "+geneName);
			final String mut=tokens[1].trim();
			int x0=0,x1=0;
			while(x1<mut.length() && Character.isLetter(mut.charAt(x1))) {
				++x1;
			}
			String substr = mut.substring(x0,x1);
			final AminoAcids.AminoAcid aa1 =  substr.length()==1?
					AminoAcids.getAminoAcidFromOneLetterCode(substr.charAt(0)):
					AminoAcids.getAminoAcidFromThreeLettersCode(substr)
					;
			if(aa1==null) throw new JvarkitException.UserError("Bad mutation "+mut+" (cannot parse left AA)");
			
			x0=x1;
			while(x1<mut.length() && Character.isDigit(mut.charAt(x1))) {
				++x1;
			}
			substr = mut.substring(x0,x1);
			if(!StringUtils.isInteger(substr)) throw new JvarkitException.UserError("Bad mutation "+mut+" (cannot parse position)");
			final int position1 = Integer.parseInt(substr);
			if(position1==0) throw new IOException("Bad position in protein ("+substr+") in "+line);
			
			substr = mut.substring(x1);
			final AminoAcids.AminoAcid aa2 =  substr.length()==1?
					AminoAcids.getAminoAcidFromOneLetterCode(substr.charAt(0)):
					AminoAcids.getAminoAcidFromThreeLettersCode(substr)
					;
			if(aa2==null) throw new JvarkitException.UserError("Bad mutation "+mut+" (cannot parse right AA)");
			
			final List<Transcript> transcripts= this.name2transcripts.get(geneName.toUpperCase());
			if(transcripts==null || transcripts.isEmpty())
				{
				LOG.warn("no transcript found for "+geneName);
				continue;
				}
			
			
			for(final Transcript transcript:transcripts)
				{
				backLocate(out,transcript, geneName, aa1, aa2, position1,tokens.length>2?tokens[2]:".");
				}
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		PrintStream out=null;
		BufferedReader in=null;
		try {
			this.referenceGenome = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			
			try(GtfReader gtfReader=new GtfReader(this.gtfPath)) {
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.referenceGenome);
				final ContigNameConverter contigNameConverter= ContigNameConverter.fromOneDictionary(dict);
				gtfReader.setContigNameConverter(contigNameConverter);
				gtfReader.getAllGenes().stream().
					filter(G->!StringUtils.isBlank(G.getGeneName())).
					flatMap(G->G.getTranscripts().stream()).
					filter(T->T.isCoding() && T.hasCDS()).
					forEach(
					T->{
						final String gn = T.getGene().getGeneName().toUpperCase();
						
						List<Transcript> L = this.name2transcripts.get(gn);
						if(L==null) {
							L=new ArrayList<>();
							this.name2transcripts.put(gn,L);
							}
						L.add(T);
						});
				
				}
			

			
			
			out = this.openPathOrStdoutAsPrintStream(this.outputFile);
			
			out.print("#User.Gene");
        	out.print('\t');
        	out.print("AA1");
        	out.print('\t');
        	out.print("petide.pos.1");
        	out.print('\t');
        	out.print("AA2");
        	out.print('\t');
        	out.print("transcript.name");
        	out.print('\t');
        	out.print("transcript.id");
        	out.print('\t');
        	out.print("transcript.strand");
        	out.print('\t');
        	out.print("transcript.AA");
        	out.print('\t');
        	out.print("index0.in.rna");
        	out.print('\t');
        	out.print("wild.codon");
        	out.print('\t');
        	out.print("potential.var.codons");
        	out.print('\t');
        	out.print("base.in.rna");
        	out.print('\t');
        	out.print("chromosome");
        	out.print('\t');
        	out.print("index0.in.genomic");
        	out.print('\t');
        	out.print("exon");
        	if(this.printSequences)
        		{
        		out.print('\t');
            	out.print("mRNA");
            	out.print('\t');
            	out.print("protein");
        		}
        	out.print('\t');
        	out.print("messages");
        	out.print('\t');
        	out.print("extra.user.data");
        	out.println();
			if(args.isEmpty())
				{
				in=super.openBufferedReader(null);
				this.run(out,in);
				CloserUtil.close(in);
				}
			else
				{
				for(final String filename:args)
					{
					in=super.openBufferedReader(filename);
					this.run(out,in);
					CloserUtil.close(in);
					}
				}
			return 0;
			}
		catch (final Exception e) {
			LOG.severe(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.referenceGenome);
			this.referenceGenome=null;
			CloserUtil.close(out);
			CloserUtil.close(in);
			}	
		}
		
	public static void main(final String[] args)
		{
		new BackLocate().instanceMainWithExit(args);
		}
	}
