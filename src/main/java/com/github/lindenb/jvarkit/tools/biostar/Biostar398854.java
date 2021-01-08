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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;
/**
 BEGIN_DOC
 ## Example
 
 ```
 java -jar ${JVARKIT_DIST}/biostar398854.jar \
 	--gtf input.gtf.gz \
 	-R ref.fasta input.vcf > out.fasta
 
 ```
 
 END_DOC
 */
@Program(name="biostar398854",
	description="Extract every CDS sequences from a VCF file",
	keywords= {"vcf","gtf"},
	biostars=398854,
	creationDate="20190916",
	modificationDate="20190916"
	)
public class Biostar398854 extends Launcher {
	private static final Logger LOG = Logger.build(Biostar398854.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-gtf","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfIn = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	
	private ReferenceSequenceFile referenceSequenceFile = null;;
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter out = null;
		try {
			this.referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.referenceSequenceFile);
			
			
			try(VCFReader in = VCFReaderFactory.makeDefault().open(Paths.get(oneAndOnlyOneFile(args)),true)) {
				out = super.openPathOrStdoutAsPrintWriter(this.outputFile);
				final PrintWriter final_out= out;
				final List<String> samples = in.getHeader().getSampleNamesInOrder();

				
				try(GtfReader gtfReader= new GtfReader(this.gtfIn)) {
					final SAMSequenceDictionary dict2 = in.getHeader().getSequenceDictionary();
					if(dict2!=null) SequenceUtil.assertSequenceDictionariesEqual(dict, dict2);
					gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					gtfReader.getAllGenes().
						stream().
						flatMap(G->G.getTranscripts().stream()).
						filter(T->T.hasCDS()).
						forEach(transcript->{
							final List<VariantContext> variants1 = in.query(transcript).
									stream().
									filter(V->V.isVariant() && AcidNucleics.isATGCN(V.getReference()) &&  V.getAlternateAlleles().stream().anyMatch(A->AcidNucleics.isATGCN(A))).
									collect(Collectors.toCollection(ArrayList::new));
							
							if(variants1.isEmpty()) return;
							
							final int positions1[]=transcript.getAllCds().
									stream().
									flatMapToInt(CDS->IntStream.rangeClosed(CDS.getStart(), CDS.getEnd())).
									toArray();

							
							final List<VariantContext> variants = variants1.
									stream().
									filter(V->{
										int insert = Arrays.binarySearch(positions1, V.getStart());
										return insert>=0 && insert < positions1.length;
										}).
									collect(Collectors.toCollection(ArrayList::new));
							if(variants.isEmpty()) return;
							
							final ReferenceSequence refSeq = this.referenceSequenceFile.getSubsequenceAt(transcript.getContig(),transcript.getStart(), transcript.getEnd());
							
							
							
							for(int nSample=0;nSample<=/* yes <= */ samples.size();nSample++)
								{
								final String fastaName= transcript.getId()+" "+transcript.getGene().getId()+" "+transcript.getGene().getGeneName()+" "+(nSample< samples.size()?samples.get(nSample):"ALL")+
											" " + transcript.getContig()+":"+transcript.getStart()+"-"+transcript.getEnd()+"("+transcript.getStrand()+")";
								final StringBuilder sb = new StringBuilder();
								int array_index=0;
								
								while(array_index< positions1.length) {
									final int x1 = positions1[array_index];
									
									final int refseqidx0 = x1-transcript.getStart();
									
									char refbase = (refseqidx0<0 || refseqidx0>=refSeq.length()?'N':(char)refSeq.getBases()[refseqidx0]);
									final int x1_final = x1;
									String base= String.valueOf(refbase);
									final VariantContext ctx=variants.stream().filter(V->V.getStart()==x1_final).findFirst().orElse(null);
									Allele alt = ctx==null?null:ctx.getAlternateAlleles().stream().filter(A->AcidNucleics.isATGCN(A)).findFirst().orElse(null);
									if(ctx!=null && nSample< samples.size()) {
										final Genotype gt = ctx.getGenotype(nSample);
										alt = gt.getAlleles().stream().filter(A->!A.isReference() &&! A.isNoCall() && AcidNucleics.isATGCN(A)).findFirst().orElse(null);
										}
									
									if(alt!=null)
										{
										base = alt.getBaseString().toUpperCase();
										int i=0;
										while(i< ctx.getReference().length() && array_index<positions1.length)
											{
											array_index++;
											i++;
											}
										}
									else
										{
										base = base.toLowerCase();
										array_index++;
										}
									sb.append(base);									
									}
								
								String fastaSeq  = transcript.isNegativeStrand()?AcidNucleics.reverseComplement(sb.toString()):sb.toString();
								final_out.print(">");
								final_out.println(fastaName);
								final_out.println(fastaSeq);
								}
							}
						);
					}
				out.flush();
				out.close();
				out=null;
				}
			return 0;
		} catch (final Throwable e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(this.referenceSequenceFile);
			}
		}	
		
public static void main(final String[] args) {
	new Biostar398854().instanceMainWithExit(args);
	}
}
