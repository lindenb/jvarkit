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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
/**

BEGIN_DOC



### Contig conversion

chromosome 23 is converted to X or chrX, chromosome 24 is converted to Y or chrY, chromosome 25 is ignored, chromosome 26 is converted to chrM or MT.


### Example


```
$ java -jar dist/bim2vcf.jar -R human_g1k_v37.fasta input.bim 

##fileformat=VCFv4.2
##INFO=<ID=MORGAN,Number=1,Type=Float,Description="Centimorgan">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variation type">
##contig=<ID=1,length=249250621,assembly=human_g1k_v37>
(...)
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
1       12      rs73422 C       .       .       .       MORGAN=0.5224;SVTYPE=NOVARIATION
1       13      rs30315 G       A       .       .       MORGAN=0.530874;SVTYPE=SNV
1       14      rs14325 C       T       .       .       MORGAN=0.532596;SVTYPE=SNV
1       15      rs31319 A       G       .       .       MORGAN=0.532682;SVTYPE=SNV
1       16      rs954   C       T       .       .       MORGAN=0.537655;SVTYPE=SNV
1       17      rs62034 G       A       .       .       MORGAN=0.548645;SVTYPE=SNV
1       18      rs25996 A       G       .       .       MORGAN=0.575595;SVTYPE=SNV
1       19      rs12117 G       A       .       .       MORGAN=0.582608;SVTYPE=SNV
(...)

```


END_DOC
*/
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;


@Program(name="bim2vcf",
	description="convert a .bim to a .vcf . For @FlorianeS44",
	keywords= {"bim","vcf"},
	modificationDate="20200116"
)
public class BimToVcf extends Launcher
	{
	private static final Logger LOG = Logger.build(BimToVcf.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;

	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();
	
	public BimToVcf() {
		}
	
	@Override
	public int doWork(final List<String> args) {
		VariantContextWriter w=null;
		BufferedReader r=null;
		ReferenceSequenceFile faidx=null;
		GenomicSequence genomic = null;
		long number_non_snvs = 0L;
		try {
			
			faidx = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			
			final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(faidx);
			
			this.writingVariants.dictionary(dict);
			
			r = super.openBufferedReader(oneFileOrNull(args));
				
			final Set<VCFHeaderLine> headerLines = new HashSet<>();
			final VCFInfoHeaderLine morgan= new VCFInfoHeaderLine("MORGAN", 1, VCFHeaderLineType.Float,"Centimorgan");
			final VCFInfoHeaderLine svtype= new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String,"Variation type");
			VCFStandardHeaderLines.addStandardInfoLines(headerLines, false, "");
			//super.addMetaData(headerLines);
			headerLines.add(morgan);
			headerLines.add(svtype);

			final List<String> genotypeSampleNames = Collections.emptyList();
			final VCFHeader header=new VCFHeader(headerLines, genotypeSampleNames);
			header.setSequenceDictionary(dict);
			w = this.writingVariants.open(this.outputFile);
			w.writeHeader(header);
			final CharSplitter tab=CharSplitter.TAB;
			String line;
			final Predicate<String> iupacATGC = S->AcidNucleics.isATGC(S);
			
			final List<String> convertCtg = Arrays.asList(
					"23", "X",
					"23", "chrX",
					"24", "Y",
					"24", "chrY",
					"26", "MT",
					"26", "chrM"
					);
			while((line=r.readLine())!=null) {
				final String tokens[]=tab.split(line);
				if(tokens.length!=6) {
					throw new JvarkitException.TokenErrors(6, tokens);
				}
				Allele a1=null;
				Allele a2=null;
				Allele ref=null;
				final String contig=tokens[0];
			
				SAMSequenceRecord ssr= null;
				ssr = dict.getSequence(contig);
				
				for(int i=0;ssr==null && i+1< convertCtg.size() ;i+=2) {
					if(!convertCtg.get(i).equals(contig)) continue;
					ssr = dict.getSequence(convertCtg.get(i+1));
					}
					
				if(ssr==null && contig.equals("25")){
					LOG.warn("ignoring "+line);
					++number_non_snvs;
					continue;
					}
				if(ssr==null){
					throw new JvarkitException.ContigNotFoundInDictionary(contig,dict);
					}
				if(genomic==null || !ssr.getSequenceName().equals(genomic.getChrom())) {
					genomic=new GenomicSequence(faidx, ssr.getSequenceName());
					}
				int pos1 = Integer.parseInt(tokens[3]);
				if(tokens[4].equals("0")) tokens[4]=tokens[5];
				if(tokens[5].equals("0")) tokens[5]=tokens[4];
				
				final VariantContextBuilder vcb = new VariantContextBuilder();
				vcb.chr(ssr.getSequenceName());
				vcb.attribute(morgan.getID(), Float.parseFloat(tokens[2]));
				
				if(iupacATGC.test(tokens[4]) && iupacATGC.test(tokens[5]))
					{
					if(tokens[4].length()!=1 || tokens[5].length()!=-1) {
						LOG.warn("Skipping "+line+" because I only handle SNVs.");
						number_non_snvs++;
						continue;
						}
					
					final String refBase=String.valueOf(genomic.charAt(pos1-1));
					ref= Allele.create(refBase,true);
					a1 = refBase.equalsIgnoreCase(tokens[4])?
							ref:
							Allele.create(tokens[4],false);
							;
					a2 = refBase.equalsIgnoreCase(tokens[5])?
							ref:
							Allele.create(tokens[5],false)
							;
					final String type;
					if(a1.isReference() && a2.isReference()) {
						type = "NOVARIATION";
						}
					else if(tokens[4].length() < tokens[5].length()) {
						type = "INS";
						}
					else if(tokens[4].length() > tokens[5].length()) {
						type = "DEL";
						}
					else
						{
						type = "SNV";
						}
					vcb.attribute(svtype.getID(),type);
					}
				else if(
					(tokens[4].equals("-") &&  iupacATGC.test(tokens[5])) ||
					(tokens[5].equals("-") &&  iupacATGC.test(tokens[4]))
					) {
					pos1--;//shift left
					String refBase=  String.valueOf(genomic.charAt(pos1-1));
					a1  = Allele.create(refBase,false);
					ref = Allele.create(refBase+tokens[tokens[4].equals("-")?5:4],true);
					a2=a1;
					vcb.attribute(svtype.getID(),"DEL");
					
					}
				else if(tokens[4].equals("-") && tokens[5].equals("-")) {
					pos1--;//shift left
					String refBase=  String.valueOf(genomic.charAt(pos1-1));
					a1= Allele.create(refBase,false);
					ref = Allele.create(refBase+genomic.charAt(pos1),true);
					a2 = a1;
					vcb.attribute(svtype.getID(),"DEL");
					}
				else
					{
					LOG.warn("Skipping "+line+".");
					number_non_snvs++;
					continue;
					}
				final Set<Allele> alleles = new HashSet<>();
				alleles.add(ref);
				alleles.add(a1);
				alleles.add(a2);
				
				vcb.start(pos1);
				vcb.stop(pos1+ref.length()-1);
				if(!tokens[1].isEmpty()) vcb.id(tokens[1]);
				vcb.alleles(alleles);
				
				w.add(vcb.make());
				}
			r.close();r=null;
			w.close();w=null;
			if(number_non_snvs>0) LOG.warn("Number lines skipped:"+number_non_snvs);
			return 0;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(faidx);
			CloserUtil.close(w);
			CloserUtil.close(r);
		}
	}
	public static void main(final String[] args) {
		new BimToVcf().instanceMainWithExit(args);
	}
	}
