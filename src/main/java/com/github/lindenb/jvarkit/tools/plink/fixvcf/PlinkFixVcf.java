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
package com.github.lindenb.jvarkit.tools.plink.fixvcf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.iterator.LineIterators;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/**
BEGIN_DOC

## Example

```
plink --recode vcf bgz -o jeter (...)
gunzip -c  jeter.vcf.gz | java -jar jvarkit.jar plinkfixvcf -R reference.fa
```

END_DOC
*/
@Program(name="plinkfixvcf",
description="Fix plink --recode vcf",
keywords={"plink","vcf"},
creationDate = "20251016",
modificationDate="20251016",
jvarkit_amalgamion =  true
)
public class PlinkFixVcf extends Launcher {
	private static final Logger LOG = Logger.of(PlinkFixVcf.class);
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path fastaRef=null;
	@Parameter(names={"-D","--discarded"},description="Save skipped/ignored variants in that file")
	private Path discardedFile =null;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile =null;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate=new WritingVariantsDelegate();

	private int scan(BufferedReader br) throws IOException {
		try(PrintWriter failedPw = (discardedFile==null?NullOuputStream.newPrintWriter():IOUtils.openPathForPrintWriter(discardedFile))) {
			GenomicSequence genomic=null;
			final LineIterator li = LineIterators.of(br);
			final List<String> header_lines=new ArrayList<>();
			while(li.hasNext() && li.peek().startsWith("#")) {
				String s = li.next();
				failedPw.println(s);
				if(s.equals("##fileformat=VCFv4.3")) s="##fileformat=VCFv4.2";//htsjdk doesn't like plink's  4.2...
				if(s.startsWith(VCFConstants.CONTIG_HEADER_START)) continue;
				header_lines.add(s);
				}
			final VCFHeader hdr = VCFUtils.parseHeader(header_lines).header;
						
			final VCFCodec codec = new VCFCodec();
			final VCFInfoHeaderLine info_swapped = new VCFInfoHeaderLine("SWAP", 1, VCFHeaderLineType.Flag, "REF/ALT and genotypes were swapped to fit the reference sequence");
			hdr.addMetaDataLine(info_swapped);
			final VCFInfoHeaderLine info_shifted = new VCFInfoHeaderLine("SHIFT", 1, VCFHeaderLineType.Flag, "POS was shiffted to left for indels to fit the VCF specification");
			hdr.addMetaDataLine(info_shifted);
			
			
			JVarkitVersion.getInstance().addMetaData(this, hdr);
			try(ReferenceSequenceFile ref= ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaRef)) {
				final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(ref);
				hdr.setSequenceDictionary(dict);
				final ContigNameConverter contig_converter = ContigNameConverter.fromOneDictionary(dict);
				try(VariantContextWriter w= writingVariantsDelegate.dictionary(dict).open(outputFile)) {
					codec.setVCFHeader(hdr, VCFHeaderVersion.VCF4_2);//htsjdk cannot write above 4.3
					w.writeHeader(hdr);
					while(li.hasNext()) {
						final String line = li.next();
						final String[] tokens = CharSplitter.TAB.split(line);
						if(hdr.hasGenotypingData() && tokens.length!=hdr.getNGenotypeSamples()+9) {
							throw new JvarkitException.TokenErrors(hdr.getNGenotypeSamples()+9, tokens);
							}
						else if(!hdr.hasGenotypingData() && tokens.length!=8) {
							throw new JvarkitException.TokenErrors(8, tokens);
							}
						final String ctg = contig_converter.apply(tokens[0]);
						if(StringUtils.isBlank(ctg)) {
							LOG.warn("skipping chromosome "+tokens[0]);
							failedPw.println(line);
							continue;
							}
						tokens[0]=ctg;
						boolean shifted_flag = false;
						int pos1 = Integer.parseInt(tokens[1]);
						if(pos1<1 || pos1 > dict.getSequence(ctg).getLengthOnReference()) {
							LOG.warn("wrong position "+pos1+" in  "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
							failedPw.println(line);
							continue;
							}
						if(tokens[3].equals("0") && tokens[3].equals(".")) {
							LOG.warn("skipping  "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
							failedPw.println(line);
							continue;
							}
						if(tokens[4].contains(",")) {
							LOG.warn("multiple ALT not handled in  "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
							failedPw.println(line);
							continue;
							}
						if(tokens[3].startsWith("<") || tokens[4].startsWith("<")) {
							LOG.warn("symbolic alleles not handled in  "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
							failedPw.println(line);
							continue;
							}
						
						if(genomic==null || !genomic.hasName(ctg)) {
							genomic = new GenomicSequence(ref, ctg);
							}
						if(tokens[3].equals("-") && AcidNucleics.isIUPAC(tokens[4])) {
							pos1--;
							shifted_flag =true;
							final  char base_before = genomic.charAt(pos1-1);
							tokens[1] = String.valueOf(pos1);
							tokens[3]= String.valueOf(base_before);
							tokens[4]= String.valueOf(base_before)+tokens[4];
							}
						else if(tokens[4].equals("-") && AcidNucleics.isIUPAC(tokens[3])) {
							pos1--;
							shifted_flag =true;
							final char base_before = genomic.charAt(pos1-1);
							tokens[1] = String.valueOf(pos1);
							tokens[3]= String.valueOf(base_before)+tokens[3];
							tokens[4]= String.valueOf(base_before);
							}
						if(tokens[3].equals("0")) {
							final char base = genomic.charAt(pos1-1);
							tokens[3]= String.valueOf(base);
							}
						if(tokens[4].equals("0")) {
							final char base = genomic.charAt(pos1-1);
							tokens[4]= String.valueOf(base);
							}
						
						if(!AcidNucleics.isATGCN(tokens[3]) && !tokens[3].equals(".")) {
							LOG.warn("REF is not ATGCN/. in "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
							failedPw.println(line);
							continue;
							}
						if(!AcidNucleics.isATGCN(tokens[4]) && !tokens[4].equals(".")) {
							LOG.warn("ALT is not ATGCN/. in "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
							failedPw.println(line);
							continue;
							}
						/*
						*/
						
						final boolean is_reference_REF = genomic
									.subSequence(pos1-1, pos1-1 + tokens[3].length())
									.toString()
									.equalsIgnoreCase(tokens[3]);
						final boolean is_reference_ALT = genomic
								.subSequence(pos1-1, pos1-1 + tokens[4].length())
								.toString()
								.equalsIgnoreCase(tokens[4]);
						
						if(is_reference_REF && is_reference_ALT && tokens[3].length()==tokens[4].length()) {
							LOG.warn("both alleles look like REF in "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
							failedPw.println(line);
							continue;
							}
						if(!is_reference_REF && !is_reference_ALT && !(tokens[3].equals("N") && tokens[4].equals("."))) {
							LOG.warn("both alleles look like ALT in "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
							failedPw.println(line);
							continue;
							}
						
						/* replace non ATGC in alleles for example if there
						 * is a N while it's another base in the REF */
						for(int x=3;x<5;x++) {
							final String allele = tokens[x];
							if(AcidNucleics.isATGC(allele)) continue;
							if(!AcidNucleics.isIUPAC(allele)) continue;
							final StringBuilder sb=new StringBuilder(allele);
							for(int i=0;i< allele.length();i++) {
								if(AcidNucleics.isATGC(allele.charAt(i))) continue;
								final char c = genomic.charAt(pos1+i-1);
								//if(!AcidNucleics.isATGC(allele)) c='N'; non pas ici
								sb.setCharAt(i,c);
								}
							tokens[x]=sb.toString();
							}
						
						boolean swapped_allele=false;
						// swap alleles
						if(!is_reference_REF && is_reference_ALT) {
							final String format=tokens[8];
							final int gt_idx= CharSplitter.COLON.splitAsStringList(format).indexOf(VCFConstants.GENOTYPE_KEY);
							if(gt_idx==-1) {
								LOG.warn("No GT field in  "+ String.join(" ",Arrays.asList(tokens).subList(0, 9)));
								failedPw.println(line);
								continue;
								}
							swapped_allele= true;
							boolean fail_swap=false;
							final String swap = tokens[3];
							tokens[3]=tokens[4];
							tokens[4]=swap;
							for(int i=9;i< tokens.length;i++) {
								final List<String> tokens2 = new ArrayList<>(CharSplitter.COLON.splitAsStringList(tokens[i]));
								if(gt_idx>=tokens2.size()) {
									LOG.warn("gt missing in "+tokens[0]+":"+tokens[1]);
									continue;
									}
								final String s=tokens2.get(gt_idx);
								if(!s.matches("[01\\.][/\\|][01\\.]") && !s.matches("[01\\.]")) {
									LOG.warn("cannot swap genotype "+tokens[i] +" in "+ String.join(" ",Arrays.asList(tokens).subList(0, 5)));
									fail_swap=true;
									}
								if(s.contains("|")) {
									LOG.warn("phased genotype "+tokens[i] +" poorly supported");
									}
								tokens2.set(gt_idx, s
									.replace('1', 'x')
									.replace('0', '1')
									.replace('x', '0')
									);
								tokens[i]=String.join(":", tokens2);
								}
							if(fail_swap==true) {
								failedPw.println(line);
								continue;
								}
							}
						
						/* replace non ATGC in alleles for compatibility with htsjdk */
						for(int x=3;x<5;x++) {
							final String allele = tokens[x];
							if(AcidNucleics.isATGC(allele)) continue;
							if(!AcidNucleics.isIUPAC(allele)) continue;
							final StringBuilder sb=new StringBuilder(allele);
							for(int i=0;i< allele.length();i++) {
								if(AcidNucleics.isATGC(allele.charAt(i))) continue;
								sb.setCharAt(i,'N');
								}
							tokens[x]=sb.toString();
							}
						
						// no variant
						if(tokens[3].equals(tokens[4])) {
							tokens[4]= ".";
							}
						final VariantContext ctx= codec.decode(String.join("\t", tokens));
						final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						if(ctx.hasAttribute("PR")) {
							vcb.rmAttribute("PR");
							}
						if(swapped_allele) {
							vcb.attribute(info_swapped.getID(), Boolean.TRUE);
							}
						if(shifted_flag) {
							vcb.attribute(info_shifted.getID(), Boolean.TRUE);
							}
						w.add(vcb.make());
						}
					}
				}
			failedPw.flush();
			}
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = super.oneFileOrNull(args);
			if(input==null || input.equals("-")) {
				try(BufferedReader br= IOUtils.openStreamForBufferedReader(stdin())) {
					return scan(br);
					}
				}
			else
				{
				try(BufferedReader br= IOUtils.openURIForBufferedReading(input)) {
					return scan(br);
					}
				}
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	public static void main(final String[] args) {
		new PlinkFixVcf().instanceMainWithExit(args);
	}

}
