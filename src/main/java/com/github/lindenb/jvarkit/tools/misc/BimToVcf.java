/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

public class BimToVcf extends AbstractBimToVcf
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(BimToVcf.class);

	public BimToVcf() {
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		VariantContextWriter w=null;
		BufferedReader r=null;
		IndexedFastaSequenceFile faidx=null;
		GenomicSequence genomic = null;
		try {
			if(super.REF==null) {
				return wrapException("Reference -"+OPTION_REF+" missing.");
			}
			faidx = new IndexedFastaSequenceFile(super.REF);
			
			final SAMSequenceDictionary dict=faidx.getSequenceDictionary();
			if(dict==null) {
				return wrapException("No dictionary in "+super.REF);
				}
			
			r = inputName==null?
					IOUtils.openStreamForBufferedReader(stdin()):
					IOUtils.openURIForBufferedReading(inputName)
					;
				
			final Set<VCFHeaderLine> headerLines = new HashSet<>();
			final VCFInfoHeaderLine morgan= new VCFInfoHeaderLine("MORGAN", 1, VCFHeaderLineType.Float,"Centimorgan");
			final VCFInfoHeaderLine svtype= new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String,"Variation type");
			VCFStandardHeaderLines.addStandardInfoLines(headerLines, false, "");
			super.addMetaData(headerLines);
			headerLines.add(morgan);
			headerLines.add(svtype);

			final List<String> genotypeSampleNames = Collections.emptyList();
			final VCFHeader header=new VCFHeader(headerLines, genotypeSampleNames);
			header.setSequenceDictionary(dict);
			w = super.openVariantContextWriter();
			w.writeHeader(header);
			final Pattern tab=Pattern.compile("[\t]");
			String line;
			final Pattern iupacATGC = Pattern.compile("[atgcATGC]");
			while((line=r.readLine())!=null) {
				String tokens[]=tab.split(line);
				if(tokens.length!=6) return wrapException("expected 6 column in "+line);
				Allele a1=null;
				Allele a2=null;
				Allele ref=null;
				String contig=tokens[0];
			
				SAMSequenceRecord ssr= null;
				ssr = dict.getSequence(contig);
				//ugly below !!
				if(ssr==null && contig.equals("23"))
					{
					ssr = dict.getSequence("X");
					}
				if(ssr==null && contig.equals("23"))
					{
					ssr = dict.getSequence("chrX");
					}
				if(ssr==null && contig.equals("24"))
					{
					ssr = dict.getSequence("Y");
					}
				if(ssr==null && contig.equals("24"))
					{
					ssr = dict.getSequence("chrY");
					}
				if(ssr==null && contig.equals("26"))
					{
					ssr = dict.getSequence("chrM");
					}
				if(ssr==null && contig.equals("26"))
					{
					ssr = dict.getSequence("MT");
					}
				if(ssr==null && contig.equals("25")){
					LOG.warn("ignoring "+line);
					continue;
				}
				if(ssr==null){
					 return wrapException("unknown chrom in "+line);
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
				
				if(iupacATGC.matcher(tokens[4]).matches() && iupacATGC.matcher(tokens[5]).matches())
					{
					String refBase=String.valueOf(genomic.charAt(pos1-1));
					ref= Allele.create(refBase,true);
					a1 = refBase.equalsIgnoreCase(tokens[4])?
							ref:
							Allele.create(tokens[4],false)
							;
					a2 = refBase.equalsIgnoreCase(tokens[5])?
							ref:
							Allele.create(tokens[5],false)
							;
					vcb.attribute(svtype.getID(),
							a1.isReference() && a2.isReference()?
							"NOVARIATION":
							"SNV"
							);
					}
				else if(
					(tokens[4].equals("-") &&  iupacATGC.matcher(tokens[5]).matches()) ||
					(tokens[5].equals("-") &&  iupacATGC.matcher(tokens[4]).matches())
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
					return wrapException("not handled: "+line);
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
			return RETURN_OK;
		} catch (Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(faidx);
			CloserUtil.close(w);
			CloserUtil.close(r);
		}
	}
	public static void main(String[] args) {
		new BimToVcf().instanceMainWithExit(args);
	}
	}
