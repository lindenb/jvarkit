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
package com.github.lindenb.jvarkit.tools.vcfrebase;

import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.Rebase;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfRebase extends AbstractVcfRebase {
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfRebase.class);

	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private Rebase rebase=Rebase.createDefaultRebase();

	private VcfRebase() {
		}
	
	@Override
	protected Collection<Throwable> doVcfToVcf(String inputName,
			VcfIterator in, VariantContextWriter out) throws IOException {
		final String ATT="ENZ";
		GenomicSequence genomicSequence=null;
		VCFHeader header=in.getHeader();
		addMetaData(header);
		header.addMetaDataLine(new VCFInfoHeaderLine(ATT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Enzyme overlapping: Format: (Name,Site,Sequence,pos-1,strand)"));
		out.writeHeader(header);
		while(in.hasNext())
			{
			VariantContext var=in.next();

			if(genomicSequence==null || !genomicSequence.getChrom().equals(var.getContig()))
				{
				LOG.info("Loading sequence "+var.getContig());
				genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile,var.getContig());
				}
			
			Set<String> hits=new HashSet<String>();
			for(Rebase.Enzyme enz:this.rebase)
				{
				int start0=Math.max(0, var.getStart() - enz.size());
				for(int y=start0;y<=var.getStart();++y)
					{
					//run each strand
					for(int strand=0;strand<2;++strand)
						{
						int x=0;
						//loop over bases of the enzyme
						for(x=0;x< enz.size() && y+x < genomicSequence.length() ;++x )
							{
							char c=(strand==0?
									enz.at(x):
									AcidNucleics.complement(enz.at((enz.size()-1)-x))
									);
							if(!Rebase.compatible(genomicSequence.charAt(y+x),c)) break;
							}
						// match found
						if(x==enz.size())
							{
							StringBuilder b=new StringBuilder("(");
							b.append(enz.getName());
							b.append("|");
							b.append(enz.getDecl());
							b.append("|");
							for(x=0;x < enz.size();++x)
								{
								char c=genomicSequence.charAt(y+x);
								if(y+x>=var.getStart()-1 && y+x<=var.getEnd()-1)
									{
									c=Character.toLowerCase(c);
									}
								b.append(c);
								}
							b.append("|");
							b.append(y+1);
							b.append("|");
							b.append(strand==0?"+":"-");
							b.append(")");
							hits.add(b.toString());
							break;
							}
						if(enz.isPalindromic()) break;
						}
					}
				}
			if(hits.isEmpty())
				{
				out.add(var);
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(var);
			vcb.attribute(ATT, hits.toArray(new String[hits.size()]));
			out.add(vcb.make());
			}
		return RETURN_OK;
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		
		if(!super.selEnzymesStr.isEmpty())
			{
			Rebase rebase2=new Rebase();
			for(String e:selEnzymesStr)
				{
				if(e.isEmpty()) continue;
				Rebase.Enzyme enz=this.rebase.getEnzymeByName(e);
				if(enz==null)
					{
					StringBuilder msg= new StringBuilder();
					msg.append("Cannot find enzyme \""+e +"\" in RE list.\n");
					msg.append("Current list is:\n");
					for(Rebase.Enzyme E: this.rebase)
						{
						msg.append("\t"+E+"\n");
						}
					return wrapException(msg.toString());
					}
				rebase2.getEnzymes().add(enz);
				}
			this.rebase=rebase2;
			}
		
		int i=0;
		while(i< rebase.size())
			{
			if(rebase.get(i).getWeight()< weight)
				{
				rebase.getEnzymes().remove(i);
				}
			else
				{
				++i;
				}
			}
		
		if(rebase.size()==0)
			{
			LOG.warn("REBASE IS EMPTY");
			}

		if(fasta==null)
			{
			return wrapException("reference.undefined");
			}
		try
			{
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(fasta);
			return doVcfToVcf(inputName);
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		}

	public static void main(String[] args)
		{
		new VcfRebase().instanceMainWithExit(args);
		}
	}
