/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.lumpysv;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import javafx.scene.control.Tab;


/**
BEGIN_DOC

## See also:

  *  https://github.com/arq5x/lumpy-sv/blob/master/scripts/l_sort.py


END_DOC
 */
@Program(name="lumpysort",
description="Java based version of Lupy l_sort.py. Sort Lumpy-SV VCFs on disk.",
keywords={"lumpy","vcf","sort"},
generate_doc=false
)
public class LumpySort 
	 extends Launcher {
	
	private static final Logger LOG = Logger.build(LumpySort.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private Launcher.WritingSortingCollection sortingArgs = new Launcher.WritingSortingCollection();
	
	private final List<VCFUtils.CodecAndHeader> codecAndHeaders = new ArrayList<>();
	
	private class LumpyVar
		implements Comparable<LumpyVar>
		{
		final int fileidx;
		final String line;
		private VariantContext ctx = null;
		LumpyVar(final int fileidx,final String line)
			{
			this.fileidx = fileidx;
			this.line = line;
			}
		VariantContext getCtx() {
			if(ctx==null)
				{
				this.ctx= LumpySort.this.codecAndHeaders.get(this.fileidx).codec.decode(this.line);
				}
			return this.ctx;
			}
		private String getSvType()
			{
			return getCtx().getAttributeAsString("SVTYPE","");
			}
		private String getContigLeft()
			{
			return getCtx().getContig();
			}
		
		private String getContigRight()
			{
			if(getSvType().equals("BND"))
				{
				final String alleleStr = getCtx().getAlternateAllele(0).getDisplayString();
				//eg: '[chr1:6660[N' or 'N]chr1:6619783]'.
				int colon = alleleStr.indexOf(':');
				if(colon==-1) throw new IllegalArgumentException("bad BND "+alleleStr);
				int x=colon-1;
				while(x>0)
					{
					if(alleleStr.charAt(x)=='[' || alleleStr.charAt(x)==']') break;
					x--;
					}	
				return alleleStr.substring(x+1, colon);
				}
			else
				{
				return getCtx().getContig();
				}
			}
		
		private String getStrands()
			{
			return getCtx().getAttributeAsString("STRANDS","").substring(0,2);
			}
		
		private int _end()
			{
			if(getSvType().equals("BND"))
				{
				final String alleleStr = getCtx().getAlternateAllele(0).getDisplayString();
				//eg: '[chr1:6660[N' or 'N]chr1:6619783]'.
				int colon = alleleStr.indexOf(':');
				if(colon==-1) throw new IllegalArgumentException("bad BND "+alleleStr);
				int x=colon+1;
				while(x<alleleStr.length())
					{
					if(alleleStr.charAt(x)=='[' || alleleStr.charAt(x)==']') break;
					x++;
					}	
				return Integer.parseInt(alleleStr.substring(colon+1,x));
				}
			else
				{
				return getCtx().getEnd();
				}
			}
		
		private int getStartLeft() {
			return  getCtx().getStart() +  getCtx().getAttributeAsIntList("CIPOS",0).get(0);
		}
		private int getEndLeft() {
			return  getCtx().getStart() +  getCtx().getAttributeAsIntList("CIPOS",0).get(1);
		}
		
		
		private int getStartRight() {
			return _end() +  getCtx().getAttributeAsIntList("CIEND",0).get(0);
		}
		private int getEndRight() {
			return _end() +  getCtx().getAttributeAsIntList("CIEND",0).get(1);
		}
		
		@Override
		public int compareTo(final LumpyVar o) {
			String s1 = this.getSvType();
			String s2 = o.getSvType();
			int i= s1.compareTo(s2);
			if(i!=0) return i;
			//
			s1 = this.getContigLeft();
			s2 = o.getContigLeft();
			i= s1.compareTo(s2);
			if(i!=0) return i;
			//
			s1 = this.getContigRight();
			s2 = o.getContigRight();
			i= s1.compareTo(s2);
			if(i!=0) return i;
			//
			s1 = this.getStrands();
			s2 = o.getStrands();
			i= s1.compareTo(s2);
			if(i!=0) return i;
			//
			int i1 = this.getStartLeft();
			int i2 = o.getStartLeft();
			i = i1 - i2;
			if(i!=0) return i;
			//
			i1 = this.getEndLeft();
			i2 = o.getEndLeft();
			i = i1 - i2;
			if(i!=0) return i;
			//
			i1 = this.getStartRight();
			i2 = o.getStartRight();
			i = i1 - i2;
			if(i!=0) return i;
			//
			i1 = this.getEndRight();
			i2 = o.getEndRight();
			i = i1 - i2;
			if(i!=0) return i;
			
			return line.compareTo(o.line);
			}
		
		}
	
	private class LumpyVarCodec
		extends AbstractDataCodec<LumpyVar>
		{
		@Override
		public LumpyVar decode(final DataInputStream dis) throws IOException {
			int i;
			try
				{
				i = dis.readInt();
				}
			catch(Exception err)
				{
				return null;
				}
			final String l= readString(dis);
			return new LumpyVar(i, l);
			}
		@Override
		public void encode(final DataOutputStream dos, final  LumpyVar v) throws IOException {
			dos.writeInt(v.fileidx);
			writeString(dos, v.line);
			}
		@Override
		public LumpyVarCodec clone() {
			return new LumpyVarCodec();
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
	VariantContextWriter vcw = null;
	LineIterator vcfIn= null;
	SortingCollection<LumpyVar> sorting = null;
	final List<File> inputs = IOUtil.unrollFiles(
			args.stream().map(S->new File(S)).collect(Collectors.toList()),
			".vcf",".vcf.gz");
	if(inputs.isEmpty()) {
		LOG.error("empty vcf list");
		return -1;
		}
	try {
		sorting = SortingCollection.newInstance(
				LumpyVar.class,
				new LumpyVarCodec(),
				(V1,V2)->V1.compareTo(V2),
				this.sortingArgs.getMaxRecordsInRam(),
				this.sortingArgs.getTmpPaths()
				);
		sorting.setDestructiveIteration(true);
		VCFHeader outHeader=null;
		for(final File vcfFile : inputs)
			{
			final Pattern tab = Pattern.compile("[\t]");
			LOG.info("Read "+vcfFile);
			vcfIn  = IOUtils.openFileForLineIterator(vcfFile);
			VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(vcfIn);
			this.codecAndHeaders.add(cah);
			if(!LumpyConstants.isLumpyHeader(cah.header))
				{
				LOG.error("doesn't look like a Lumpy-SV vcf header "+vcfFile);
				return -1;
				}
			if(outHeader==null) 
				{
				outHeader= new VCFHeader(
						cah.header.getMetaDataInInputOrder().
							stream().filter(H->!H.getKey().equals("fileDate")).
							collect(Collectors.toSet()),
						Collections.singletonList("VARIOUS")
						);
				outHeader.addMetaDataLine(new VCFInfoHeaderLine(
						"SNAME",
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Source sample name"
						));
				outHeader.addMetaDataLine(new VCFInfoHeaderLine(
						"ALG",
						1,
						VCFHeaderLineType.String,
						"Evidence PDF aggregation algorithm"
						));
				}
			
			for(final String sample: cah.header.getSampleNamesInOrder())
				{
				outHeader.addMetaDataLine(new VCFHeaderLine("SAMPLE", "<ID="+sample+">"));
				}
			while(vcfIn.hasNext()) {
				final String tokens[]=tab.split(vcfIn.next());
				String infoStr = tokens[7];
				if(infoStr.contains("SECONDARY")) continue;
				if(!cah.header.getSampleNamesInOrder().isEmpty()) {
					final StringBuilder sb=new StringBuilder(infoStr);
					sb.append(";SNAME=");
					sb.append(String.join(",", cah.header.getSampleNamesInOrder()));
					infoStr = sb.toString();
					}
				if(infoStr.startsWith("SVTYPE=BND")) {
					int colon = tokens[4].indexOf(':');
					int m = colon-1;
					while(m>=0)
						{
						char c=infoStr.charAt(m);
						if(c=='[' || c==']') break;
						--m;
						}
					int n = colon+1;
					while(n<tokens[4].length())
						{
						char c=infoStr.charAt(n);
						if(c=='[' || c==']') break;
						++n;
						}
					final String o_chr = infoStr.substring(m+1, colon);
					final int o_pos = Integer.parseInt( infoStr.substring(colon+1,n));
					if(o_chr.equals(tokens[0]) && 
							(infoStr.contains("--:")!=infoStr.contains("++:"))
						)
						{
						final int neg_s = infoStr.indexOf("--:");

                        if (neg_s > 0) {
                        	final int neg_e = infoStr.indexOf(";",neg_s);
                            final String pre = infoStr.substring(0,neg_s);
                            final String mid = infoStr.substring(neg_s,neg_e);
                            final String post=infoStr.substring(neg_e);
                            infoStr = pre + "++:0," + mid + post;
                        	}
                        else
                        	{
                        	final int pos_s = infoStr.indexOf("++:");
                            int pos_e = infoStr.indexOf(";",pos_s) ;
                            final String pre = infoStr.substring(0,pos_s);
                            final String  mid = infoStr.substring(pos_s,pos_e);
                            final String post = infoStr.substring(pos_e);
                            infoStr = pre + mid + ",--:0" + post;
                        	}
                        infoStr = "SVTYPE=INV" + infoStr.substring(10) + ";END=" + o_pos;
                        tokens[4] = "<INV>";
						}	
					
					}
				tokens[7] = infoStr;
				
				sorting.add(new LumpyVar(
						this.codecAndHeaders.size()-1,
						String.join("\t", tokens)
						));
				}
			CloserUtil.close(vcfIn);vcfIn=null;
			}
		sorting.doneAdding();
		CloseableIterator<LumpyVar> iter = sorting.iterator();
		
		vcw = super.openVariantContextWriter(this.outputFile);
		vcw.writeHeader(outHeader);
		while(iter.hasNext())
			{
			vcw.add(iter.next().getCtx());
			}
		vcw.close();vcw=null;
		iter.close();
		
		return 0;
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(vcfIn);
		CloserUtil.close(vcw);
		}
	}
	 
	public static void main(final String[] args) {
		new LumpySort().instanceMainWithExit(args);
	}
}
