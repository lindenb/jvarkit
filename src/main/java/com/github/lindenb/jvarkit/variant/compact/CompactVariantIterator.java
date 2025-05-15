package com.github.lindenb.jvarkit.variant.compact;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

public class CompactVariantIterator implements VCFIterator {
	private SAMSequenceDictionary dict;
	private int n_samples=0;
	private BinaryCodec codec;
	public CompactVariantIterator(OutputStream out) {
		codec=new BinaryCodec(out);
		}

	
	private List<String> decodeVariants() {
		int n=0;
		switch(n) {
			case 0 : return Arrays.asList("A","C");
			case 1 : return Arrays.asList("A","G");
			case 2 : return Arrays.asList("A","T");
			case 3 : return Arrays.asList("C","A");
			case 4 : return Arrays.asList("C","G");
			case 5 : return Arrays.asList("C","T");
			case 6 : return Arrays.asList("G","A");
			case 7 : return Arrays.asList("G","C");
			case 8 : return Arrays.asList("G","T");
			case 9 : return Arrays.asList("T","A");
			case 10 : return Arrays.asList("T","C");
			case 11 : return Arrays.asList("T","G");

			}
		return null;
		}

	@Override
	public void setHeader(VCFHeader header) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void writeHeader(VCFHeader header) {
		
	}

	@Override
	public void close() {
		try {
			this.codec.getOutputStream().flush();
		} catch (IOException e) {
		}
		this.codec.close();
	}

	@Override
	public boolean checkError() {
		return false;
	}

	@Override
	public void add(VariantContext vc) {
		
	}

	public void addVariant(String contig,int position,String id,List<CharSequence> alleles,List<LightGenotype> genotypes) {
		SAMSequenceRecord ssr = this.dict.getSequence(contig);
		if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(contig, dict);
		final int tid = ssr.getSequenceIndex();
		if(this.dict.size() <= BinaryCodec.MAX_UBYTE) 
			{
			this.codec.writeUByte((short)tid);
			}
		else if(this.dict.size() <= BinaryCodec.MAX_USHORT) 
			{
			this.codec.writeUShort(tid);
			}
		else
			{
			this.codec.writeUInt(tid);
			}
		if(ssr.getSequenceLength() <= BinaryCodec.MAX_UBYTE) 
			{
			this.codec.writeUByte((short)position);
			}
		else if(ssr.getSequenceLength() <= BinaryCodec.MAX_USHORT) 
			{
			this.codec.writeUShort(position);
			}
		else
			{
			this.codec.writeUInt(position);
			}
		boolean write_alt=false;
		if(genotypes.size()==2) {
			String A1 = genotypes.get(0).toString();
			String A2 = genotypes.get(1).toString();
		    if(A1.equals("A") && A2.equals("C")) {this.codec.writeUByte((short)1);}
			else if(A1.equals("A") && A2.equals("G")) {this.codec.writeUByte((short)2);}
			else if(A1.equals("A") && A2.equals("T")) {this.codec.writeUByte((short)3);}
			else if(A1.equals("C") && A2.equals("A")) {this.codec.writeUByte((short)4);}
			else if(A1.equals("C") && A2.equals("G")) {this.codec.writeUByte((short)5);}
			else if(A1.equals("C") && A2.equals("T")) {this.codec.writeUByte((short)6);}
			else if(A1.equals("G") && A2.equals("A")) {this.codec.writeUByte((short)7);}
			else if(A1.equals("G") && A2.equals("C")) {this.codec.writeUByte((short)8);}
			else if(A1.equals("G") && A2.equals("T")) {this.codec.writeUByte((short)9);}
			else if(A1.equals("T") && A2.equals("A")) {this.codec.writeUByte((short)10);}
			else if(A1.equals("T") && A2.equals("C")) {this.codec.writeUByte((short)11);}
			else if(A1.equals("T") && A2.equals("G")) {this.codec.writeUByte((short)12);}
			else
				{
				write_alt=true;
				}
			}
		else
			{
			write_alt = true;
			}
		}
}
