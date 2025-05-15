package com.github.lindenb.jvarkit.variant.compact;

import java.io.IOException;
import java.io.OutputStream;
import java.util.BitSet;
import java.util.List;

import com.github.lindenb.jvarkit.io.ByteBufferSequence;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

public class CompactVariantContextWriter extends VCFCompactBase implements VariantContextWriter {
	private SAMSequenceDictionary dict;
	private int n_samples=0;
	private BinaryCodec codec;
	private ByteBufferSequence m_buff1=new ByteBufferSequence();
	public CompactVariantContextWriter(OutputStream out) {
		codec=new BinaryCodec(out);
	}

	private void writeU16String(final String s) {
		byte[] a=s.getBytes(ENCODING);
		codec.writeUShort(a.length);
		codec.writeBytes(a);
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

	private static void writeUInt(BinaryCodec bc,int value,int max) {
		if(value>max) throw new IllegalArgumentException();
		if(value<0) throw new IllegalArgumentException();
		if(max <= BinaryCodec.MAX_UBYTE)
			{
			bc.writeUByte((short)value);
			}
		else if(max <= BinaryCodec.MAX_USHORT)
			{
			bc.writeUShort(value);
			}
		else
			{
			bc.writeUInt(value);
			}
		}
	
	private void writeUInt(int value,int max) {
		writeUInt(this.codec, value, max);
		}
	
	@Override
	public void add(VariantContext vc) {
		
	}

	public void addVariant(String contig,int position,String id,List<CharSequence> alleles,List<LightGenotype> genotypes) {
		// skip non-informative variant
		if(genotypes.stream().noneMatch(G->G.hasAlt())) return;
		
		SAMSequenceRecord ssr = this.dict.getSequence(contig);
		if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(contig, dict);
		final int tid = ssr.getSequenceIndex();
		
		writeUInt(tid,this.dict.size());
		writeUInt(position,ssr.getSequenceLength());
		
		boolean write_alt=false;
		if(alleles.size()==2) {
			final String A1 = alleles.get(0).toString().toUpperCase();
			final String A2 = alleles.get(1).toString().toUpperCase();
				 if(A1.equals("A") && A2.equals("C")) {this.codec.writeUByte((short)0);}
			else if(A1.equals("A") && A2.equals("G")) {this.codec.writeUByte((short)1);}
			else if(A1.equals("A") && A2.equals("T")) {this.codec.writeUByte((short)2);}
			else if(A1.equals("C") && A2.equals("A")) {this.codec.writeUByte((short)3);}
			else if(A1.equals("C") && A2.equals("G")) {this.codec.writeUByte((short)4);}
			else if(A1.equals("C") && A2.equals("T")) {this.codec.writeUByte((short)5);}
			else if(A1.equals("G") && A2.equals("A")) {this.codec.writeUByte((short)6);}
			else if(A1.equals("G") && A2.equals("C")) {this.codec.writeUByte((short)7);}
			else if(A1.equals("G") && A2.equals("T")) {this.codec.writeUByte((short)8);}
			else if(A1.equals("T") && A2.equals("A")) {this.codec.writeUByte((short)9);}
			else if(A1.equals("T") && A2.equals("C")) {this.codec.writeUByte((short)10);}
			else if(A1.equals("T") && A2.equals("G")) {this.codec.writeUByte((short)11);}
			else
				{
				write_alt=true;
				}
			}
		else
			{
			write_alt = true;
			}
		if(write_alt) {
			this.codec.writeUByte((short)(alleles.size()+11));
			for(int i=0;i< alleles.size();i++) {
				this.writeU16String(StringUtils.compactDNA(alleles.get(i)));
				}
			}
		
		boolean any_missing = genotypes.stream().anyMatch(G->G.isMissing());
		int min_ploidy =  genotypes.stream().mapToInt(G->G.getPloidy()).min().orElse(0);
		int max_ploidy =  genotypes.stream().mapToInt(G->G.getPloidy()).max().orElse(0);
		byte[] array1=null;
		int method=0;
		if(min_ploidy==2 && max_ploidy==2) {
			if(any_missing) 
				{
				array1 = encodeDiploidNoMissingMethod1(genotypes);
				method = 1;
				byte[] array2 = encodeDiploidNoMissingMethod2(genotypes);
				if(array2.length< array1.length) {
					array1= array2;
					method=2;
					}
				array2 = encodeDiploidNoMissingMethod3(genotypes);
				if(array2.length< array1.length) {
					array1= array2;
					method=3;
					}
				}
			else 
				{
				array1=null;
				}
			}
		this.codec.writeByte((byte)method);
		this.codec.writeUInt(array1.length);
		this.codec.writeBytes(array1);
		}
	
	
		private byte[] encodeDiploidNoMissingMethod1(List<LightGenotype> genotypes) {
			/* most common case diploid and no missing 
			 * use two bits per genotype. Put in a bitset an convert to array
			 * */
			final BitSet bits =new BitSet(2*genotypes.size());
			for(int i=0;i< genotypes.size();i++) {
				LightGenotype gt = genotypes.get(i);
				if(gt.isHomRef()) {
					//skip
					}
				else if(gt.isHet()) {
					bits.set(i*2);
					}
				else if(gt.isHomVar()) {
					bits.set(i*2+0);
					bits.set(i*2+1);
					}
				else
					{
					throw new IllegalStateException();
					}
				}
			m_buff1.clear();
			try(BinaryCodec bc2=new BinaryCodec(this.m_buff1) ){
				byte[] array= bits.toByteArray();
				bc2.writeUInt(array.length);
				bc2.writeBytes(array);
				return m_buff1.toByteArray();
				}
			}
		
		private byte[] encodeDiploidNoMissingMethod2(List<LightGenotype> genotypes) {
			int index_beg=0;
			while(index_beg < genotypes.size()) {
				if(!genotypes.get(index_beg).isHomRef()) break;
				index_beg++;
				}
			int index_end=genotypes.size();
			while(index_end >0) {
				if(!genotypes.get(index_beg-1).isHomRef()) break;
				index_end--;
				}
			
			
		
			final BitSet bits =new BitSet(2*genotypes.size());
			for(int i=index_beg;i< index_end;i++) {
				LightGenotype gt = genotypes.get(i);
				if(gt.isHomRef()) {
					//skip
					}
				else if(gt.isHet()) {
					bits.set(i*2);
					}
				else if(gt.isHomVar()) {
					bits.set(i*2+0);
					bits.set(i*2+1);
					}
				else
					{
					throw new IllegalStateException();
					}
				}
			m_buff1.clear();
			try(BinaryCodec bc2=new BinaryCodec(this.m_buff1) ){
				writeUInt(bc2,index_beg, genotypes.size());
				byte[] array= bits.toByteArray();
				bc2.writeUInt(array.length);
				bc2.writeBytes(array);
				
				return m_buff1.toByteArray();
				}
			}
		
		private byte[] encodeDiploidNoMissingMethod3(List<LightGenotype> genotypes) {
			m_buff1.clear();
			try(BinaryCodec bc2=new BinaryCodec(this.m_buff1)) {
				int prev_index=0;
				for(int i=0;i< genotypes.size();i++) {
					LightGenotype gt = genotypes.get(i);
					if(gt.isHomRef()) {
						continue;
						}
					writeUInt(bc2,i-prev_index,genotypes.size());
					prev_index=i;
					if(gt.isHet()) {
						bc2.writeUByte((byte)1);
						}
					else if(gt.isHomVar()) {
						bc2.writeUByte((byte)3);
						}
					else
						{
						throw new IllegalStateException();
						}
					}
				return m_buff1.toByteArray();
				}
			}
	}
