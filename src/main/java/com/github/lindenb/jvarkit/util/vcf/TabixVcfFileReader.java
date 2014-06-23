package com.github.lindenb.jvarkit.util.vcf;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;


public class TabixVcfFileReader extends AbstractTabixObjectReader<VariantContext>
	//,Iterable<VariantContext> NO, not a true iterator
	{
    private AbstractVCFCodec vcfCodec;
    private VCFHeader vcfHeader;
   
    
    public TabixVcfFileReader(String uri) throws IOException
    	{
    	super(uri);
    	
    	List<String> stack=new ArrayList<String>();
    	String line;
    	while((line=super.tabix.readLine())!=null && line.startsWith("#"))
    		{
    		stack.add(line);
    		if(line.startsWith("#CHROM\t")) break;
    		}
    	VCFUtils.CodecAndHeader cah=VCFUtils.parseHeader(stack);
    	this.vcfHeader=cah.header;
    	this.vcfCodec=cah.codec;
    	}
    
    
	public AbstractVCFCodec getCodec()
		{
		return this.vcfCodec;
		}
	
	public VCFHeader getHeader()
		{
	    return this.vcfHeader;
		}
    
    @Override
    protected  Iterator<VariantContext> iterator(Iterator<String> delegate)
		{
		return new MyIterator(delegate);
		}
    
    
    private class MyIterator
    	extends AbstractMyIterator
    	{
    	MyIterator(Iterator<String> delegate)
    		{
    		super(delegate);
    		}
    	@Override
    	public VariantContext next() {
    		return getCodec().decode(delegate.next());
    		}
    	}	
    
	}
