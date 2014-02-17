package com.github.lindenb.jvarkit.util.vcf;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;

import javax.xml.ws.WebServiceException;


import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;


public class TabixVcfFileReader implements Closeable
	//,Iterable<VariantContext> NO, not a true iterator
	{
    private final VCFCodec vcfCodec = new VCFCodec();
    private final VCFHeader vcfHeader;
    private TabixFileReader tabix=null;
    private String uri;
    private static class LR implements LineReader
    	{
    	LinkedList<String> stack=null;
    	LR(LinkedList<String> stack)
    		{
    		this.stack=stack;
    		}
    	@Override
    	public String readLine() throws IOException {
    		return (stack.isEmpty()?null:stack.removeFirst());
    		}
    	@Override
    	public void close() {
    		
    		}
    	}
    
    public TabixVcfFileReader(String uri) throws IOException
    	{
    	this.uri=uri;
    	this.tabix=new TabixFileReader(uri);
    	
    	LinkedList<String> stack=new LinkedList<>();
    	String line;
    	while((line=tabix.readLine())!=null && line.startsWith("#"))
    		{
    		stack.add(line);
    		if(line.startsWith("#CHROM\t")) break;
    		}
    	this.vcfHeader=(VCFHeader)vcfCodec.readActualHeader(
    			new LineIteratorImpl(new LR(stack))
    			);
    	}
    
    public String getURI()
    	{
    	return this.uri;
    	}
    
	public VCFCodec getCodec()
		{
		return this.vcfCodec;
		}
	
	public VCFHeader getHeader()
		{
	    return this.vcfHeader;
		}
    public Iterator<VariantContext> iterator(String chrom)
		{
    	return iterator(tabix.iterator(chrom));
		}

    public Iterator<VariantContext> iterator(String chrom,int start)
		{
    	return iterator(tabix.iterator(chrom+":"+start));
		}
    public Iterator<VariantContext> iterator(String chrom,int start,int end)
    	{
    	return iterator(tabix.iterator(chrom+":"+start+"-"+end));
    	}
    private  Iterator<VariantContext> iterator(Iterator<String> delegate)
		{
		return new MyIterator(delegate);
		}
    
    @Override
    public void close() throws WebServiceException {
    	this.tabix.close();
    	}
    
    private class MyIterator
    	implements Iterator<VariantContext>
    	{
    	Iterator<String> delegate;
    	MyIterator(Iterator<String> delegate)
    		{
    		this.delegate=delegate;
    		}
    	@Override
    	public boolean hasNext()
    		{
    		return delegate.hasNext();
    		}
    	@Override
    	public VariantContext next() {
    		return getCodec().decode(delegate.next());
    		}
    	@Override
    	public void remove() {
    		throw new UnsupportedOperationException();
    		}
    	}	
    
    @Override
    public String toString() {
    	return getURI();
    	}
	}
