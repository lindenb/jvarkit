package com.github.lindenb.jvarkit.util.vcf;

import java.io.Closeable;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;

import javax.xml.ws.WebServiceException;

import net.sf.picard.PicardException;

import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.TabixReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;


public class TabixFileReader implements Closeable
	//,Iterable<VariantContext> NO, not a true iterator
	{
    private final VCFCodec vcfCodec = new VCFCodec();
    private final VCFHeader vcfHeader;
    private TabixReader tabix=null;
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
    
    public TabixFileReader(String uri) throws IOException
    	{
    	this.uri=uri;
    	this.tabix=new TabixReader(uri);
    	
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
    	return iterator(tabix.parseReg(chrom));
		}

    public Iterator<VariantContext> iterator(String chrom,int start)
		{
    	return iterator(tabix.parseReg(chrom+":"+start));
		}
    public Iterator<VariantContext> iterator(String chrom,int start,int end)
    	{
    	return iterator(tabix.parseReg(chrom+":"+start+"-"+end));
    	}
    private  Iterator<VariantContext> iterator(int parseReg[])
		{
		if(parseReg==null || parseReg.length!=3 ||
				parseReg[0]==-1 || parseReg[1]>parseReg[2])
			{
			return Collections.emptyIterator();
			}
		TabixReader.Iterator titer=this.tabix.query(parseReg[0], parseReg[1],parseReg[2]);
		if(titer==null)
			{
			return Collections.emptyIterator();
			}
		
		return new MyIterator(titer);
		}
    
    @Override
    public void close() throws WebServiceException {
    	this.tabix.close();
    	}
    
    private class MyIterator
    	implements Iterator<VariantContext>
    	{
    	TabixReader.Iterator delegate;
    	String _next=null;
    	MyIterator(TabixReader.Iterator delegate)
    		{
    		this.delegate=delegate;
    		}
    	@Override
    	public boolean hasNext()
    		{
    		if(delegate==null) return false;
    		if(_next==null)
    			{
    			try
    				{
    				_next=delegate.next();
    				}
    			catch(IOException err)
    				{
    				throw new PicardException("Tabix",err);
    				}
    			if(_next==null) delegate=null;
    			}
    		return _next!=null;
    		}
    	@Override
    	public VariantContext next() {
    		if(!hasNext()) throw new IllegalStateException();
    		String s=_next;
    		_next=null;
    		return getCodec().decode(s);
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
