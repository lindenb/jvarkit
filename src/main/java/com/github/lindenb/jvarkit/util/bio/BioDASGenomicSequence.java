package com.github.lindenb.jvarkit.util.bio;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.StringReader;
import java.net.URL;
import java.net.URLEncoder;
import java.util.logging.Logger;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;

public class BioDASGenomicSequence extends AbstractCharSequence implements
		ChromosomeSequence
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private boolean _length_searched=false;
	private SAXParser parser;
	private String chrom=null;
	private Integer _length=null;
	private byte buffer[]=null;
	private int buffer_pos=-1;
	private int half_buffer_capacity=1000000;
	
	public BioDASGenomicSequence(String chrom)
		{
		this.chrom=chrom;
		SAXParserFactory f=SAXParserFactory.newInstance();
		f.setSchema(null);
		f.setNamespaceAware(false);
		f.setValidating(false);
		try
			{
			this.parser=f.newSAXParser();
			}
		catch(Exception err)
			{
			throw new RuntimeException(err.getMessage());
			}
	
		}
	
	@Override
	public char charAt(int index0)
		{
		if(index0 >= length())
			{
			throw new IndexOutOfBoundsException("index:"+index0);
			}
		if(buffer!=null && index0>=buffer_pos && index0-buffer_pos < buffer.length)
			{
			return (char)buffer[index0-buffer_pos];
			}
		int minStart=Math.max(0, index0-half_buffer_capacity);
		int maxEnd=Math.min(minStart+2*half_buffer_capacity,this.length());
		try
			{
			
			String uri= "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="+
					URLEncoder.encode(getChrom()+":"+(minStart+1)+","+(maxEnd),
					"UTF-8");
			DasHandler handler=new DasHandler();
			LOG.info("getting "+uri);
			InputStream in=new URL(uri).openStream();
			this.parser.parse(in, handler);
			in.close();
			this.buffer=handler.baos.toByteArray();
			this.buffer_pos=minStart;
			}
		catch (Exception err)
			{
			throw new RuntimeException("BioDas::get",err);
			}
		return (char)buffer[index0-minStart];
		}

	@Override
	public int length()
		{
		if(_length==null && _length_searched==false)
			{
			_length_searched=true;
			LOG.info("getting length");
			try
				{
				DasHandler handler=new DasHandler();
				String uri="http://genome.ucsc.edu/cgi-bin/das/hg19/entry_points";
				LOG.info("getting length "+uri);
				InputStream in=new URL(uri).openStream();
				this.parser.parse(in, handler);
				in.close();
				this._length=handler.length;
				LOG.info("length("+getChrom()+")="+this._length);
				}
			catch (Exception err)
				{
				throw new RuntimeException("BioDas::length",err);
				}
			
			}
		return _length==null?Integer.MAX_VALUE-1:_length;
		}

	@Override
	public String getChrom()
		{
		return this.chrom;
		}
	

	
	private class DasHandler
		extends DefaultHandler
		{
		private ByteArrayOutputStream baos=null;
		int reserve=100000;
		Integer length=null;
		public DasHandler()
			{
			}
		
		 @Override
		public void startDocument() throws SAXException
			{
			baos=null;
			}
		
		 public InputSource resolveEntity (String publicId, String systemId)
	         {
	         return new InputSource(new StringReader(""));
	         }
		
		@Override
	    public void startElement(String uri, String localName, String name,
	            Attributes attributes) throws SAXException
	        {
	        if(name.equals("DNA"))
	            {
	            this.baos=new ByteArrayOutputStream(this.reserve);
	            }
	        else if(name.equals("SEGMENT") )
	        	{
	        	String c=attributes.getValue("id");
	        	if(c==null) return;
	        	if(!(getChrom().equals(c) || getChrom().equals("chr"+c))) return;
	        	c=attributes.getValue("stop");
	        	if(c==null) return;
	        	this.length=Integer.parseInt(c);
	        	}
	        }
		
	    @Override
	    public void characters(char[] ch, int start, int length)
	            throws SAXException
	        {
	        if(this.baos==null) return;
	        for(int i=0;i< length;++i)
	            {
	            char c= Character.toUpperCase(ch[start+i]);
	            if(Character.isWhitespace(c)) continue;
	            this.baos.write((byte)c);
	            }
	        }
	
	
		
		
		}
	public static void main(String[] args)
		{
		BioDASGenomicSequence seq=new BioDASGenomicSequence("chrM");
		for(int i=0;i< seq.length() ;++i) System.out.println(""+i+" : "+seq.charAt(i)+" "+seq.length());
		}

	}
