package com.github.lindenb.jvarkit.tools.cgi;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.prefs.InvalidPreferencesFormatException;
import java.util.prefs.Preferences;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;

public abstract class AbstractCGI extends Launcher
	{
	//private static final Logger LOG=Logger.build(AbstractCGI.class).make();
	private static final String PROPERTY_PREFFILE="prefs.file.xml";
	
	protected StringBuilder logStream=new StringBuilder();
	private List<Parameter> parameters=new ArrayList<Parameter>();
	private int contentMaxLength=2048;
	private Preferences prefs=null;
	private boolean mimeHeaderPrinted=false;
	
	
	
	/** Interface for a Parameter */
	static public interface Parameter
		{
		public String getKey();
		public boolean isFile();
		}
		
	static public class DefaultParameter
		implements Parameter
		{
		private String key;
		private String value;
		public DefaultParameter(String key,String value)
			{
			this.key=key;
			this.value=value;
			}
		@Override
		public String getKey() {
			return key;
			}
		public String getValue() { return value;}
		@Override
		public final boolean isFile() {
			return false;
			} 
		}
	
	protected AbstractCGI()
		{
		super();
		final PrintStream originalStderr=System.err;
		OutputStream redirect=new OutputStream()
			{
			@Override
			public void write(int c) throws IOException
				{
				if(c==-1) return;
				originalStderr.write(c);
				if(logStream.length() < 2048)
					{
					logStream.append((char)c);
					}
				}
			};
		System.setErr(new PrintStream(redirect));
		//
		}
	
	protected List<File> getPreferenceFiles()
		{		
		String s = System.getProperty(PROPERTY_PREFFILE);
		if(s==null) return Collections.emptyList();
		return Arrays.asList(new File(s));	
		}
	
	protected InputStream openPreferences() throws IOException
		{
		for(File f:getPreferenceFiles())
			{
			if(f.exists() && f.isFile() && f.canRead())
				{
				return new FileInputStream(f);
				}
			}
		StringBuilder xml=new StringBuilder("<!DOCTYPE preferences SYSTEM 'http://java.sun.com/dtd/preferences.dtd'>");
		xml.append("<preferences EXTERNAL_XML_VERSION=\"1.0\">");
		xml.append("<root type='user'>");
		xml.append("<map/>");
		xml.append("</root>");
		xml.append("</preferences>");
		return new ByteArrayInputStream(xml.toString().getBytes());
		}
	
	
	protected Preferences getPreferences() throws IOException
		{
		if(this.prefs==null)
			{
			InputStream is=null;
			try
				{
				is=openPreferences();
				Preferences.importPreferences(is);
				}
			catch(InvalidPreferencesFormatException er)
				{
				throw new IOException(er);
				}
			finally
				{
				CloserUtil.close(is);
				}
			this.prefs=Preferences.userRoot();
			}
		return this.prefs;
		}
	
	
	public void setContentMaxLength(int contentMaxLength) {
		this.contentMaxLength = contentMaxLength;
		}
	
	public int getContentMaxLength() {
		return contentMaxLength;
		}
	
	public String getenv(String key)
		{
		return System.getenv(key);
		}
	
	public Set<String> getParameterNames()
		{
		Set<String> keys= new HashSet<String>(getParameters().size());
		for(Parameter p:getParameters())
			{
			keys.add(p.getKey());
			}
		return keys;
		}
	
	public List<Parameter> getParameters()
		{
		return this.parameters;
		}

	public List<Parameter> getParameters(String key)
		{
		List<Parameter> L= new ArrayList<Parameter>();
		for(Parameter p:getParameters())
			{
			if(p.getKey().equals(key)) L.add(p);
			}
		return L;
		}
	
	/** find simple parameter and return its value as a String */
	public String getString(String key)
		{
		Parameter p= getParameter(key);
		if(p==null) return null;
		if(p instanceof DefaultParameter)
			{
			return DefaultParameter.class.cast(p).getValue();
			}
		return null;
		}
	
	public Parameter getParameter(String key)
		{
		for(Parameter p:getParameters())
			{
			if(p.getKey().equals(key))return p;
			}
		return null;
		}
	
	
	
	public boolean hasParameter(String key)
		{
		return getParameter(key)!=null;
		}
	
	
	/** par the HTTP request */
	public void parse() throws IOException
		{
		String requestMethod=null;
	
	    if((requestMethod=getenv("REQUEST_METHOD"))==null) throw new IOException("Cannot find REQUEST_METHOD");
	    if(requestMethod.equals("POST"))
	    	{
	    	parsePOST();
	    	}
	    else if(requestMethod.equals("GET"))
	    	{
	    	parseGET();
	    	}
	    else
		    {
		    throw new IOException("Unknown REQUEST_METHOD \""+requestMethod+"\"");
			}
		}
	
	
	
	private void parseGET()  throws IOException
		{
		String queryString = getenv("QUERY_STRING");
		if(queryString==null) throw new IOException("Cannot find QUERY_STRING");
		if(queryString.length()>getContentMaxLength()) throw new IOException("Content too large");
		parse(new StringReader(queryString),queryString.length());
		}
	
	private  boolean isMultipart()
		{
		String contentType= getenv("CONTENT_TYPE");
		if(contentType==null) return false;
		return contentType.indexOf("multipart/form-data")!=-1;
		}
	
	private void parsePOST()  throws IOException
		{
		String contentLengthStr= getenv("CONTENT_LENGTH");
		if(contentLengthStr==null) throw new IOException("CONTENT_LENGTH missing");
		int contentLength=0;
		try {
			contentLength= Integer.parseInt(contentLengthStr);
		} catch (Exception e) {
			throw new IOException("Bad content Length "+contentLength);
			}
		if(contentLength>getContentMaxLength()) throw new IOException("Content too large");
	 	
		if(!isMultipart())
			{
			parse(new BufferedReader(
					new InputStreamReader(System.in)),
					contentLength);
			}
		else
			{
			throw new IOException("Cannot parse multipart actions");
			}
		}
	
	

	
	private void parse(Reader in,int maxCharRead) throws IOException
		{
		int c;
		StringBuilder key=new StringBuilder();
		StringBuilder value=null;
		int count=0;
		while((c=in.read())!=-1 && count < maxCharRead)
			{
			++count;
			if(c=='+') c=' ';
			if(c=='&')
				{
				if(key.length()!=0)
					{
					addParameter(key.toString(),value);
					}	
				key=new StringBuilder();
				value=null;
				}
			else if(c=='=' && key.length()>0 && value==null)
				{
				value = new StringBuilder();
				}
			else
				{
				if(c=='%' && count+2<=maxCharRead)
					{
					int c2= in.read();
					if(c2==-1) throw new IOException("Bad Input");
					int c3= in.read();
					if(c3==-1) throw new IOException("Bad Input");	
					c=x2c(c2,c3);
					count+=2;
					}
				if(value!=null)
					{
					value.append((char)c);
					}
				else
					{
					key.append((char)c);
					}
				}
			}
		if(key.length()!=0)
			{
			addParameter(key.toString(),value);
			}
		
		}
	
	private static int x2c(int c1,int c2) throws IOException
		{
		try {
			return Integer.parseInt(String.valueOf((char)c1)+(char)c2, 16);
			}
		catch (NumberFormatException e)
			{
			throw new IOException("cannot cast hexa int("+(c1)+")int("+c2+")",e);
			}
		}
	
	private void addParameter(String key,CharSequence value)
		{
		if(key==null || key.isEmpty()) return;
		this.parameters.add(new DefaultParameter(key.toString(),value==null?"":value.toString()));
		}
	
	
	protected boolean isMimeHeaderPrinted()
		{
		return mimeHeaderPrinted;
		}
	
	protected void setMimeHeaderPrinted(boolean mimeHeaderPrinted) {
		this.mimeHeaderPrinted = mimeHeaderPrinted;
		}	
	
	protected void writeHTMLException(XMLStreamWriter w,Throwable err) throws XMLStreamException
		{
		if(err==null) return;
		w.writeStartElement("pre");
		w.writeCharacters(String.valueOf(err.getMessage()));w.writeEmptyElement("br");
		for(StackTraceElement e:err.getStackTrace())
			{
			w.writeCharacters(" ");
			w.writeCharacters(String.valueOf(e.getClassName()));
			w.writeCharacters(" ");
			w.writeCharacters(String.valueOf(e.getMethodName()));
			w.writeCharacters(" ");
			w.writeCharacters(String.valueOf(e.getLineNumber()));
			w.writeEmptyElement("br");
			}
		w.writeEndElement();
		}
	protected void writeHtmlFooter(XMLStreamWriter w) throws XMLStreamException
		{
		w.writeEmptyElement("hr");
		w.writeStartElement("div");
		
		w.writeCharacters("Author : ");

		w.writeStartElement("a");
		w.writeAttribute("href","#");
		w.writeAttribute("title","my mail");
		w.writeCharacters("Pierre Lindenbaum");
		w.writeEndElement();
		/*
		String s=getOnlineDocUrl();
		if(s!=null && !s.isEmpty())
			{
			w.writeCharacters(" Documentation : ");
			w.writeStartElement("a");
			w.writeAttribute("href", s);
			w.writeAttribute("title", s);
			w.writeCharacters(s);
			w.writeEndElement();
			}*/
		
		w.writeCharacters(" Version: ");
		w.writeStartElement("i");
		w.writeCharacters(String.valueOf(getVersion()));
		w.writeEndElement();
		
		w.writeCharacters(" Compilation: ");
		w.writeStartElement("i");
		w.writeCharacters(String.valueOf(getCompileDate()));
		w.writeEndElement();
		
		
		w.writeEndElement();
		}
	
	abstract protected void doCGI();

	
	@Override
	final public int doWork(List<String> args)
		{
		try {
			parse();
			doCGI();
			}
		catch (Exception e)
			{
			if(!isMimeHeaderPrinted())
				{
				setMimeHeaderPrinted(true);
				System.out.print("Content-type: text/plain\n");
				System.out.println();
				}
			e.printStackTrace(System.out);
			}
		finally
			{
			System.out.flush();
			System.out.close();
			}
		return 0;
		}
	protected void writeHTMLFooter(XMLStreamException out) throws XMLStreamException
		{
		
		}
	}
