package com.github.lindenb.jvarkit.util.biomart;


import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.StringWriter;
import java.io.Writer;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineReader;



public class BiomartQuery
	{
	private String url="http://www.biomart.org/biomart/martservice/result";
	private boolean uniqRows=true;
	private String charset="UTF-8";
	private String dataSetConfigVersion="0.6";
	private String dataSetName="hsapiens_gene_ensembl";
	private List<String> attributes=Collections.emptyList();
	private Map<String,String> filters=new LinkedHashMap<String,String>();
	private boolean showHeader=false;
	public BiomartQuery()
		{
		
		}
	
	public void setShowHeader(boolean showHeader)
		{
		this.showHeader = showHeader;
		}
	
	public boolean isShowingHeader()
		{
		return showHeader;
		}
	
	public void setServiceUrl(String url)
		{
		this.url = url;
		}
	
	public String getServiceUrl()
		{
		return url;
		}
	
	public String getDataSetConfigVersion()
		{
		return dataSetConfigVersion;
		}
	public String getDataSetName()
		{
		return dataSetName;
		}
	
	public void setDataSetConfigVersion(String dataSetConfigVersion)
		{
		this.dataSetConfigVersion = dataSetConfigVersion;
		}
	
	public void setDataSetName(String dataSetName)
		{
		this.dataSetName = dataSetName;
		}
	
	public boolean isUniqRows()
		{
		return uniqRows;
		}
	
	public void setUniqRows(boolean uniqRows)
		{
		this.uniqRows = uniqRows;
		}
	
	public String getCharSet()
        {
        return charset;
        }
	
	public void setCharset(String charset)
		{
		this.charset = charset;
		}
	
	public List<String> getAttributes()
		{
		return attributes;
		}
	
	public void setAttributes(List<String> attributes)
		{
		this.attributes = new ArrayList<String>(attributes);
		}
	
	public void setAttributes(String...atts)
		{
		this.attributes = Arrays.asList(atts);
		}
	
	public void addFilter(String name,String value)
		{
		this.filters.put(name, value);
		}
	
	public Map<String,String> getFilters()
		{
		return this.filters;
		}
	
	protected void writeFilters(XMLStreamWriter w) throws XMLStreamException
		{
		Map<String,String> m=this.getFilters();
		for(String f:m.keySet())
	         {
	         w.writeEmptyElement("Filter");
	         w.writeAttribute("name", f);
	         w.writeAttribute("value", f);
	         }
		}
	
	protected void writeAttributes(XMLStreamWriter w) throws XMLStreamException
		{
		for(String att:getAttributes())
             {
             w.writeEmptyElement("Attribute");
             w.writeAttribute("name", att);
             }
		}
	
	 public void write(XMLStreamWriter w) throws XMLStreamException
	 	{
	 	w.writeStartDocument(getCharSet(), "1.0");
	     w.writeStartElement("Query");
	     w.writeAttribute("virtualSchemaName", "default");
	     w.writeAttribute("formatter", "TSV");
	     w.writeAttribute("header", isShowingHeader()?"1":"0");
	     w.writeAttribute("uniqueRows",isUniqRows()?"0":"1");
	     w.writeAttribute("count", "0");
	     w.writeAttribute("datasetConfigVersion",getDataSetConfigVersion());
	
	     w.writeStartElement("Dataset");
	     w.writeAttribute("name",getDataSetName());
	     w.writeAttribute("interface", "default");
	     writeFilters(w);
	     writeAttributes(w);
	    
	
	     w.writeEndElement();
	
	     w.writeEndElement();
	     w.writeEndDocument();
	     w.flush();
	 	}
	
	 public final void write(Writer stream) throws XMLStreamException
	     {
	     XMLOutputFactory xof=XMLOutputFactory.newFactory();
	     XMLStreamWriter w=xof.createXMLStreamWriter(stream);
	     write(w);
	     w.flush();
	     w.close();
	     }
	
	public String getXmlString()
		{
		try
			{
			StringWriter sw=new StringWriter();
			write(sw);
			return sw.toString();
			}
		catch (XMLStreamException e)
			{
			throw new RuntimeException(e);
			}
		}
	 
	@Override
	public String toString()
		{
		return getXmlString();
		} 
	
	
	  public String getMd5() throws Exception
	      {
	      MessageDigest complete = MessageDigest.getInstance("MD5");
	      complete.update(getXmlString().getBytes());
	      complete.update(getServiceUrl().getBytes());
	      byte b[]=complete.digest();
	      StringBuilder result = new StringBuilder();
	       for (int i=0; i < b.length; i++) {
	       result.append( Integer.toString( ( b[i] & 0xff ) + 0x100, 16).substring( 1 ));
	       }
	      return result.toString();
	      }
	
	public LineReader  execute() throws IOException
	     {
	     URLConnection connection = new URL(this.getServiceUrl()).openConnection();
	     connection.setDoOutput(true); 
	     connection.setRequestProperty("Accept-Charset", this.getCharSet());
	     connection.setRequestProperty("Content-Type", "application/x-www-form-urlencoded;charset=" + this.getCharSet());
	     if(connection instanceof HttpURLConnection )
	             {
	             HttpURLConnection httpConnection = (HttpURLConnection)connection;
	             httpConnection.setRequestMethod("POST");
	             httpConnection.setInstanceFollowRedirects(true);
	             }
	
	     String q="query="+URLEncoder.encode(this.getXmlString(),this.getCharSet());
	     OutputStream output = null;
	     try {
	          output = connection.getOutputStream();
	          output.write(q.getBytes(this.getCharSet()));
	          output.flush();
	     } finally {
	          if (output != null) try { output.close(); } catch (IOException logOrIgnore) {logOrIgnore.printStackTrace();}
	     }
	     InputStream response = connection.getInputStream();
	     AsciiLineReader r=new AsciiLineReader(response);
	     return r;
	     }


	}
