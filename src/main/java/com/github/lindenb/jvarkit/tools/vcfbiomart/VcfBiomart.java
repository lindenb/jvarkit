package com.github.lindenb.jvarkit.tools.vcfbiomart;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.util.CloserUtil;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfBiomart extends AbstractVCFFilter2
	{
	private int batchSize=200;
	private String TAG="BIOMART";
	private List<Element> attributes=new ArrayList<Element>();
	private Document dom=null;
	private Element queryElement=null;
	private Element dataSetElement=null;
	private Element filterElementChromosomalLocation=null;
	private int chromColumn1=-1;
	private int startColumn1=-1;
	private int endColumn1=-1;
	private Set<Integer> visibleIndexes0=new HashSet<Integer>();
	private String serviceUrl="http://www.biomart.org/biomart/martservice/result";
	
	@Override
	public String getProgramDescription() {
		return "BiomartQueries with VCF.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfBiomart";
		}
	
	private String escapeInfo(String s)
		{
		if(s==null || s.isEmpty()) return "";
		return s.replaceAll("[ =;\t]","_").replaceAll("[_]+", " ");
		}

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		final String encoding="UTF-8";
		TransformerFactory factory=TransformerFactory.newInstance();
		Transformer transformer=null;
		try
			{
			transformer=factory.newTransformer();
			
			transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
			}
		catch(Exception err)
			{
			error(err);
			throw new IOException(err);
			}

		
		VCFHeader header=in.getHeader();
		StringBuilder desc=new StringBuilder("Biomart query.");
		if(!this.visibleIndexes0.isEmpty())
			{
			boolean first=true;
			desc.append(" Format:");
			for(Integer col:this.visibleIndexes0)
				{
				if(!first) desc.append("|"); first=false;
				desc.append(this.attributes.get(col).getAttribute("name"));
				}
			}
		
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		if(this.visibleIndexes0.isEmpty())
			{
			header.addMetaDataLine(new VCFInfoHeaderLine(TAG,
					0,
					VCFHeaderLineType.Flag,
					desc.toString()
					));
			}
		else
			{
			header.addMetaDataLine(new VCFInfoHeaderLine(TAG,
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					desc.toString()
					));
			}
		
		out.writeHeader(header);
		List<VariantContext> buffer=new ArrayList<VariantContext>(this.batchSize);
		
		for(;;)
			{
			if(!in.hasNext() || buffer.size()>=this.batchSize)
				{
				if(!buffer.isEmpty())
					{
					Set<String> locations=new HashSet<String>();
					for(VariantContext ctx:buffer)
						{
						locations.add(
							ctx.getChr()+":"+ctx.getStart()+":"+ctx.getEnd()+":1"	
							);
						}
					StringBuilder sb=new StringBuilder();
					for(String loc:locations)
						{
						if(sb.length()!=0) sb.append(",");
						sb.append(loc);
						}
					locations=null;
					
					IntervalTreeMap<String> treemap=new IntervalTreeMap<String>();
					this.filterElementChromosomalLocation.setAttribute("value",sb.toString());
					
					StringWriter xmlToSend=new StringWriter();
					try
						{
						transformer.transform(new DOMSource(this.dom), new StreamResult(xmlToSend));
						}
					catch (Exception e)
						{
						error(e);
						throw new IOException(e);
						}
					
					 info("POSTing to "+this.serviceUrl+" buffer.size="+buffer.size());
				     URLConnection connection = new URL(this.serviceUrl).openConnection();
				     connection.setDoOutput(true); 
				     connection.setRequestProperty("Accept-Charset",encoding);
				     connection.setRequestProperty("Content-Type", "application/x-www-form-urlencoded;charset=" +encoding);
				     if(connection instanceof HttpURLConnection )
				             {
				             HttpURLConnection httpConnection = (HttpURLConnection)connection;
				             httpConnection.setRequestMethod("POST");
				             httpConnection.setInstanceFollowRedirects(true);
				             }
				     
				     String q="query="+URLEncoder.encode(xmlToSend.toString(),encoding);
				     debug(q);
				     OutputStream output = null;
				     try {
				          output = connection.getOutputStream();
				          output.write(q.getBytes(encoding));
				          output.flush();
				     	}
				     finally
				     	{
				        CloserUtil.close(output);
				     	}
				    InputStream response = connection.getInputStream();
				    LineReader r=new AsciiLineReader(response);
					LineIterator li=new  LineIteratorImpl(r);
					while(li.hasNext())
						{
						String line=li.next();
						debug(line);
						String tokens[]=tab.split(line);
						debug(line+" L="+tokens.length );
						Interval interval=new Interval(
								tokens[this.chromColumn1-1],
								Integer.parseInt(tokens[this.startColumn1-1]),
								Integer.parseInt(tokens[this.endColumn1-1])
								);
						boolean foundSomething=false;
						StringBuilder content=new StringBuilder();
						for(Integer col: this.visibleIndexes0)
							{
							if(content.length()!=0) content.append("|");
							String s2=col>=tokens.length?"":tokens[col];
							if(!s2.trim().isEmpty()) foundSomething=true;
							content.append(escapeInfo(s2)); 
							}
						if(foundSomething || this.visibleIndexes0.isEmpty())
							{
							treemap.put(interval,content.toString());
							}
						}
					CloserUtil.close(r);
					CloserUtil.close(response);
					
					for(VariantContext ctx:buffer)
						{
						List<String> array=new ArrayList<String>(new HashSet<String>(treemap.getOverlapping(new Interval(ctx.getChr(),ctx.getStart(),ctx.getEnd()))));
						if(!array.isEmpty())
							{
							VariantContextBuilder vcb=new VariantContextBuilder(ctx);
							if(this.visibleIndexes0.isEmpty())//just a TAG
								{
								vcb.attribute(TAG, true);
								}
							else
								{
								vcb.attribute(TAG, array);
								}
							ctx=vcb.make();
							}
						out.add(ctx);
						}
					
					}
				if(!in.hasNext()) break;
				buffer.clear();
				}
			buffer.add(in.next());
			}

		}
	
	private int findColumn1(String tag)
		{
		int column=-1;
		for(int i=0;i< attributes.size();++i)
			{
			Element e=attributes.get(i);
			Attr att=e.getAttributeNode(tag);
			if(att!=null && att.getValue().equals("true"))
				{
				info("Attribute @"+tag+" was specified in the XML");
				if(column!=-1)
					{
					error("XML: Two @"+tag+"=true ?");
					return -1;
					}
				column=(i+1);
				}
			}
		
		for(int i=0; column==-1 && i< attributes.size();++i)
			{
			Element e=attributes.get(i);
			Attr att=e.getAttributeNode("name");
			if(att==null) throw new IllegalStateException("@name ?");
			if(tag.equals("chrom"))
				{
				if(att.getValue().equals("chromosome_name") ||
					att.getValue().equals("chrom_name"))
					{
					column=(i+1);
					}
				}
			else if(tag.equals("start"))
				{
				if(att.getValue().equals("start_position"))
					{
					column=(i+1);
					}
				}
			else if(tag.equals("end"))
				{
				if(att.getValue().equals("end_position"))
					{
					column=(i+1);
					}
				}
			}
			
		if(column<1 || column> this.attributes.size())
			{
			error("Cannot use column for \""+tag+"\".");
			return -1;
			}
		return column;
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -X (XML-file) XML biomart template.");
		out.println(" -n (int) batch size. default:"+this.batchSize);
		out.println(" -T (string) VCF output tag.");
		out.println(" -C (int) column index (1-based) for chromosome . Optional.");
		out.println(" -S (int) column index (1-based) for start . Optional.");
		out.println(" -E (int) column index (1-based) for end . Optional.");
		out.println(" -u (url) biomart service url. default:"+this.serviceUrl);
		
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		String xmlTemplate=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "n:T:X:C:S:E:u:"))!=-1)
			{
			switch(c)
				{
				case 'u': this.serviceUrl=opt.getOptArg();break;
				case 'C': this.chromColumn1=Integer.parseInt(opt.getOptArg());break;
				case 'S': this.startColumn1=Integer.parseInt(opt.getOptArg());break;
				case 'E': this.endColumn1=Integer.parseInt(opt.getOptArg());break;
				case 'X': xmlTemplate=opt.getOptArg();break;
				case 'n': this.batchSize=Math.max(1, Integer.parseInt(opt.getOptArg())); break;
				case 'T': this.TAG=opt.getOptArg(); break;
				default: 
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(xmlTemplate==null)
			{
			error("Undefined XML template");
			return -1;
			}
		try
			{	
			DocumentBuilderFactory f=DocumentBuilderFactory.newInstance();
			f.setCoalescing(true);
			f.setNamespaceAware(false);
			f.setValidating(false);
			f.setExpandEntityReferences(true);
			f.setIgnoringComments(false);
			f.setIgnoringElementContentWhitespace(true);
			DocumentBuilder docBuilder= f.newDocumentBuilder();
			info("Parsing xml "+xmlTemplate);
			InputStream in=IOUtils.openURIForReading(xmlTemplate);
			this.dom=docBuilder.parse(in);
			in.close();
			this.queryElement=this.dom.getDocumentElement();
			if(this.queryElement==null || !this.queryElement.getNodeName().equals("Query"))
				{
				error("XML root is not <Query/> but "+this.queryElement.getNodeName());
				return -1;
				}
			this.queryElement.setAttribute("formatter", "TSV");
			this.queryElement.setAttribute("header", "0");
			this.queryElement.setAttribute("uniqueRows", "1");
			this.queryElement.setAttribute("count", "");
			for(Node c1=this.queryElement.getFirstChild();c1!=null;c1=c1.getNextSibling())
				{
				if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
				if(c1.getNodeName().equals("Dataset"))
					{
					if(this.dataSetElement!=null)
						{
						error("XML: Two <DataSet> elements ?");
						return -1;
						}
					this.dataSetElement=Element.class.cast(c1);
					for(Node c2=this.dataSetElement.getFirstChild();c2!=null;c2=c2.getNextSibling())
						{
						if(c2.getNodeType()!=Node.ELEMENT_NODE) continue;
						Element e2=Element.class.cast(c2);
						if(e2.getNodeName().equals("Attribute"))
							{
							if(e2.getAttributeNode("name")==null)
								{
								error("XML: Attribute without @name ?");
								return -1;
								}
							this.attributes.add(Element.class.cast(c2));
							}
						else if(c2.getNodeName().equals("Filter"))
							{
							Attr att=null;
							if((att=e2.getAttributeNode("name"))!=null && att.getValue().equals("chromosomal_region"))
								{
								if(this.filterElementChromosomalLocation!=null)
									{
									error("XML: two Filter@chromosomal_region ?");
									return -1;
									}
								this.filterElementChromosomalLocation=e2;
								}
							}
						}
					if(this.filterElementChromosomalLocation==null)
						{
						this.filterElementChromosomalLocation=this.dom.createElement("Filter");
						this.filterElementChromosomalLocation.setAttribute("name","chromosomal_region");
						this.dataSetElement.insertBefore(
								this.filterElementChromosomalLocation,
								this.dataSetElement.getFirstChild()
								);
						}
					}
				
				}
			if(this.dataSetElement==null)
				{
				error("DataSet element missing");
				return -1;
				}
			if(this.attributes.size()<2)
				{
				error("Expected at least two elementss <Attribute>");
				return -1;
				}
			
			/* find attribute column for CHROMOSOME */
			if(chromColumn1==-1)
				{
				this.chromColumn1=findColumn1("chrom");
				if(this.chromColumn1<1) return -1;
				}
			if(this.chromColumn1< 1 || this.chromColumn1> attributes.size())
				{
				error("CHROM column index out of range");
				return -1;
				}
			if(startColumn1==-1)
				{
				this.startColumn1=findColumn1("start");
				if(this.startColumn1<1) return -1;
				}
			if(this.startColumn1< 1 || this.startColumn1> attributes.size())
				{
				error("START column index out of range");
				return -1;
				}
			if(endColumn1==-1)
				{
				this.endColumn1=findColumn1("end");
				if(this.startColumn1<1) endColumn1=startColumn1;
				}
			if(this.endColumn1< 1 || this.endColumn1> attributes.size())
				{
				error("END column index out of range");
				return -1;
				}
			
			
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		
		for(int i=0;i< this.attributes.size();++i)
			{
			Element E=this.attributes.get(i);
			Attr att=E.getAttributeNode("visible");
			if(att!=null && att.getValue().equals("false")) continue;
			this.visibleIndexes0.add(i);
			}
		
		return doWork(opt.getOptInd(), args);
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfBiomart().instanceMainWithExit(args);
		}

	}
