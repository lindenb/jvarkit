package com.github.lindenb.jvarkit.tools.evs2bed;

import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.soap.SOAPConstants;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;


public class DumpExomeVariantServerData
	extends AbstractCommandLineProgram
	{
	private int STEP_SIZE=25000;
    private long LIMIT=-1L;
	private long count_records=0L;
	private long genome_total_size=0L;
	private long genome_curr_size=0L;
	
	public static final String EVS_NS="http://webservice.evs.gs.washington.edu/";
	private DocumentBuilder documentBuilder;
	private Transformer transformer;
	private DumpExomeVariantServerData()
		{
		}
	
	private static Element first(Element root,String namespaceuri,String localName)
		{
		if(root==null) return null;
		for(Node n=root.getFirstChild();n!=null;n=n.getNextSibling())
			{
			if(n.getNodeType()!=Node.ELEMENT_NODE) continue;
			if(namespaceuri!=null && !namespaceuri.equals(n.getNamespaceURI())) continue;
			if(namespaceuri!=null && !localName.equals(n.getLocalName())) continue;
			if(namespaceuri==null && !localName.equals(n.getNodeName())) continue;
			return Element.class.cast(n);
			}
		return null;
		}
	private final int MAX_TRY=10; 
	private Element fetchEvsData(String chrom,int start,int end)
		{
		double ratio=100.0*(this.genome_curr_size+start)/(double)this.genome_total_size;
		
		info(chrom+":"+start+"-"+end+ " N="+count_records+" "+(int)ratio+"%");
		try
			{
		    URL url = new URL("http://gvs-1.gs.washington.edu/wsEVS/EVSDataQueryService");
	
		    // Send data
		    URLConnection conn = null;
		    for(int n_try=0;n_try<MAX_TRY;++n_try)
			    {
		    	try
		    		{
		    		conn=url.openConnection();
		    		}
		    	catch(java.net.ConnectException err)
		    		{
		    		if(n_try+1==MAX_TRY) throw err;
		    		System.err.println(
		    			"Error: trying "+(n_try)+"/"+MAX_TRY+" "+url
		    			);
		    		}	
			    }
		    conn.setDoOutput(true);
		    PrintStream wr=new PrintStream(conn.getOutputStream());
		    wr.print("<?xml version='1.0' ?>"+
		    		"<S:Envelope xmlns:S='http://schemas.xmlsoap.org/soap/envelope/'>"+
		    		  "<S:Body>"+
		    		    "<ns2:getEvsData xmlns:ns2='http://webservice.evs.gs.washington.edu/'>"+
		    		      "<arg0>"
		    		);
		    wr.print(chrom);
		    wr.print(":");
		    wr.print(String.valueOf(start));
		    wr.print("-");
		    wr.print(String.valueOf(end));
		    wr.print("</arg0>"+
	    		    "</ns2:getEvsData>"+
	    		  "</S:Body>"+
	    		"</S:Envelope>"
	    		);
		    wr.flush();
		    InputStream rd = conn.getInputStream();
		    Document dom=this.documentBuilder.parse(rd);
		    wr.close();
		    rd.close();
		    Element e=first(dom.getDocumentElement(), SOAPConstants.URI_NS_SOAP_ENVELOPE, "Body");
		    e=first(e, EVS_NS, "getEvsDataResponse");
		    e=first(e, null, "return");
		    if(e==null) return null;
			return e;
			}
		catch(Exception err)
			{
			err.printStackTrace();
			return null;
			}
		}
	
	private class Fetcher
		{
		String chrom;
		int length;
		Fetcher(String chrom,int length)
			{
			this.chrom=chrom;
			this.length=length;
			}
		
		public void run() throws Exception
			{
			if(DumpExomeVariantServerData.this.LIMIT>0 && DumpExomeVariantServerData.this.count_records>=DumpExomeVariantServerData.this.LIMIT) return;
			final int step=DumpExomeVariantServerData.this.STEP_SIZE;
			int start=1;
			do
				{
				
				Element root=fetchEvsData(chrom,start,start+step+10);
				for(Node n=(root==null?null:root.getFirstChild());n!=null;n=n.getNextSibling())
					{
					if(n.getNodeType()!=Node.ELEMENT_NODE) continue;
					if(!n.getNodeName().equals("snpList")) continue;
					String chromosome=null;
					String chrPosition=null;
					for(Node n2=n.getFirstChild();n2!=null;n2=n2.getNextSibling())
						{
						if(n2.getNodeType()!=Node.ELEMENT_NODE) continue;
						if(n2.getNodeName().equals("chromosome"))
							{
							chromosome=n2.getTextContent();
							}
						else if(n2.getNodeName().equals("chrPosition"))
							{
							chrPosition=n2.getTextContent();
							}
						}
					count_records++;
					if(LIMIT>0 && count_records>=LIMIT) break;
					
					StringWriter sw=new StringWriter();
					transformer.transform(
							new DOMSource(n),
							new StreamResult(sw)
							);
					sw.flush();
					String xml=sw.toString().replace("\n", "");
					
					System.out.print(chromosome);
					System.out.print('\t');
					System.out.print(Integer.parseInt(chrPosition)-1);
					System.out.print('\t');
					System.out.print(chrPosition);
					System.out.print('\t');
					System.out.println(xml);
					}
				
				
				start+=step;
				if(LIMIT>0 && count_records>=LIMIT) break;
				} while(start<=length);
			
			}
		}
	
	private Fetcher fetch(String chrom,int length)
		throws Exception
		{
		return new Fetcher(chrom, length);
		}
	
	
	

	private int doWork()
		{
		try {
			DocumentBuilderFactory f=DocumentBuilderFactory.newInstance();
			f.setCoalescing(true);
			f.setNamespaceAware(true);
			f.setValidating(false);
			f.setExpandEntityReferences(true);
			f.setIgnoringComments(false);
			this.documentBuilder= f.newDocumentBuilder();
			
			TransformerFactory factory=TransformerFactory.newInstance();
			this.transformer=factory.newTransformer();
			this.transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
			
			List<Fetcher> fetchers=new ArrayList<Fetcher>(24);
			
			fetchers.add( fetch("1",249250621) );
			fetchers.add( fetch("2",243199373) );
			fetchers.add( fetch("3",198022430) );
			fetchers.add( fetch("4",191154276) );
			fetchers.add( fetch("5",180915260) );
			fetchers.add( fetch("6",171115067) );
			fetchers.add( fetch("7",159138663 ));		
			fetchers.add( fetch("8",146364022) );
			fetchers.add( fetch("9",141213431) );
			fetchers.add( fetch("10",135534747) );
			fetchers.add( fetch("11",135006516) );
			fetchers.add( fetch("12",133851895) );
			fetchers.add( fetch("13",115169878) );
			fetchers.add( fetch("14",107349540) );
			fetchers.add( fetch("15",102531392) );
			fetchers.add( fetch("16",90354753) );
			fetchers.add( fetch("17",81195210) );
			fetchers.add( fetch("18",78077248) );
			fetchers.add( fetch("19",59128983) );
			fetchers.add( fetch("20",63025520) );
			fetchers.add( fetch("21",48129895) );
			fetchers.add( fetch("22",51304566) );
			fetchers.add( fetch("X",155270560) );
			//fetch("Y",59373566); not in evs
			//fetch("M",16571);

			this.genome_total_size=0L;
			this.genome_curr_size=0L;
			for(Fetcher fetcher: fetchers)
				{
				this.genome_total_size += fetcher.length; 
				}
			
			for(Fetcher fetcher: fetchers)
				{
				fetcher.run();
				this.genome_curr_size += fetcher.length; 
				}
			
			} 
		catch (Exception e)
			{
			e.printStackTrace();
			return -1;
			}
		return 0;
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/EVS2Bed";
		}
	
	@Override
	public String getProgramName() {
		return "EVS2Bed";
		}
	
	@Override
	public String getProgramDescription() {
		return "Download data from EVS http://evs.gs.washington.edu/EVS as a BED chrom/start/end/XML For later use, see VCFTabixml.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-N (integer) download using a step of  'N' bases. Optional. Default:"+STEP_SIZE);
		out.println("-L (integer) limit to L records (for debugging). Optional. ");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"N:L:"))!=-1)
			{
			switch(c)
				{
				case 'N': this.STEP_SIZE=Math.max(1, Integer.parseInt(opt.getOptArg())); break;
				case 'L': this.LIMIT=  Integer.parseInt(opt.getOptArg()); break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		
		try
			{
			return doWork();
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	public static void main(String[] args)
		{
		new DumpExomeVariantServerData().instanceMainWithExit(args);
		}
	}
