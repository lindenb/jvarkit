/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.evs2bed;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.index.DynamicIndexCreator;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringReader;
import java.io.StringWriter;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.stream.StreamSource;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.LocationAwareOutputStream;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import edu.washington.gs.evs.SnpData;


public class EvsDumpXml
	extends AbstractCommandLineProgram
	{
	private int STEP_SIZE=1000000;
    private long LIMIT=-1L;
	private long count_records=0L;
	private long genome_total_size=0L;
	private long genome_curr_size=0L;
	private final int MAX_TRY=10; 
	private File outfilename = null;
	public static final String EVS_NS="http://webservice.evs.gs.washington.edu/";
	private XMLInputFactory xmlInputFactory;
	private Transformer transformer;
	private SortingCollection<String> sortingCollection;
	private boolean doSort=false;
	private LocationAwareOutputStream outputstream=null;
	
	
	private EvsDumpXml()
		{
		
		}
	
	private static class SnpDataBinding
		{
		private Marshaller marshaller;
		private Unmarshaller unmarshaller;
		protected SnpDataBinding()
			{
			try
				{
				JAXBContext jc = JAXBContext.newInstance(SnpData.class);
				this.marshaller = jc.createMarshaller();
				this.marshaller.setProperty(Marshaller.JAXB_ENCODING, "UTF-8");
				this.marshaller.setProperty(Marshaller.JAXB_FRAGMENT,true);
				this.unmarshaller= jc.createUnmarshaller();
				
				
				}
			catch(Exception err)
				{
				throw new RuntimeException(err);
				}
			}
		SnpData convert(String s)
			{
			try {
				return this.unmarshaller.unmarshal(new StreamSource(new StringReader(
						s)), SnpData.class).getValue();
				}
			catch (Exception err)
				{
				throw new RuntimeException(err);
				}
			}
		}
	
	
	private class SnpStringCodec 
	extends AbstractDataCodec<String>
		{
		@Override
		public String decode(DataInputStream dis) throws IOException
			{
			try
				{
				return IOUtils.readString(dis);
				}
			catch(java.io.EOFException err)
				{
				return null;
				}
			}
		@Override
		public void encode(DataOutputStream dos, String s) throws IOException {
			IOUtils.writeString(dos,s);
			}
		@Override
		public AbstractDataCodec<String> clone() {
			return new SnpStringCodec();
			}
		}

	private class SnpDataComparator
	extends SnpDataBinding
	implements Comparator<String>
		{
		
		
		@Override
		public int compare(final String s1,final  String s2)
			{
			SnpData o1=convert(s1);
			SnpData o2=convert(s2);
			int i=o1.getChromosome().compareTo(o2.getChromosome());
			if(i!=0) return i;
			i=o1.getChrPosition()-o2.getChrPosition();
			if(i!=0) return i;
			return 0;
			}
		}
	
	
	
	private void fetchEvsData(String chrom,int start,int end)
		{
		SnpDataBinding dataBinding=new SnpDataBinding();
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
		    		warning(
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
		    XMLEventReader xmlr=this.xmlInputFactory.createXMLEventReader(rd);
		 
		    while(xmlr.hasNext())
		    	{
		    	XMLEvent evt=xmlr.peek();
		    	if(!evt.isStartElement() || !evt.asStartElement().getName().getLocalPart().equals("snpList"))
		    		{
		    		xmlr.nextEvent();
		    		continue;
		    		}
		    	
		    	SnpData snpData=dataBinding.unmarshaller.unmarshal(xmlr,SnpData.class).getValue();
		    	StringWriter sw=new StringWriter();
		    	dataBinding.marshaller.marshal(new JAXBElement<SnpData>(new QName("snpList"), SnpData.class,snpData), sw);
				
		    	if(this.sortingCollection!=null)
		    		{
		    		this.sortingCollection.add(sw.toString());
		    		}
		    	else
		    		{
		    		this.outputstream.write(sw.toString().getBytes());
		    		this.outputstream.write('\n');
		    		}
				++count_records;
				if(LIMIT>0 && count_records>=LIMIT) break;
		    	}
		    xmlr.close();
		    wr.close();
		    rd.close();
			}
		catch(Exception err)
			{
			err.printStackTrace();
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
			if(EvsDumpXml.this.LIMIT>0 && EvsDumpXml.this.count_records>=EvsDumpXml.this.LIMIT) return;
			final int step=EvsDumpXml.this.STEP_SIZE;
			int start=1;
			do
				{
				fetchEvsData(chrom,start,start+step+1);
				
				
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
			this.xmlInputFactory =XMLInputFactory.newFactory();
			
			TransformerFactory factory=TransformerFactory.newInstance();
			this.transformer=factory.newTransformer();
			this.transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
			
			if(this.doSort)
				{
				this.sortingCollection=SortingCollection.newInstance(
						String.class,
						new SnpStringCodec(),
						new SnpDataComparator(),
						100000,
						getTmpDirectories()
						);
				this.sortingCollection.setDestructiveIteration(true);
				}
			
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
			
			DynamicIndexCreator indexer=null;
			if(this.outfilename!=null)
				{
				info("Opening "+this.outfilename);
				this.outputstream = new LocationAwareOutputStream(new FileOutputStream(this.outfilename));
				indexer =new DynamicIndexCreator(
						this.outfilename,
						 IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME
						);
				}
			else
				{
				this.outputstream = new LocationAwareOutputStream(System.out);
				}
			//print header
			final String xml_header = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"+
					"<evsData xmlns=\"http://webservice.evs.gs.washington.edu/\">\n"
					;
			this.outputstream.write(xml_header.getBytes());

			
			for(Fetcher fetcher: fetchers)
				{
				fetcher.run();
				this.genome_curr_size += fetcher.length; 
				}
			
			if(this.sortingCollection!=null)
				{
				SnpDataBinding snpDataBinding=new SnpDataBinding();
				this.sortingCollection.doneAdding();
				String prev = null;
				CloseableIterator<String> iter=sortingCollection.iterator();
				 while(iter.hasNext())
				 	{
					 String s= iter.next();
					 if( prev !=null && prev.equals(s))
					 	{
						continue; 
					 	}
					 
				    long position = outputstream.getPosition();
				    outputstream.write(s.getBytes());
				    outputstream.write('\n');//important SnpDataCodec needs separate lines
					if(indexer!=null)
						{
						 SnpData sd=snpDataBinding.convert(s);
						indexer.addFeature(new SnpDataFeature(sd),position);
						}
					 prev=s;
					 }
				 iter.close();
				}
			 long last_index=  this.outputstream.getPosition();
			 final String xml_footer = "</evsData>\n";
			 this.outputstream.write(xml_footer.getBytes());
			 this.outputstream.flush();
			 this.outputstream.close();

			 if(indexer!=null)
			 	{
				 info("Writing index");
				 final Index index = indexer.finalizeIndex(
						 last_index
						 );
				 index.writeBasedOnFeatureFile(this.outfilename);
			 	}
			} 
		catch (Exception e)
			{
			e.printStackTrace();
			return -1;
			}
		finally
			{
			if(this.sortingCollection!=null)
				this.sortingCollection.cleanup();
			}
		return 0;
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/Evs2Xml";
		}
	
	@Override
	public String getProgramName() {
		return "Evs2Xml";
		}
	
	@Override
	public String getProgramDescription() {
		return "Download data from EVS http://evs.gs.washington.edu/EVS as XML file. ";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-N (integer) download using a step of  'N' bases. Optional. Default:"+STEP_SIZE);
		out.println("-L (integer) limit to L records (for debugging). Optional. ");
		out.println("-o (filename) output filename. must end with"+SnpDataCodec.XML_FILE_SUFFIX+" will be indexed with tribble ");
		out.println("-s sort data. (default if filename specified with -o)");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"N:L:o:s"))!=-1)
			{
			switch(c)
				{
				case 'o': this.outfilename = new File(opt.getOptArg()); break;
				case 'N': this.STEP_SIZE=Math.max(1, Integer.parseInt(opt.getOptArg())); break;
				case 'L': this.LIMIT =  Integer.parseInt(opt.getOptArg()); break;
				case 's': this.doSort=true;break;
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
		
		if(this.outfilename!=null)
			{
			if(!this.outfilename.getName().endsWith(SnpDataCodec.XML_FILE_SUFFIX))
				{
				error("Ouput filename should end with "+SnpDataCodec.XML_FILE_SUFFIX+": "+outfilename);
				return -1;
				}
			this.addTmpDirectory(this.outfilename.getParentFile());
			this.doSort=true;
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
		/* https://blogs.oracle.com/joew/entry/jdk_7u45_aws_issue_123 */
		System.setProperty("jdk.xml.entityExpansionLimit","0");
		new EvsDumpXml().instanceMainWithExit(args);
		}
	}
