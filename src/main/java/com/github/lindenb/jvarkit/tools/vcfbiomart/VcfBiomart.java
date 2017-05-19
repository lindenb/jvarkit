package com.github.lindenb.jvarkit.tools.vcfbiomart;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
BEGIN_DOC

##Example

Let's annotate a VCF with the gene-ncbi/gene-start/gene-end/gene-ncbiGeneId .

From Ensembl, we've downloaded the following XML file:

```xml
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "entrezgene" />
	</Dataset>
</Query>
```
This file is used  as follow:

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  sed 's/chr//g' |\
 java -jar dist/vcfbiomart.jar -X mart.xml -T NCBIGENE  -C 1 -S 2 -E 3 |\
 grep NCBIGENE

##INFO=<ID=NCBIGENE,Number=.,Type=String,Description="Biomart query. Format:chromosome_name|start_position|end_position|entrezgene">
1	145273345	.	T	C	289.85	.	NCBIGENE=1|145209119|145291972|388677,1|145209145|145319796|	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
10	1142208	.	T	C	3404.30	.	NCBIGENE=10|1095478|1178237|22884	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
10	135210791	.	T	C	65.41	.	NCBIGENE=10|135204338|135233999|,10|135207598|135234811|92170	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
```
###XML attributes: @chrom/@start/@end.
By adding the attributes chrom=true, start=true end=true the column for chrom/start/end can be specified in the XML
```XML
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "chromosome_name" chrom="true"/>
		<Attribute name = "start_position" start="true"/>
		<Attribute name = "end_position" end="true" />
		<Attribute name = "entrezgene" />
	</Dataset>
</Query>
```
```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  sed 's/chr//g' |\
 java -jar dist/vcfbiomart.jar -X mart.xml -T NCBIGENE  |\
 grep NCBIGENE

##INFO=<ID=NCBIGENE,Number=.,Type=String,Description="Biomart query. Format:chromosome_name|start_position|end_position|entrezgene">
1	145273345	.	T	C	289.85	.	NCBIGENE=1|145209119|145291972|388677,1|145209145|145319796|	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
10	1142208	.	T	C	3404.30	.	NCBIGENE=10|1095478|1178237|22884	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
10	135210791	.	T	C	65.41	.	NCBIGENE=10|135204338|135233999|,10|135207598|135234811|92170	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
```
###Hiding columns: @visible=false
By adding the attribute visible=false, some columns can be removed from the result.
```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "chromosome_name" chrom="true" visible="false"/>
		<Attribute name = "start_position" start="true" visible="false"/>
		<Attribute name = "end_position" end="true" visible="false"/>
		<Attribute name = "entrezgene" />
	</Dataset>
</Query>
```

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  sed 's/chr//g' |\
 java -jar dist/vcfbiomart.jar -X mart.xml -T NCBIGENE  |\
 grep NCBIGENE

##INFO=<ID=NCBIGENE,Number=.,Type=String,Description="Biomart query. Format:entrezgene">
1	145273345	.	T	C	289.85	.	NCBIGENE=388677	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
10	1142208	.	T	C	3404.30	.	NCBIGENE=22884	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
10	135210791	.	T	C	65.41	.	NCBIGENE=92170	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
```
if all the attributes are set to visible=false, the INFO is set to 'Flag'
```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "chromosome_name" chrom="true" visible="false"/>
		<Attribute name = "start_position" start="true" visible="false"/>
		<Attribute name = "end_position" end="true" visible="false"/>
		<Attribute name = "entrezgene" visible="false"/>
	</Dataset>
</Query>
```

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  sed 's/chr//g' |\
 java -jar dist/vcfbiomart.jar -X mart.xml -T NCBIGENE  |\
 grep NCBIGENE

##INFO=<ID=NCBIGENE,Number=0,Type=Flag,Description="Biomart query.">
1	145273345	.	T	C	289.85	.	NCBIGENE	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
10	1142208	.	T	C	3404.30	.	NCBIGENE	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
10	135210791	.	T	C	65.41	.	NCBIGENE	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
```
END_DOC
 */
@Program(name="vcfbiomart",description="BiomartQueries with VCF",keywords={"vcf","ensembl","biomart"})
public class VcfBiomart extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBiomart.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	private List<Element> attributes=new ArrayList<Element>();
	private Document dom=null;
	private Element queryElement=null;
	private Element dataSetElement=null;
	private Element filterElementChromosomalLocation=null;
	private Set<Integer> visibleIndexes0=new HashSet<Integer>();
	
	
	private String escapeInfo(String s)
		{
		if(s==null || s.isEmpty()) return "";
		return s.replaceAll("[ =;\t]","_").replaceAll("[_]+", " ");
		}

	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		try {
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
			LOG.error(err);
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
							ctx.getContig()+":"+ctx.getStart()+":"+ctx.getEnd()+":1"	
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
						LOG.error(e);
						throw new IOException(e);
						}
					
					LOG.info("POSTing to "+this.serviceUrl+" buffer.size="+buffer.size());
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
				     LOG.debug(q);
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
					@SuppressWarnings("resource")
					LineIterator li=new  LineIteratorImpl(r);
					while(li.hasNext())
						{
						String line=li.next();
						LOG.debug(line);
						String tokens[]=tab.split(line);
						LOG.debug(line+" L="+tokens.length );
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
						List<String> array=new ArrayList<String>(new HashSet<String>(treemap.getOverlapping(new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd()))));
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
		return 0;
		}
		catch(Exception err) {
			LOG.error(err);
			return -1;
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
				LOG.info("Attribute @"+tag+" was specified in the XML");
				if(column!=-1)
					{
					LOG.error("XML: Two @"+tag+"=true ?");
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
			LOG.error("Cannot use column for \""+tag+"\".");
			return -1;
			}
		return column;
		}
	
	@Parameter(names="-X",description=" (XML-file) XML biomart template.")
	String xmlTemplate=null;
	@Parameter(names="-n",description=" (int) batch size.);")
	private int batchSize=200;

	@Parameter(names="-T",description=" (string) VCF output tag.")
	private String TAG="BIOMART";

	@Parameter(names="-C",description=" (int) column index (1-based) for chromosome . Optional")
	private int chromColumn1=-1;
	

	@Parameter(names="-S",description=" (int) column index (1-based) for start . Optional")
	private int startColumn1=-1;
	@Parameter(names="-E",description=" (int) column index (1-based) for end . Optional")
	private int endColumn1=-1;
	@Parameter(names="-u",description=" (url) biomart service url")
	private String serviceUrl="http://www.biomart.org/biomart/martservice/result";

	
	@Override
	public int doWork(List<String> args) {
		
		if(xmlTemplate==null)
			{
			LOG.error("Undefined XML template");
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
			LOG.info("Parsing xml "+xmlTemplate);
			InputStream in=IOUtils.openURIForReading(xmlTemplate);
			this.dom=docBuilder.parse(in);
			in.close();
			this.queryElement=this.dom.getDocumentElement();
			if(this.queryElement==null || !this.queryElement.getNodeName().equals("Query"))
				{
				LOG.error("XML root is not <Query/> but "+this.queryElement.getNodeName());
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
						LOG.error("XML: Two <DataSet> elements ?");
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
								LOG.error("XML: Attribute without @name ?");
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
									LOG.error("XML: two Filter@chromosomal_region ?");
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
				LOG.error("DataSet element missing");
				return -1;
				}
			if(this.attributes.size()<2)
				{
				LOG.error("Expected at least two elementss <Attribute>");
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
				LOG.error("CHROM column index out of range");
				return -1;
				}
			if(startColumn1==-1)
				{
				this.startColumn1=findColumn1("start");
				if(this.startColumn1<1) return -1;
				}
			if(this.startColumn1< 1 || this.startColumn1> attributes.size())
				{
				LOG.error("START column index out of range");
				return -1;
				}
			if(endColumn1==-1)
				{
				this.endColumn1=findColumn1("end");
				if(this.startColumn1<1) endColumn1=startColumn1;
				}
			if(this.endColumn1< 1 || this.endColumn1> attributes.size())
				{
				LOG.error("END column index out of range");
				return -1;
				}
			
			
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		
		for(int i=0;i< this.attributes.size();++i)
			{
			Element E=this.attributes.get(i);
			Attr att=E.getAttributeNode("visible");
			if(att!=null && att.getValue().equals("false")) continue;
			this.visibleIndexes0.add(i);
			}
		
		return doVcfToVcf(args, outputFile);
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfBiomart().instanceMainWithExit(args);
		}

	}
