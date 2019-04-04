/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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


*/
package com.github.lindenb.jvarkit.tools.ensembl;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.util.CloserUtil;
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
import com.github.lindenb.jvarkit.io.TeeInputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.utils.URIBuilder;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;


/**
BEGIN_DOC

## History

* rewritten 2018-02-19

## Example

the XML query:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "end" value = "10000000"/>
		<Filter name = "start" value = "1000000"/>
		<Filter name = "chromosome_name" value = "4"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "external_transcript_name" />
		<Attribute name = "transcription_start_site" />
		<Attribute name = "transcript_start" />
		<Attribute name = "transcript_end" />
	</Dataset>
</Query>
```

the XML file above contains three fields that will be used/replaced to query the position of a variant:

```xml
	<Filter name = "end" value = "10000000"/>
	<Filter name = "start" value = "1000000"/>
	<Filter name = "chromosome_name" value = "4"/>
```

running the query:

```xml
java  dist/vcfbiomart.jar -X query.xml input.vcf
```

```
##INFO=<ID=BIOMART,Number=.,Type=String,Description="Biomart query. Format: ensembl_gene_id|ensembl_transcript_id|start_position|end_position|external_transcript_name|transcription_start_site|transcript_start|transcript_end">
(...)
(...)A|||||4113|;BIOMART=ENSG00000183873|ENST00000414099|38589548|38691164|SCN5A-004|38674840|38589548|38674840,ENSG00000183873|ENST00000423572|38589548|38691164|SCN5A-003|38674853|38589553|38674853,ENSG00000183873|ENST00000413689|38589548|38691164|SCN5A-001|38691163|38589553|38691163,ENSG00000183873|ENST00000333535|38589548|38691164|SCN5A-014|38691119|38589557|38691119,ENSG00000183873|ENST00000455624|38589548|38691164|SCN5A-002|38674823|38590619|38674823,ENSG00000183873|ENST00000450102|38589548|38691164|SCN5A-008|38674840|38591459|38674840,ENSG00000183873|ENST00000449557|38589548|38691164|SCN5A-010|38674807|38591812|38674807,ENSG00000183873|ENST00000464652|38589548|38691164|SCN5A-011|38596040|38594472|38596040,ENSG00000183873|ENST00000491944|38589548|38691164|SCN5A-005|38691164|38655519|38691164,ENSG00000183873|ENST00000476683|38589548|38691164|SCN5A-013|38674711|38671333|38674711,ENSG00000183873|ENST00000327956|38589548|38691164|SCN5A-006|38687267|38674602|38687267,ENSG00000183873|ENST00000451551|38589548|38691164|SCN5A-203|38691164|38589553|38691164,ENSG00000183873|ENST00000443581|38589548|38691164|SCN5A-202|38691163|38589553|38691163,ENSG00000183873|ENST00000425664|38589548|38691164|SCN5A-201|38691163|38589553|38691163;(...)
(...)
```

END_DOC
 */
@Program(
		name="vcfbiomart",
		description="BiomartQueries with VCF",
		keywords={"vcf","ensembl","biomart","annotation"})
public class VcfBiomart extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBiomart.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names= {"-X","--xml"},description=" (XML-file) XML biomart template. <Query>...</Query>",required=true)
	private File xmlTemplate=null;
	@Parameter(names= {"-T","--tag"},description=" (string) VCF output tag.")
	private String TAG="BIOMART";
	@Parameter(names= {"-C","--contig","--chrom"},description="@name attribute <Filter> in the <Query> xml document.")
	private String chromColumn1= "chromosome_name";
		@Parameter(names= {"-S","--start"},description=" (int) column index (1-based) for start . Optional")
	private String startColumn1 = "start";
	@Parameter(names={"-E","--end"},description=" (int) column index (1-based) for end . Optional")
	private String endColumn1 = "end";
	@Parameter(names= {"-u","--url","--mart"},description=" (url) biomart service url. See http://grch37.ensembl.org/info/data/biomart/biomart_restful.html")
	private String serviceUrl="http://grch37.ensembl.org/biomart/martservice";
	@Parameter(names={"-tee","--tee"},description="'Tee' response to stderr")
	private boolean teeResponse = false;
	@Parameter(names={"-label","--label","--labels"},description="Add the field label in the INFO attribute 'label1|value1|label2|value2'")
	private boolean printLabels = false;

	private class FilterColumn
		{
		final String name;
		final Element element;
		FilterColumn(final String name,final Element element) {
			this.name= name;
			if(StringUtil.isBlank(name)) throw new IllegalArgumentException("empty name");
			this.element = element;
			}
		public void set(final String s)
			{
			this.element.setAttribute("value", s);
			}
		@Override
		public String toString() {
			return name;
			}
		}
	
	
	
	/** xml document <Query> .. </Query> ... */
	private Document domQuery = null;
	private FilterColumn filterColumnContig = null;
	private FilterColumn filterColumnStart = null;
	private FilterColumn filterColumnEnd = null;
	private List<String> attributes = new ArrayList<>();
	private CloseableHttpClient httpClient = null;

	
	private String escapeInfo(final String s)
		{
		if(StringUtil.isBlank(s)) return "";
		return VCFUtils.escapeInfoField(s);
		}

	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator iter,
		final VariantContextWriter out) 
		{
		HttpGet httpGet = null;
		final Pattern tab=Pattern.compile("[\t]");
		try
			{
			final TransformerFactory factory=TransformerFactory.newInstance();		
			final Transformer transformer = factory.newTransformer();
			//transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");

			final VCFHeader header= iter.getHeader();
			StringBuilder desc=new StringBuilder("Biomart query. Format: ");
			
			desc.append(this.attributes.stream().map(S->this.printLabels?S+"|"+S:S).collect(Collectors.joining("|")));
			
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			header.addMetaDataLine(new VCFInfoHeaderLine(
				this.TAG,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				desc.toString()
				));
			out.writeHeader(header);
			
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
			while(iter.hasNext())
				{
				final VariantContext ctx = progress.watch(iter.next());
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				vcb.rmAttribute(this.TAG);
				
				this.filterColumnContig.set(ctx.getContig());
				this.filterColumnStart.set(String.valueOf(ctx.getStart()));
				this.filterColumnEnd.set(String.valueOf(ctx.getEnd()));
				final StringWriter domToStr = new StringWriter();
				transformer.transform(
						new DOMSource(this.domQuery),
						new StreamResult(domToStr)
						);
								
				final URIBuilder builder = new URIBuilder(this.serviceUrl);
				builder.addParameter("query", domToStr.toString());
				
				//System.err.println("\nwget -O - 'http://grch37.ensembl.org/biomart/martservice?query="+escapedQuery+"'\n");
				//escapedQuery = URLEncoder.encode(escapedQuery,"UTF-8");
				
				httpGet = new HttpGet(builder.build());
				
				final CloseableHttpResponse httpResponse = httpClient.execute(httpGet);
				int responseCode = httpResponse.getStatusLine().getStatusCode();
				if(responseCode != 200)
				 	{
					throw new RuntimeIOException("Response code was not 200. Detected response was "+responseCode);
				 	}
				 InputStream response =httpResponse.getEntity().getContent();
				 if(this.teeResponse)
					 {
					 response = new TeeInputStream(response,stderr(),false);
					 }
				 final BufferedReader br = new BufferedReader(new InputStreamReader(response));
				 final Set<String> infoAtts = 
						 br.
						 lines().
						 filter(L->!StringUtil.isBlank(L)).
						 filter(L->!L.equals("[success]")).
						 map(L->tab.split(L)).
						 map(T->{
							final StringBuilder sb=new StringBuilder();
							for(int i=0; i< this.attributes.size();i++)
								{
								if(i>0) sb.append("|");
								if(this.printLabels) sb.append(escapeInfo(this.attributes.get(i))).append("|");
								sb.append(i< T.length?escapeInfo(T[i]):"");
								}	
							return sb.toString();
						 }).
						 collect(Collectors.toCollection(LinkedHashSet::new));
				CloserUtil.close(br);
				CloserUtil.close(response);
				CloserUtil.close(httpResponse);
				if(!infoAtts.isEmpty())
					{
					vcb.attribute(this.TAG, new ArrayList<>(infoAtts));
					}
				out.add(vcb.make());
				}
			progress.finish();
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			throw new RuntimeIOException(err);
			}
	
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.xmlTemplate==null)
			{
			LOG.error("Undefined XML template");
			return -1;
			}
		try
			{	
			final DocumentBuilderFactory f=DocumentBuilderFactory.newInstance();
			f.setCoalescing(true);
			f.setNamespaceAware(false);
			f.setValidating(false);
			f.setExpandEntityReferences(true);
			f.setIgnoringComments(false);
			f.setIgnoringElementContentWhitespace(true);
			final DocumentBuilder docBuilder= f.newDocumentBuilder();
			LOG.info("Parsing xml "+xmlTemplate);
			this.domQuery = docBuilder.parse(this.xmlTemplate);

			final Element queryElement=this.domQuery.getDocumentElement();
			if(queryElement==null || !queryElement.getNodeName().equals("Query"))
				{
				LOG.error("XML root is not <Query/> but "+ queryElement.getNodeName());
				return -1;
				}
			queryElement.setAttribute("formatter", "TSV");
			queryElement.setAttribute("header", "0");
			queryElement.setAttribute("uniqueRows", "1");
			queryElement.setAttribute("count", "");
			/*
			*"If you want to make sure you are getting all the data from your BioMart query, you can add a "CompletionStamp" to the xml file. To do this, just open the previously obtained xml file in the "Obtaining the BioMart xml" section and add the following text in the query tag"
			 */
			queryElement.setAttribute("completionStamp","1");
			
			Element dataSetElement=null;
			for(Node c1= queryElement.getFirstChild();
					c1!=null;c1=c1.getNextSibling())
				{
				if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
				if(c1.getNodeName().equals("Dataset"))
					{
					if(dataSetElement!=null)
						{
						LOG.error("XML: Two <DataSet> elements ?");
						return -1;
						}
					dataSetElement=Element.class.cast(c1);
					for(Node c2= dataSetElement.getFirstChild();c2!=null;c2=c2.getNextSibling())
						{
						if(c2.getNodeType()!=Node.ELEMENT_NODE) continue;
						final Element e2=Element.class.cast(c2);
						if(e2.getNodeName().equals("Attribute"))
							{
							final Attr att;
							if((att=e2.getAttributeNode("name"))==null)
								{
								LOG.error("XML: <Attribute> without @name ?");
								return -1;
								}
							this.attributes.add(att.getValue());
							}
						else if(c2.getNodeName().equals("Filter"))
							{
							final Attr att;
							if((att=e2.getAttributeNode("name"))==null)
								{
								LOG.error("XML: <Filter>  without @name ?");
								return -1;
								}
							final String attName= att.getValue();
							if(attName.equals(this.chromColumn1))
								{
								this.filterColumnContig = new FilterColumn(this.chromColumn1,e2);
								}
							else if(attName.equals(this.startColumn1))
								{
								this.filterColumnStart = new FilterColumn(this.startColumn1,e2);
								}
							else if(attName.equals(this.endColumn1))
								{
								this.filterColumnEnd = new FilterColumn(this.endColumn1,e2);
								}
							}
						}
					}
				
				}
			if(this.filterColumnContig==null)
				{
				LOG.error("<Filter @name='"+this.chromColumn1+"'> is missing");
				return -1;
				}
			if(this.filterColumnStart==null)
				{
				LOG.error("<Filter  @name='"+this.startColumn1+"'> is missing");
				return -1;
				}
			if(this.filterColumnEnd==null)
				{
				LOG.error("<Filter @name='"+this.endColumn1+"'> is missing");
				return -1;
				}
			
			if(this.attributes.isEmpty())
				{
				LOG.error("no <Attribute>");
				return -1;
				}
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		try
			{
			/** create http client */
			this.httpClient = HttpClients.createSystem();
			
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.httpClient);
			this.httpClient=null;
			}
		
		}
	
	public static void main(final String[] args)
		{
		new VcfBiomart().instanceMainWithExit(args);
		}

	}
