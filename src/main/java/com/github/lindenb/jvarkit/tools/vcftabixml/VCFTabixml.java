package com.github.lindenb.jvarkit.tools.vcftabixml;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.transform.OutputKeys;
import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.TabixReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;


public class VCFTabixml extends AbstractVCFFilter
	{
	/*
	 * 				System.out.println("Usage: java -jar vcftabix.jar -f src.vcf.gz -x style.xsl (file.vcf|stdin) " );
				System.out.println(" -f (BED indexed with tabix. The 4th column is a XML string.) REQUIRED.");
				System.out.println(" -H (String like '##INFO=...') append extra-info header");
				System.out.println(" -x xslt-stylesheet. REQUIRED. Should produce a valid set of INFO field.");
				System.out.println("4th column of the BED indexed with TABIX is a XML string." +
						"It will be processed with the xslt-stylesheet and should procuce a valdid set of" +
						" INFO fields. Carriage returns will be removed." +
						"Parameters to be passed to the stylesheet: vcfchrom (string) vcfpos(int) vcfref(string) vcfalt(string). ");

	 */
	public String BEDFILE=null;
	public File STYLESHEET=null;
	private Set<String> extraInfoFields=new LinkedHashSet<String>();
	private Templates stylesheet;
	
	private VCFTabixml()
		{
		

		}
	
	private boolean isEmpty(String s)
		{
		return s==null || s.equals(".") || s.isEmpty();
		}
	
	@Override
	protected void doWork(LineReader in, VariantContextWriter w)
			throws IOException
		{
		TabixReader tabixReader =null;

		try {
			tabixReader=new TabixReader(this.BEDFILE);
			
			Pattern tab=Pattern.compile("[\t]");
			Transformer transformer=this.stylesheet.newTransformer();
			transformer.setOutputProperty(OutputKeys.METHOD,"text");

			VCFCodec codeIn=new VCFCodec();		
			VCFHeader header=(VCFHeader)codeIn.readHeader(in);
			String line;
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			//TODO h2.addMetaDataLine(new VCFInfoHeaderLine(TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "metadata added from "+TABIXFILE+" . Format was "+FORMAT));
			w.writeHeader(h2);
			
			while((line=in.readLine())!=null)
				{
				VariantContext ctx=codeIn.decode(line);
				
				

				TabixReader.Iterator iter= tabixReader.query(
						ctx.getChr()+":"+(ctx.getStart())+"-"+(ctx.getEnd()+1));
				String line2=null;
				
				String infoToAppend=null;
				while(iter!=null && (line2=iter.next())!=null)
					{

					String tokens2[]=tab.split(line2,5);
					
					if(tokens2.length<4)
						{
						System.err.println("[VCFTabixml] VCF. Error not enough columns in tabix.line "+line2);
						return;
						}
					
					int chromStart=Integer.parseInt(tokens2[1]);
					int chromEnd=Integer.parseInt(tokens2[2]);
					if(chromStart+1!=chromEnd)
						{
						System.err.println("Error in "+this.BEDFILE+" extected start+1=end int "+tokens2[0]+":"+tokens2[1]+"-"+tokens2[2]);
						continue;
						}
					
					

					if(ctx.getStart()-1!=chromStart) continue;

					
					transformer.setParameter("vcfchrom",ctx.getChr());
					transformer.setParameter("vcfpos",ctx.getStart());
					transformer.setParameter("vcfref",ctx.getReference().getBaseString());
					transformer.setParameter("vcfalt",ctx.getAltAlleleWithHighestAlleleCount().getBaseString());
					
					StringWriter sw=new StringWriter();
					try {
						StreamSource src=new StreamSource(new StringReader(tokens2[3]));
						StreamResult rez=new StreamResult(sw);
						transformer.transform(src, rez);
						}
					catch (Exception e)
						{
						continue;
						}
					
					infoToAppend=sw.toString().replace('\n',';').replaceAll("[;]+",";");
					if(infoToAppend.isEmpty() || infoToAppend.equals(";"))
						{
						infoToAppend=null;
						}
					
					}
				/*
				String newInfo;
				if(isEmpty(tokens[7]))
					{
					newInfo=infoToAppend;
					if(isEmpty(newInfo)) newInfo=".";
					}
				else
					{
					newInfo=tokens[7];
					if(!newInfo.endsWith(";")) newInfo+=";";
					if(!isEmpty(infoToAppend)) newInfo+=infoToAppend;
					}*/

				}
			
			
		} catch (TransformerConfigurationException err)
			{
			throw new IOException(err);
			}
		}
	

	public static void main(String[] args) throws Exception
		{
		new VCFTabixml().instanceMainWithExit(args);
		}
	}
