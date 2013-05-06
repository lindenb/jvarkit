package fr.inserm.umr1087.jvarkit.tools.vcftabixml;

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


import net.sf.samtools.tabix.TabixReader;


public class VCFTabixml
	{
	private PrintStream out=System.out;
	private String bedFile;
	private TabixReader tabixReader =null;
	private Set<String> extraInfoFields=new LinkedHashSet<String>();
	private Templates stylesheet;
	
	private VCFTabixml()
		{
		

		}
	
	private boolean isEmpty(String s)
		{
		return s==null || s.equals(".") || s.isEmpty();
		}
	
	
	private void run(BufferedReader in) throws IOException,TransformerConfigurationException
		{
		Transformer transformer=this.stylesheet.newTransformer();
		transformer.setOutputProperty(OutputKeys.METHOD,"text");
		Pattern tab=Pattern.compile("[\t]");
		String line;
		while((line=in.readLine())!=null)
			{

			if(line.isEmpty()) continue;
			
			if(line.startsWith("#"))
				{
				if(line.startsWith("#CHROM"))
					{
					out.println("##Annotated with "+getClass().getCanonicalName()+":"+bedFile);
					for(String h:this.extraInfoFields)
						{
						if(!h.startsWith("##")) out.print("##");
						out.println(h);
						}
					out.println(line);
					continue;
					}
				out.println(line);
				continue;
				}
			
			String tokens[]=tab.split(line,9);
			if(tokens.length<8)
				{
				System.err.println("[VCFTABIXml] Error not enough columns in "+line);
				continue;
				}
			String chrom=tokens[0];
			Integer pos1=Integer.parseInt(tokens[1]);
			
			
			TabixReader.Iterator iter=this.tabixReader.query(chrom+":"+(pos1)+"-"+(pos1+1));
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
					System.err.println("Error in "+this.bedFile+" extected start+1=end int "+tokens2[0]+":"+tokens2[1]+"-"+tokens2[2]);
					continue;
					}
				
				

				if(pos1-1!=chromStart) continue;

				
				transformer.setParameter("vcfchrom",tokens[0]);
				transformer.setParameter("vcfpos",pos1);
				transformer.setParameter("vcfref",tokens[3].toUpperCase());
				transformer.setParameter("vcfalt",tokens[4].toUpperCase());
				
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
				}
			
			
			for(int i=0;i< tokens.length;++i)
				{
				if(i>0) out.print('\t');
				out.print(i==7?newInfo:tokens[i]);
				}
			out.println();
			}
		
		}
	public int run(String[] args) throws Exception
		{
		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				System.out.println("VCF Tabixml. Author: Pierre Lindenbaum PhD. 2013.");
				System.out.println("Usage: java -jar vcftabix.jar -f src.vcf.gz -x style.xsl (file.vcf|stdin) " );
				System.out.println(" -f (BED indexed with tabix. The 4th column is a XML string.) REQUIRED.");
				System.out.println(" -H (String like '##INFO=...') append extra-info header");
				System.out.println(" -x xslt-stylesheet. REQUIRED. Should produce a valid set of INFO field.");
				System.out.println("4th column of the BED indexed with TABIX is a XML string." +
						"It will be processed with the xslt-stylesheet and should procuce a valdid set of" +
						" INFO fields. Carriage returns will be removed." +
						"Parameters to be passed to the stylesheet: vcfchrom (string) vcfpos(int) vcfref(string) vcfalt(string). ");
				return 0;
				}
			else if(args[optind].equals("-H")&& optind+1< args.length)
				{
				this.extraInfoFields.add(args[++optind]);
				}
			else if(args[optind].equals("-x")&& optind+1< args.length)
				{
				TransformerFactory tFactory =  TransformerFactory.newInstance();
				this.stylesheet=tFactory.newTemplates(new StreamSource(new File(args[++optind])));
				}			
			else if(args[optind].equals("-f") && optind+1< args.length)
				{
				this.bedFile=args[++optind];
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unnown option: "+args[optind]);
				return -1;
				}
			else
				{
				break;
				}
			++optind;
			}
		if(bedFile==null)
			{
			System.err.println("Undefined bed+xml tabix File");
			return -1;
			}
		if(this.stylesheet==null)
			{
			System.err.println("No stylesheet provided");
			return -1;
			}
		this.tabixReader=new TabixReader(this.bedFile);
		
		if(optind==args.length)
			{
			this.run(new BufferedReader(new InputStreamReader(System.in)));
			}
		else if(optind+1==args.length)
			{
			String inputName=args[optind++];
			BufferedReader in=new BufferedReader(new FileReader(inputName));
			this.run(in);
			in.close();
			}
		else
			{
			System.err.println("Illegal Number of arguments");
			return -1;
			}
		return 0;
		}
	
	public static void main(String[] args) throws Exception
		{
		VCFTabixml app=new VCFTabixml();
		int ret=app.run(args);
		if(ret!=0) System.exit(ret);
		}
}
