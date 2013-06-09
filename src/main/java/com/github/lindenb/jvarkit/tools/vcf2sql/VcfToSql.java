package com.github.lindenb.jvarkit.tools.vcf2sql;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.CommonInfo;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;

public class VcfToSql extends CommandLineProgram
	{
	
	private static final String COLUMN_ID="id INTEGER PRIMARY KEY AUTOINCREMENT,";

	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Creates the code to insert one or more VCF into a SQL database. ";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF files to process.",minElements=0)
	public List<File> IN=new ArrayList<File>();
    
    private PrintWriter out=new PrintWriter(System.out);
    
    @Override
    public String getVersion()
    	{
    	return "1.0";
    	}
    
	@Override
	protected int doWork()
		{
		try
			{
		    out.print( "create table if not exists VCFFILE("+
		    		COLUMN_ID+
		    		"filename TEXT NOT NULL"+
		    		")");
		    out.print( "create table if not exists SAMPLE("+
			    COLUMN_ID+
			    "name TEXT NOT NULL UNIQUE"+
			   ")");
		    out.print( "create table if not exists VARIATION("+
				COLUMN_ID+
				"CHROM VARCHAR(20) NOT NULL,"+
				"POS INT NOT NULL,"+
				"RS VARCHAR(50),"+
				"REF VARCHAR(10) NOT NULL"+
			    ")");
		    out.print(     "create table if not exists VCFROWINFO("+
				COLUMN_ID+
				"vcfrow_id INT NOT NULL REFERENCES VCFROW(id) ON DELETE CASCADE,"+
				"prop VARCHAR(10) NOT NULL,"+
				"value TEXT"+
				   ")");

			
		    
			out.println("begin transaction;");
			if(IN.isEmpty())
				{
				read(System.in,"<stdin>");
				}
			else
				{
				for(File input: IN)
					{
					InputStream in=IoUtil.openFileForReading(input);
					read(in,input.toString());
					in.close();
					}
				}
			out.println("commit;");
			}
		catch(IOException err)
			{
			err.printStackTrace();
			return -1;
			}
		return 0;
		}
	
	private void read(InputStream in,String filename)
		throws IOException
		{
		out.println("insert into FILE(filename) values ("+quote(filename)+");");
		AsciiLineReader r=new AsciiLineReader(in);
		

		VCFCodec codec=new VCFCodec();
		VCFHeader header=(VCFHeader)codec.readHeader(r);
		for(String S:header.getSampleNamesInOrder())
			{
			out.println("insert or ignore into SAMPLE(name) values ("+quote(S)+");");
			}
		String line;
		
		while((line=r.readLine())!=null)
			{
			out.println("/* "+line +" */");
			
			VariantContext var=codec.decode(line);
			
			//"create table if not exists FILE(id,filename text)";
			//"create table if not exists VARIATION(id,file_id,chrom,pos,start0,end0,rs_id,ref,qual)";
		
			
			out.println(
					"insert into VARIATION(file_id,chrom,pos,start0,end0,rs_id,ref,qual) values "+
					"(select max(id) from FILE),"+
					quote(var.getChr())+","+
					var.getStart()+","+
					(var.getStart()-1)+","+
					var.getEnd()+","+
					quote(var.getID())+","+
					quote(var.getReference().getDisplayString())+","+
					(var.getPhredScaledQual()<0?"NULL":var.getPhredScaledQual())+";"
					);
			//"create table if not exists ALT(id,var_id,alt)";

			for(Allele alt: var.getAlternateAlleles())
				{
				out.println(
						"insert into ALT(var_id,alt) values "+
						"(select max(id) from VARIATION),"+
						quote(alt.getDisplayString())+";"
						);
				}
			//"create table if not exists FILTER(id,var_id,filter)";

			for(String filter:var.getFilters())
				{
				out.println(
						"insert into FILTER(var_id,filter) values "+
						"(select max(id) from VARIATION),"+
						quote(filter)+";"
						);
				}
			CommonInfo infos=var.getCommonInfo();
			for(String key:infos.getAttributes().keySet())
				{
				Object val=infos.getAttribute(key);
				
				//"create table if not exists INFO(id,var_id,k,v)";

				
				out.println(
						"insert into INFO(var_id,k,v) values "+
						"(select max(id) from VARIATION),"+
						quote(key)+","+
						quote(infotoString(val))+";"
						);
				}
			GenotypesContext genotypesCtx =var.getGenotypes();
			for(Genotype g:genotypesCtx)
				{
				//"create table if not exists GENOTYPE(id,var_id,k,v)";

				
				out.println(
						"insert into GENOTYPE(var_id,sample_id,dp,ad,gq,is_phased) values "+
						"(select max(id) from VARIATION),"+
						"(select id from SAMPLE where name="+quote(g.getSampleName())+"),"+
						(g.hasDP()?"null":g.getDP())+","+
						quote(infotoString(g.getAD()))+","+
						g.getGQ()+","+
						(g.isPhased()?1:0)+","+
						(g.isHom()?1:0)+","+
						(g.isHomRef()?1:0)+","+
						(g.isHomVar()?1:0)+","+
						(g.isMixed()?1:0)+","+
						(g.isNoCall()?1:0)+","+
						(g.isNonInformative()?1:0)+","+
						(g.isAvailable()?1:0)+","+
						(g.isNoCall()?1:0)+","+
						(g.isCalled()?1:0)+","+
						(g.isPhased()?1:0)+","+						
						(g.isCalled()?1:0)+","+						
						g.getGenotypeString()+"/"+g.getDP()+"/"+g.getAlleles()
						);
				
				for(String key:g.getExtendedAttributes().keySet())
					{
					Object val=g.getExtendedAttribute(key);
					if(val==null) continue;
					out.println(
							"insert into GTPROP(genotype_id,k,v) values "+
							"(select max(id) from GENOTYPE),"+
							quote(key)+","+
							quote(infotoString(val))+";"
							);
					}
				if(!g.isAvailable()) break;
				if(!g.isCalled()) break;
				}
			
			}
		}
	
	private String quote(String s)
		{
		if(s==null) return "NULL";
		StringBuilder b=new StringBuilder();
		b.append("\'");
		for(int i=0;i<s.length();++i)
			{
			char c=s.charAt(i);
			switch(c)
				{
				case '\'': b.append("'''"); break;
				default: b.append(c); break;
				}
			}
		b.append("\'");
		return b.toString();
		}
	
	private String infotoString(Object o)
		{
		if(o instanceof List)
			{
			List<?> L=List.class.cast(o);
			StringBuilder b=new StringBuilder();
			for(int i=0;i< L.size();++i)
				{
				if(i>0) b.append(",");
				b.append(infotoString(L.get(i)));
				}
			return b.toString();
			}
		return o.toString();
		}
	public static void main(String[] args)
		{
		new VcfToSql().instanceMainWithExit(args);
		}
	}
