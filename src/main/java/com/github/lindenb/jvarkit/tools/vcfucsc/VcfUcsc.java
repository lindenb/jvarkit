package com.github.lindenb.jvarkit.tools.vcfucsc;

import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashSet;
import java.util.Set;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.bio.bin.BinIterator;
import com.github.lindenb.jvarkit.util.printf.Printf;

public class VcfUcsc extends AbstractVCFFilter2
	{
	private Printf expression=null;
	private String database="hg19";
	private String table=null;
	private Connection connection=null;
	private String jdbcuri="jdbc:mysql://genome-mysql.cse.ucsc.edu";
	private boolean has_bin_column=false;
	private String chromColumn=null;
	private String startColumn=null;
	private String endColumn=null;
	
	private void select(Set<String> atts,PreparedStatement pstmt) throws SQLException
		{
		ResultSet row=pstmt.executeQuery();
		while(row.next())
			{
			String s=this.expression.eval(row);
			if(s==null || s.isEmpty()) continue;
			atts.add(s);
			}
		row.close();
		}
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		String TAG="UCSC_"+database.toUpperCase()+"_"+table.toUpperCase();
		PreparedStatement pstmt=null;
		ResultSet row=null;
		VCFHeader header=in.getHeader();
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				TAG,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				database+"."+table)
				);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));

		
		out.writeHeader(h2);
		try	{
			StringBuilder b=new StringBuilder(
					"select * from "+database+"."+table+" where ");
			b.append(chromColumn).append("=? and NOT(");
			b.append(endColumn).append("<=? or ");
			b.append(startColumn).append(">=? ) ");
			if(has_bin_column) b.append(" and bin=?");
			
			pstmt=connection.prepareStatement(b.toString());
			
			while(in.hasNext())
				{
				VariantContext ctx=in.next();
				Set<String> atts=new HashSet<String>();
				pstmt.setString(1, ctx.getChr());
				pstmt.setInt(2, ctx.getStart());
				pstmt.setInt(3, ctx.getEnd());
				if(this.has_bin_column)
					{
					BinIterator biter=new BinIterator(ctx.getStart()-1, ctx.getEnd());
					while(biter.hasNext())
						{
						pstmt.setInt(4, biter.next());
						select(atts,pstmt);
						}
					}
				else
					{
					select(atts,pstmt);
					}
				if(atts.isEmpty())
					{
					out.add(ctx);
					continue;
					}
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(TAG,atts.toArray());
				out.add(vcb.make());
				}
			
			}
		catch(SQLException err)
			{
			throw new PicardException("SQLError", err);
			}
		finally
			{
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
			}
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "e:T:D:"))!=-1)
			{
			switch(c)
				{
				case 'D': this.database=opt.getOptArg();break;
				case 'T': this.table=opt.getOptArg();break;
				case 'e':
					{
					try
						{
						this.expression=Printf.compile(opt.getOptArg()); break;
						}
					catch(Exception err)
						{
						error(err);
						return -1;
						}
					}
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		if(this.table==null)
			{
			error("Table undefined.");
			return -1;
			}
		if(this.expression==null)
			{
			error("Expression undefined.");
			return -1;
			}
		try
			{
			info("Getting jdbc-driver");
			Class.forName("com.mysql.jdbc.Driver");
			this.connection=DriverManager.getConnection(
					jdbcuri+"/"+database+"?user=genome&password=");
			Statement stmt=this.connection.createStatement();
			ResultSet row=stmt.executeQuery("describe "+database+"."+table);
			Set<String> cols=new HashSet<String>();
			while(row.next())
				{
				String colName=row.getString("Field");
				cols.add(colName);
				}
			row.close();
			stmt.close();
			this.has_bin_column=cols.contains("bin");
			
			
			for(String col:new String[]{"chrom"})
				{
				if(this.chromColumn==null && cols.contains(col))
					{
					this.chromColumn=col;
					}
				}
			if(chromColumn==null)
				{
				error("cannot find chromColumn in "+cols);
				return -1;
				}
			
			for(String col:new String[]{"txStart","cdsStart","chromStart"})
				{
				if(this.startColumn==null && cols.contains(col))
					{
					this.startColumn=col;
					}
				}
			if(startColumn==null)
				{
				error("cannot find startColumn in "+cols);
				return -1;
				}
			for(String col:new String[]{"txEnd","cdsEnd","chromEnd"})
				{
				if(this.endColumn==null && cols.contains(col))
					{
					this.endColumn=col;
					}
				}
			if(endColumn==null)
				{
				error("cannot find endColumn in "+cols);
				return -1;
				}
			return doWork(opt.getOptInd(), args);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.connection);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfUcsc().instanceMainWithExit(args);
		}

	}
