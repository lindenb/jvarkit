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
* 2014 creation
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.Reader;
import java.sql.Clob;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;


public class VcfDerby01
	extends AbstractVcfDerby01
	{
	private long ID_GENERATOR = System.currentTimeMillis();
	public static final String VCF_HEADER_FILE_IDENTIFIER="VcfFileIdentifier";
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfDerby01.class);
	private Connection conn=null;
	public VcfDerby01()
		{
		}
	 
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.derbyFilePath==null ||super.derbyFilePath.isEmpty())
			{
			return wrapException("Undefined Derby DB Directory: option -"+OPTION_DERBYFILEPATH);
			}
		try {
			Class.forName("org.apache.derby.jdbc.ClientDriver").newInstance();
		} catch(Exception err ){
			LOG.error("Cannot get derby driver");
			return wrapException(err);
		}
		return super.initializeKnime();
	 	}
	
	
	private File getDerbyDirectory() {
		return new File(this.derbyFilePath);
	}
	
	private void openDerby() {
		try {
			boolean create;
			final Properties props = new Properties();
			final File derbyDir = getDerbyDirectory();
			LOG.info("open derby :" + getDerbyDirectory());
			if(derbyDir.exists()) {
				if(!derbyDir.isDirectory()) {
					throw new RuntimeIOException("derby database is not a directory : "+derbyDir);
					}
				create=false;
				}
			else
				{
				create=true;
				}
			props.setProperty("create", String.valueOf(create));
			this.conn = DriverManager.getConnection("jdbc:derby:"+derbyDir,props);
			
			if(create) {
				final String tableId = "ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY";
				final Statement stmt= this.conn.createStatement();
				final String sqls[]={
						"CREATE TABLE ROWCONTENT("+tableId+",MD5SUM CHAR(32) UNIQUE,CONTENT CLOB,CONTIG VARCHAR(20),START INT,STOP INT,ALLELE_REF VARCHAR(50))",
						"CREATE TABLE VCF("+tableId+",NAME VARCHAR(255))",
						"CREATE TABLE VCFROW("+tableId+",VCF_ID INTEGER CONSTRAINT row2vcf REFERENCES VCF,ROW_ID INTEGER CONSTRAINT row2content REFERENCES ROWCONTENT)"
						};
				for(final String sql:sqls) {
					LOG.warn(sql);
					stmt.execute(sql);
				}
				stmt.close();
			}
			this.conn.setAutoCommit(true);
		} catch (Exception e) {
			CloserUtil.close(this.conn);
			this.conn = null;
			throw new RuntimeException(e);
		}
	}
	
	private static long getLastGeneratedId(final PreparedStatement pstmt) throws SQLException {
			ResultSet keys = null;
			long id = -1L;
			try {
			keys = pstmt.getGeneratedKeys();
		    while (keys.next()) {
		       if(id!=-1L) throw new IllegalStateException();
		       id= keys.getLong(1);
		       }	 
			if(id==-1L) throw new IllegalStateException("No SQL ID was found");
			return id;
			}
			finally {
				CloserUtil.close(keys);
			}
		}
	
	private void closeDerby() {
		CloserUtil.close(this.conn);
		this.conn = null;
		try {
			final Properties props = new Properties();
			props.setProperty("shutdown", "true");
			final File derbyDir = getDerbyDirectory();
			DriverManager.getConnection("jdbc:derby:"+derbyDir,props);
		} catch (Exception e) {
			
			}
		}
	private Collection<Throwable> doCommandDumpAll(){
		Statement stmt = null;
		ResultSet row = null;
		final Set<Long> vcfids = new TreeSet<>();
		if(!getInputFiles().isEmpty()) {
			return wrapException("Too many arguments");
		}
		try {
			stmt = this.conn.createStatement();
			row = stmt.executeQuery("SELECT ID from VCF");
			while(row.next()) {
				vcfids.add(row.getLong(1));
			}
			row.close();row=null;
			stmt.close();stmt=null;
			return dump(vcfids);
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(row);
			CloserUtil.close(stmt);
		}
	}

	private Collection<Throwable> dump(final Set<Long> vcfIds){
		PreparedStatement pstmt = null;
		PreparedStatement pstmt2 = null;
		ResultSet row = null;
		PrintWriter pwOut = null;
		FileOutputStream fout = null;
		ZipOutputStream zout= null;
		int num_vcf_exported=0;
		Set<String> seen_names=null;
		try {
			pstmt = this.conn.prepareStatement("SELECT NAME from VCF where ID=?");
			pstmt2 = this.conn.prepareStatement("SELECT ROWCONTENT.CONTENT FROM VCF,VCFROW,ROWCONTENT WHERE VCFROW.VCF_ID=VCF.ID AND VCFROW.ROW_ID = ROWCONTENT.ID AND VCF.ID=? ORDER BY VCFROW.ID ");
			
			if(getOutputFile()!=null && getOutputFile().getName().endsWith(".zip")) {
				fout = new FileOutputStream(getOutputFile());
				zout =new ZipOutputStream(fout);
				seen_names=new HashSet<>();
			} else
			{
			pwOut = openFileOrStdoutAsPrintWriter();
			}
			
			for(final long vcf_id : vcfIds)
				{
				LOG.info("Getting VCF "+vcf_id);
				
				pstmt.setLong(1, vcf_id);
				String vcfName=null;
				row = pstmt.executeQuery();
				while(row.next()) {
					vcfName = row.getString(1);
				}
				row.close();
				if(vcfName==null) 
					{
					return wrapException("Cannot find VCF ID "+vcf_id);
					}
				LOG.info("dumping "+vcfName+" ID:"+vcf_id);
				pstmt2.setLong(1, vcf_id);
				row =  pstmt2.executeQuery();
				
				ZipEntry entry=null;
				if(zout!=null)
					{
					int t=0;
					String entryName=vcfName+(vcfName.endsWith(".vcf.")?"":".vcf");
					while(seen_names.contains(entryName)) {
						entryName=String.format("%05d.ID%d.vcf",++t, vcf_id);
					}
					seen_names.add(entryName);
					entry = new ZipEntry(entryName);
					zout.putNextEntry(entry);
					pwOut = new PrintWriter(zout);
					}
				
				while(row.next()) {
					final Clob clob = row.getClob(1);
					final Reader r=clob.getCharacterStream();
					IOUtils.copyTo(r,pwOut);
					r.close();
					pwOut.println();
				}
				
				if(zout!=null) {
					pwOut.flush();
					pwOut=null;
					zout.closeEntry();
				}
				
				row.close();
				num_vcf_exported++;
				}
				
			pstmt.close();
			pstmt2.close();
			if(zout!=null){
				zout.finish();
				zout.flush();
				fout.close();
			} else
			{
				pwOut.flush();
				pwOut.close();
			}
			if(num_vcf_exported==0 && zout==null) {
				return wrapException("NO VCF WAS EXPORTED");
			}
			LOG.info("count(VCF) exported "+num_vcf_exported);
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(fout);
			CloserUtil.close(zout);
			CloserUtil.close(pwOut);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
			CloserUtil.close(pstmt2);
		}
	}

	
	
	private Collection<Throwable> doCommandDump(){
		final Pattern comma = Pattern.compile("[,]");
		final List<String> args = getInputFiles();
		final Set<Long> vcfIds = new TreeSet<>();
		try {
			for(final String token: args) {
				for(final String idStr: comma.split(token))
					{
					if(idStr.trim().isEmpty()) continue;
					LOG.info("Getting VCF "+idStr);
					long vcf_id = -1L;
					try {
						vcf_id = Long.parseLong(idStr);
					} catch(final NumberFormatException err) {
						return wrapException("Bad VCF ID :"+idStr);
					}
					vcfIds.add(vcf_id);
					}
				}
			return dump(vcfIds);
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
		}
	}

	
	private Collection<Throwable> doReadConcatenatedVcf(){
		PreparedStatement pstmt = null;
		PreparedStatement pstmt2 = null;
		PreparedStatement pstmt3 = null;
		ResultSet row = null;
		PrintWriter pw = null;
		final List<String> args = new ArrayList<>(IOUtils.unrollFiles(getInputFiles()));
		LOG.info(args.toString());
		LineIterator lineIter=null;
		try {
			int fileidx=0;
			
			pw = openFileOrStdoutAsPrintWriter();
			pw.println("#ID\tNAME");

			do
			{
				if(fileidx==0 && args.isEmpty()) {
					lineIter = IOUtils.openStreamForLineIterator(stdin());
				} else
				{
					lineIter = IOUtils.openURIForLineIterator(args.get(fileidx));
				}
				int num_vcf_in_this_stream = 0;
				while(lineIter.hasNext()) {
					String filename= "vcf"+(++ID_GENERATOR);
					if(num_vcf_in_this_stream==0 && !args.isEmpty()) {
						filename = args.get(fileidx);
					}
					
					final List<String> headerLines = new ArrayList<>();
					while(lineIter.hasNext() && lineIter.peek().startsWith("#")) {
						headerLines.add(lineIter.next());
					}
					final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(headerLines);
					final VCFHeaderLine nameLine=cah.header.getOtherHeaderLine(VCF_HEADER_FILE_IDENTIFIER);
					if(nameLine!=null && !nameLine.getValue().trim().isEmpty()) {
						filename = nameLine.getValue().trim();
					}
					pstmt = this.conn.prepareStatement("INSERT INTO VCF(NAME) VALUES(?)",PreparedStatement.RETURN_GENERATED_KEYS);
					pstmt.setString(1, filename);
					if(pstmt.executeUpdate()!=1) {
						return wrapException("Cannot insert VCF ?");
					}
					final long vcf_id =getLastGeneratedId(pstmt);
					pstmt.close();
					
					pw.print(vcf_id);
					pw.print("\t");
					pw.println(filename);
					pw.flush();
					
					pstmt = this.conn.prepareStatement("SELECT ID FROM ROWCONTENT WHERE MD5SUM=?");
					pstmt2 = this.conn.prepareStatement("INSERT INTO ROWCONTENT(MD5SUM,CONTENT,CONTIG,START,STOP,ALLELE_REF) VALUES (?,?,?,?,?,?)",PreparedStatement.RETURN_GENERATED_KEYS);
					pstmt3 = this.conn.prepareStatement("INSERT INTO VCFROW(VCF_ID,ROW_ID) VALUES (?,?)");
					pstmt3.setLong(1, vcf_id);
					
					/* insert VCF header lines */
					for(final String line:headerLines) {
						final String md5=super.md5(line);
						long content_id = -1L;
						pstmt.setString(1, md5);
						row = pstmt.executeQuery();
						while(row.next()) {
							content_id = row.getLong(1);
						}
						row.close();				
						
						/* vcf content was not found, create it */
						if(content_id==-1L) {
							pstmt2.setString(1, md5);
							pstmt2.setString(2, line);
							pstmt2.setNull(3,Types.VARCHAR);
							pstmt2.setNull(4,Types.INTEGER);
							pstmt2.setNull(5,Types.INTEGER);
							pstmt2.setNull(6,Types.VARCHAR);
							if(pstmt2.executeUpdate()!=1) {
								return wrapException("Cannot insert ROWCONTENT ?");
							}
							content_id =getLastGeneratedId(pstmt2);
						}
						
						/* insert new VCF row */
						pstmt3.setLong(2, content_id);
						if(pstmt3.executeUpdate()!=1) {
							return wrapException("Cannot insert VCFROW ?");
						}
					}
					
					LOG.info("Inserted "+filename+" ID="+vcf_id);
					while(lineIter.hasNext() && !lineIter.peek().startsWith("#")) {
						final String line = lineIter.next();
						final String md5 = super.md5(line);
						
						long content_id = -1L;
						pstmt.setString(1, md5);
						row = pstmt.executeQuery();
						while(row.next()) {
							content_id = row.getLong(1);
						}
						row.close();
						/* vcf variants content was not found, create it */
						if(content_id==-1L) {
							/* decode to get chrom/start/end/ref */
							final VariantContext ctx = cah.codec.decode(line);
							pstmt2.setString(1, md5);
							pstmt2.setString(2, line);
							pstmt2.setString(3, ctx.getContig());
							pstmt2.setInt(4, ctx.getStart());
							pstmt2.setInt(5, ctx.getEnd());
							String refBase =ctx.getReference().getBaseString();
							/* sql table for Ref_allele is a varchar(50) */
							if(refBase.length()>50) {
								LOG.warn("Warning: TRUNCATING LARGE REF BASE TO 50 characters:"+refBase);
								refBase = refBase.substring(0,50);
							}
							pstmt2.setString(6,refBase );
							if(pstmt2.executeUpdate()!=1) {
								return wrapException("Cannot insert ROWCONTENT ?");
							}
							content_id =getLastGeneratedId(pstmt2);
						}
						
						/* insert new VCF row */
						pstmt3.setLong(2, content_id);
						if(pstmt3.executeUpdate()!=1) {
							return wrapException("Cannot insert VCFROW ?");
						}
					}
					pstmt2.close();
					pstmt3.close();
					pstmt.close();
					num_vcf_in_this_stream++;
					} /* end of while iter has next */
				CloserUtil.close(lineIter);
				lineIter=null;
				fileidx++;
			} while(fileidx < args.size());
			
			pw.flush();
			pw.close();
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(pw);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
			CloserUtil.close(pstmt2);
			CloserUtil.close(pstmt3);
			CloserUtil.close(lineIter);
		}
	}

			
	
	private Collection<Throwable> doCommandList(){
		Statement pstmt = null;
		ResultSet row = null;
		PrintWriter pw = null;
		if(!getInputFiles().isEmpty()) {
			return wrapException("Too many arguments");
		}
		try {
			/** user can add extra columns in the VCF, we collect the names of the columns */
			final List<String> cols =new ArrayList<>();
			pstmt = this.conn.createStatement();
			row = pstmt.executeQuery("SELECT  sys.SYSCOLUMNS.COLUMNNAME FROM  sys.SYSCOLUMNS, sys.SYSTABLES WHERE  sys.SYSTABLES.TABLEID= sys.SYSCOLUMNS.REFERENCEID AND  sys.SYSTABLES.TABLENAME=\'VCF\'");
			while(row.next()) {
				cols.add("VCF."+row.getString(1));
				}
			row.close();
			pstmt.close();
			final String sql=new StringBuilder("SELECT ").
					append(String.join(",", cols)).
					append(",COUNT(VCFROW.ID)  as \"COUNT_VARIANTS\" " 
						+ "FROM VCF,VCFROW,ROWCONTENT "
						+ "WHERE VCFROW.VCF_ID=VCF.ID AND VCFROW.ROW_ID = ROWCONTENT.ID AND ROWCONTENT.CONTIG IS NOT NULL "
						+ "GROUP BY "
						).
				append(String.join(",", cols)).
				toString()
				;
			
			LOG.info(sql);
			pstmt = this.conn.createStatement();
			row = pstmt.executeQuery(sql);
			final ResultSetMetaData meta = row.getMetaData();
			pw = super.openFileOrStdoutAsPrintWriter();
			for(int i=0;i< meta.getColumnCount();i++)
				{
				pw.print((i==0?"#":"\t"));
				pw.print(meta.getColumnLabel(i+1));
				}
			pw.println();
			while(row.next()) {
				for(int i=0;i< meta.getColumnCount();i++)
				{
				pw.print((i==0?"":"\t"));
				pw.print(row.getString(i+1));
				}
			pw.println();
			}
			
			pw.flush();
			return RETURN_OK;
		} catch (Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(pw);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
		}
	}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		try {
			openDerby();
			final String command  = String.valueOf(super.actionStr);
			
			if(command.equals("read")) {
				return doReadConcatenatedVcf();
			} else if(command.equals("list")) {
				return doCommandList();
			}else if(command.equals("dump")) {
				return doCommandDump();
			}
			else if(command.equals("dumpall")) {
				return doCommandDumpAll();
			}
			else {
				return wrapException("unknown command : "+command);
			}			
		} catch (Exception e) {
			return wrapException(e);
		} finally {
			closeDerby();
		}
		}
	
	
	public static void main(String[] args)
		{
		new VcfDerby01().instanceMainWithExit(args);
		}
	}
