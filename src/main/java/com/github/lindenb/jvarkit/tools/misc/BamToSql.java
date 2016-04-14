
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
* 2015 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMSequenceDictionaryHelper;

public class BamToSql
	extends AbstractBamToSql
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BamToSql.class);
	public BamToSql()
			{
			}

	private String quote(final String s ){
		if(s==null) return "NULL";
		final StringBuilder sb=new StringBuilder(s.length()+2);
		sb.append("'");
		sb.append(s);
		sb.append("'");
		return sb.toString();
	}
	
	@Override
	public Collection<Throwable> call() throws Exception 
		{				
		if(super.faidxFile==null) {
			return wrapException("ref sequence not defined (-"+OPTION_FAIDXFILE+")");
		}
		SAMRecordIterator iter=null;
		SamReader sfr=null;
		PrintWriter out =null;
		GenomicSequence genomicSequence=null;
		IndexedFastaSequenceFile indexedFastaSequenceFile=null;
		final List<String> args = new ArrayList<String>(IOUtils.unrollFiles(super.getInputFiles()));
		try
			{		
			
			out = super.openFileOrStdoutAsPrintWriter();
			indexedFastaSequenceFile=new IndexedFastaSequenceFile(super.faidxFile);

			
			out.println("CREATE TABLE IF NOT EXISTS SamFile");
			out.println("(");
			out.println("id INTEGER PRIMARY KEY,");
			out.println("filename TEXT");
			out.println(");");
			
			out.println("CREATE TABLE IF NOT EXISTS Dictionary");
			out.println("(");
			out.println("id INTEGER PRIMARY KEY,");
			out.println("name TEXT NOT NULL,");
			out.println("length INT NOT NULL,");
			out.println("tid INT NOT NULL,");
			out.println("samfile_id INT NOT NULL,");
			out.println("FOREIGN KEY(samfile_id) REFERENCES SamFile(id)");
			out.println(");");
			
			out.println("CREATE TABLE IF NOT EXISTS ReadGroup");
			out.println("(");
			out.println("id INTEGER PRIMARY KEY,");
			out.println("groupId TEXT NOT NULL,");
			out.println("sample TEXT NOT NULL,");
			out.println("samfile_id INT NOT NULL,");
			out.println("FOREIGN KEY(samfile_id) REFERENCES SamFile(id)");
			out.println(");");
			
			out.println("CREATE TABLE IF NOT EXISTS Read");
			out.println("(");
			out.println("id INTEGER PRIMARY KEY,");
			out.println("name TEXT NOT NULL,");
			out.println("flag INTEGER NOT NULL,");
			if(super.printflag){
				for(final SAMFlag flg: SAMFlag.values()) {
					out.println(flg.name()+" INTEGER NOT NULL,");
				}
			}
			out.println("rname TEXT,");
			out.println("pos INTEGER,");
			out.println("mapq INTEGER NOT NULL,");
			out.println("cigar TEXT,");
			out.println("rnext TEXT,");
			out.println("pnext INTEGER,");
			out.println("tlen INTEGER,");
			out.println("sequence TEXT NOT NULL,");
			out.println("qualities TEXT NOT NULL,");
			out.println("samfile_id INT NOT NULL,");
			out.println("group_id INT,");
			out.println("FOREIGN KEY(samfile_id) REFERENCES SamFile(id),");
			out.println("FOREIGN KEY(group_id) REFERENCES ReadGroup(id)");
			out.println(");");
			
			
			out.println("CREATE TABLE IF NOT EXISTS Cigar");
			out.println("(");
			out.println("id INTEGER PRIMARY KEY,");
			out.println("read_pos INT ,");
			out.println("read_base TEXT,");
			out.println("read_qual INT ,");
			out.println("ref_pos INT ,");
			out.println("ref_base TEXT,");
			out.println("operator TEXT NOT NULL,");
			out.println("read_id INT NOT NULL,");
			out.println("FOREIGN KEY(read_id) REFERENCES Read(id)");
			out.println(");");

			
			out.println("begin transaction;");


			
			int samIndex=0;
			do {
				final String inputName;
			
				if(samIndex==0 && args.isEmpty()) {
					sfr = openSamReader(null);
					inputName="<stdin>";
					}
				else
					{
					inputName= args.get(samIndex);
					sfr = openSamReader(inputName);
					}
				final SAMFileHeader header1=sfr.getFileHeader();
				if(header1==null)
					{
					return wrapException("File header missing");
					}
				final SAMSequenceDictionary dict=header1.getSequenceDictionary();
				if(dict==null)  {
					return wrapException("No Dictionary in input");
					}
				final SAMSequenceDictionaryHelper dix = new SAMSequenceDictionaryHelper(dict);
				
				final Interval userInterval;
				iter= null;
				if(super.regionStr==null || super.regionStr.isEmpty()) {
					LOG.warn("You're currently scanning the whole BAM ???!!!");
					iter = sfr.iterator();
					userInterval = null;
					}
				else
					{
					final Optional<Interval> interval = dix.parseInterval(super.regionStr);
					if(!interval.isPresent()) {
						return wrapException("cannot parse interval "+super.regionStr);
					}
					userInterval = interval.get();
					iter = sfr.query(
							interval.get().getContig(),
							interval.get().getStart(),
							interval.get().getEnd(),
							false
							);
					}
				out.println(String.join(" ",
						"insert into SamFile(filename) values(",
						quote(inputName),
						");"
						));
				
				for(int i=0;i< dict.size();++i) {
					final SAMSequenceRecord ssr = dict.getSequence(i);
					out.println(
						"insert into Dictionary(name,length,tid,samfile_id) select "+
						quote(inputName) + ","+
						ssr.getSequenceLength()+","+
						i+",max(id) from SamFile;"
						);
					}
				for(final SAMReadGroupRecord g:header1.getReadGroups()){
					out.println(
							"insert into ReadGroup(groupId,sample,samfile_id) select "+
							quote(g.getId()) + ","+
							quote(g.getSample())+","+
							"max(id) from SamFile;"
							);
					}
				
				
				
				final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1);

				while(iter.hasNext())
					{
					final SAMRecord rec= progress.watch(iter.next());
					final StringBuilder sql = new StringBuilder();
					sql.append("insert into Read("
							+ "name,flag,");
					if(super.printflag){
						for(final SAMFlag flg: SAMFlag.values()) {
							sql.append(flg.name()).append(",");
						}
					}
					sql.append("rname,pos,mapq,cigar,rnext,pnext,tlen,sequence,qualities,group_id,samfile_id) select ");

					
					
				    sql.append(quote(rec.getReadName())).append(",");
				    sql.append(rec.getFlags()).append(",");
				    
				    if(super.printflag){
						for(final SAMFlag flg: SAMFlag.values()) {
							sql.append(flg.isSet(rec.getFlags())?1:0);
							sql.append(",");
						}
					}
				    
				    if(rec.getReferenceName()==null || rec.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
				    	sql.append("NULL,NULL");
				    } else
				    {
				      sql.append(quote(rec.getReferenceName()));
			    	  sql.append(",");
			    	  sql.append(rec.getAlignmentStart());
				    }
				    sql.append(",");
					sql.append(rec.getMappingQuality());
					sql.append(",");
					
					//cigar
					if(rec.getCigarString()==null || rec.getCigarString().equals(SAMRecord.NO_ALIGNMENT_CIGAR)) {
				    	sql.append("NULL");
				    } else
				    {
				      sql.append(quote(rec.getCigarString()));
				    }
					sql.append(",");
					
					//rnext
					 if(rec.getMateReferenceName()==null || rec.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
					    	sql.append("NULL,NULL");
					    } else
					    {
					      sql.append(quote(rec.getMateReferenceName()));
				    	  sql.append(",");
				    	  sql.append(rec.getMateAlignmentStart());
					    }
					sql.append(",");
				    
					//tlen
					sql.append(rec.getInferredInsertSize());
					sql.append(",");
					//sequence
					sql.append(quote(rec.getReadString()));
					sql.append(",");
					//qualities
					sql.append(quote(rec.getBaseQualityString()));
					sql.append(",");
					
					if(rec.getReadGroup()==null)
						{
						sql.append("NULL");
						}
					else
						{
						sql.append("G.id");
						}
					sql.append(",F.id FROM SamFile as F");
					if(rec.getReadGroup()!=null)
						{
						sql.append(" , ReadGroup as G where G.groupId=").
							append(quote(rec.getReadGroup().getId())).
							append(" and F.id = G.samfile_id ");
						
	 					}
					sql.append("  ORDER BY F.id DESC LIMIT 1;");
					out.println(sql.toString());
					
					if(super.printcigar && !rec.getReadUnmappedFlag() && rec.getCigar()!=null) {
						if(genomicSequence==null || !genomicSequence.getChrom().equals(rec.getReferenceName())) {
							genomicSequence = new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
						}
						int ref = rec.getUnclippedStart();
						final byte bases[]=rec.getReadBases();
						final byte quals[]=rec.getBaseQualities();
						int read  = 0;
						for(final CigarElement ce:rec.getCigar()) {
							final CigarOperator op=ce.getOperator();
							if(op.equals(CigarOperator.P)) continue;
							for(int i=0;i< ce.getLength();++i) {
								sql.setLength(0);
								boolean in_user_interval=true;
								sql.append("insert into Cigar(operator,read_pos,read_base,read_qual,ref_pos,ref_base,read_id) ");
								sql.append("select '");
								sql.append(op.name());
								sql.append("',");
								
								if(userInterval!=null && 
									!(rec.getReferenceName().equals(userInterval.getContig()) &&
											ref>=userInterval.getStart() && ref<=userInterval.getEnd()))
									{
									in_user_interval = false;
									}
								
								switch(op){
									case I: {
										sql.append(read);	
										sql.append(",");
										sql.append("'"+(char)bases[read]+"',");
										sql.append(""+quals[read]+"");
										sql.append(",");
										sql.append("NULL,NULL");	
										read++;
										break;
										}
									case D:case N:case H://yes H (hard clip)
										{
										sql.append("NULL,NULL,NULL,");
										sql.append(ref);	
										sql.append(",'");
										sql.append((ref<1 || ref-1>=genomicSequence.length())?'*':genomicSequence.charAt(ref-1));
										sql.append("'");
										ref++;
										break;
										}
									case M:case X:case EQ:case S: //yes S, soft clip
										{
										sql.append(read);	
										sql.append(",");
										sql.append("'"+(char)bases[read]+"',");
										sql.append(""+quals[read]+"");
										sql.append(",");
										sql.append(ref);	
										sql.append(",'");
										sql.append((ref<1 || ref-1>=genomicSequence.length())?'*':genomicSequence.charAt(ref-1));
										sql.append("'");
										ref++;
										read++;
										break;
										}
									default: throw new IllegalStateException();
									}
								
								
								
								sql.append(", id from Read ORDER BY id DESC LIMIT 1;");
								if(in_user_interval) out.println(sql.toString());
			
							}
						}
					}
					
					}
				iter.close();iter=null;
				sfr.close();sfr=null;
				progress.finish();
				samIndex++;
			} while(samIndex< args.size());
			
			
		
			
			
			
			out.println("COMMIT;");
			out.flush();
			out.close();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(out);
			CloserUtil.close(indexedFastaSequenceFile);			
			}
		}
				
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BamToSql().instanceMainWithExit(args);
		}
	}
