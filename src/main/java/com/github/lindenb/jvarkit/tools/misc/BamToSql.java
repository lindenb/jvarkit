
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


History:
* 2015 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
/**
 BEGIN_DOC
 
## Motivation

Inserting a BAM in a SQL is not a good idea of course !

But it might be interesting to get some informations about the bases in a segment of bam.

## Schema
The schema can change if some options (-c , -f) are used.
At the time of writing the schema is :

```sql
CREATE TABLE IF NOT EXISTS SamFile
(
id INTEGER PRIMARY KEY,
filename TEXT
);

CREATE TABLE IF NOT EXISTS Dictionary
(
id INTEGER PRIMARY KEY,
name TEXT NOT NULL,
length INT NOT NULL,
tid INT NOT NULL,
samfile_id INT NOT NULL,
FOREIGN KEY(samfile_id) REFERENCES SamFile(id)
);

CREATE TABLE IF NOT EXISTS ReadGroup
(
id INTEGER PRIMARY KEY,
groupId TEXT NOT NULL,
sample TEXT NOT NULL,
samfile_id INT NOT NULL,
FOREIGN KEY(samfile_id) REFERENCES SamFile(id)
);

CREATE TABLE IF NOT EXISTS Read
(
id INTEGER PRIMARY KEY,
name TEXT NOT NULL,
flag INTEGER NOT NULL,
rname TEXT,
pos INTEGER,
mapq INTEGER NOT NULL,
cigar TEXT,
rnext TEXT,
pnext INTEGER,
tlen INTEGER,
sequence TEXT NOT NULL,
qualities TEXT NOT NULL,
samfile_id INT NOT NULL,
group_id INT,
FOREIGN KEY(samfile_id) REFERENCES SamFile(id),
FOREIGN KEY(group_id) REFERENCES ReadGroup(id)
);

CREATE TABLE IF NOT EXISTS Cigar
(
id INTEGER PRIMARY KEY,
read_pos INT ,
read_base TEXT,
read_qual INT ,
ref_pos INT ,
ref_base TEXT,
operator TEXT NOT NULL,
read_id INT NOT NULL,
FOREIGN KEY(read_id) REFERENCES Read(id)
);
```

## Example

Build a sqlite3 database for a set of BAM files in the region "rotavirus:1-10""

```
$java -jar dist/bam2sql.jar -r 'rotavirus:1-10' -R  ref.fa -c S*.bam |\
sqlite3 database.sqlite
```

Select data from sqlite database where the genomic position is "rotavirus:5" 

```sql
select  SamFile.filename,
		ReadGroup.sample,
		Read.flag,
		Read.rname,
		Cigar.operator,
		Cigar.read_pos,
		Cigar.read_base,
		Cigar.read_qual,
		Cigar.ref_pos,
		Cigar.ref_base
from
		SamFile,Read,Cigar,ReadGroup
where
		SamFile.id = Read.samfile_id AND
		ReadGroup.id = Read.group_id AND 
		Cigar.read_id = Read.id and
		Read.rname = "rotavirus" and 
		Cigar.ref_pos= 5
		;
```

query:

```
$ sqlite3 -header -separator '   ' database.sqlite &lt; query.sql  | column -t 
```

output:

```
filename  sample  flag  rname      operator  read_pos  read_base  read_qual  ref_pos  ref_base
S1.bam    S1      99    rotavirus  M         4         T          10         5        T
S1.bam    S1      163   rotavirus  M         4         T          10         5        T
S1.bam    S1      163   rotavirus  M         4         T          10         5        T
S1.bam    S1      99    rotavirus  M         4         T          10         5        T
S1.bam    S1      163   rotavirus  M         4         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      163   rotavirus  M         3         T          10         5        T
S1.bam    S1      163   rotavirus  M         3         T          10         5        T
S1.bam    S1      163   rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      163   rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         2         T          10         5        T
S1.bam    S1      99    rotavirus  M         2         T          10         5        T
S1.bam    S1      163   rotavirus  M         2         T          10         5        T
S1.bam    S1      163   rotavirus  M         2         T          10         5        T
S1.bam    S1      99    rotavirus  M         1         T          10         5        T
S1.bam    S1      99    rotavirus  M         1         T          10         5        T
S1.bam    S1      163   rotavirus  M         1         T          10         5        T
S1.bam    S1      163   rotavirus  M         1         T          10         5        T
S1.bam    S1      163   rotavirus  M         4         T          10         5        T
S1.bam    S1      163   rotavirus  M         0         T          10         5        T
S1.bam    S1      163   rotavirus  M         0         T          10         5        T
S1.bam    S1      163   rotavirus  M         0         T          10         5        T
S1.bam    S1      163   rotavirus  M         0         T          10         5        T
S1.bam    S1      99    rotavirus  S         3         T          10         5        T
S1.bam    S1      163   rotavirus  S         0         T          10         5        T
S1.bam    S1      99    rotavirus  S         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         4         T          10         5        T
S2.bam    S2      163   rotavirus  M         4         T          10         5        T
S2.bam    S2      99    rotavirus  M         4         T          10         5        T
S2.bam    S2      99    rotavirus  M         4         T          10         5        T
S2.bam    S2      163   rotavirus  M         4         T          10         5        T
S2.bam    S2      163   rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         2         A          10         5        T
S2.bam    S2      163   rotavirus  M         2         T          10         5        T
S2.bam    S2      99    rotavirus  M         1         T          10         5        T
S2.bam    S2      99    rotavirus  M         1         T          10         5        T
S2.bam    S2      163   rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         1         T          10         5        T
S2.bam    S2      99    rotavirus  M         1         T          10         5        T
S2.bam    S2      99    rotavirus  M         0         T          10         5        T
S2.bam    S2      99    rotavirus  M         0         T          10         5        T
S2.bam    S2      163   rotavirus  S         4         T          10         5        T
S2.bam    S2      99    rotavirus  S         2         A          10         5        T
S3.bam    S3      99    rotavirus  M         4         A          10         5        T
S3.bam    S3      163   rotavirus  M         4         T          10         5        T
S3.bam    S3      99    rotavirus  M         4         T          10         5        T
S3.bam    S3      99    rotavirus  M         3         T          10         5        T
S3.bam    S3      99    rotavirus  M         3         T          10         5        T
S3.bam    S3      99    rotavirus  M         3         T          10         5        T
S3.bam    S3      163   rotavirus  M         3         T          10         5        T
S3.bam    S3      163   rotavirus  M         3         T          10         5        T
S3.bam    S3      163   rotavirus  M         2         T          10         5        T
S3.bam    S3      163   rotavirus  M         2         T          10         5        T
S3.bam    S3      99    rotavirus  M         2         T          10         5        T
S3.bam    S3      163   rotavirus  M         2         T          10         5        T
S3.bam    S3      99    rotavirus  M         1         T          10         5        T
S3.bam    S3      163   rotavirus  M         1         A          10         5        T
S3.bam    S3      99    rotavirus  M         1         A          10         5        T
S3.bam    S3      99    rotavirus  M         1         A          10         5        T
S3.bam    S3      99    rotavirus  M         1         T          10         5        T
S3.bam    S3      163   rotavirus  M         1         T          10         5        T
S3.bam    S3      99    rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         A          10         5        T
S3.bam    S3      99    rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         T          10         5        T
S3.bam    S3      99    rotavirus  S         2         A          10         5        T
S4.bam    S4      163   rotavirus  M         4         T          10         5        T
S4.bam    S4      163   rotavirus  M         4         T          10         5        T
S4.bam    S4      99    rotavirus  M         4         T          10         5        T
S4.bam    S4      163   rotavirus  M         4         T          10         5        T
S4.bam    S4      163   rotavirus  M         3         T          10         5        T
S4.bam    S4      163   rotavirus  M         3         T          10         5        T
S4.bam    S4      99    rotavirus  M         3         T          10         5        T
S4.bam    S4      163   rotavirus  M         2         T          10         5        T
S4.bam    S4      99    rotavirus  M         1         T          10         5        T
S4.bam    S4      99    rotavirus  M         0         T          10         5        T
S4.bam    S4      99    rotavirus  M         4         T          10         5        T
S4.bam    S4      163   rotavirus  M         0         T          10         5        T
S4.bam    S4      163   rotavirus  M         0         T          10         5        T
```

 
 END_DOC
 */
@Program(name="bam2sql",
	description="Convert a SAM/BAM to sqlite statements",
	keywords={"bam","sam","sql","sqlite"}
		)
public class BamToSql
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BamToSql.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-r","--region"},description=IntervalParser.OPT_DESC)
	private String regionStr = "";

	@Parameter(names={"-c","--cigar"},description="print cigar data")
	private boolean printcigar = false;

	@Parameter(names={"-f","--flag"},description="expands details about sam flag")
	private boolean printflag = false;

	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidxFile=null;
	
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
	public int doWork(List<String> args) {				
		if(this.faidxFile==null) {
			LOG.error("ref sequence faidx not defined");
			return -1;
			}
		SAMRecordIterator iter=null;
		SamReader sfr=null;
		PrintWriter out =null;
		GenomicSequence genomicSequence=null;
		IndexedFastaSequenceFile indexedFastaSequenceFile=null;
		args = new ArrayList<String>(IOUtils.unrollFiles(args));
		try
			{		
			
			out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			indexedFastaSequenceFile=new IndexedFastaSequenceFile(this.faidxFile);

			
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
			if(this.printflag){
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
					throw new JvarkitException.FileFormatError("File header missing");
					}
				final SAMSequenceDictionary dict=header1.getSequenceDictionary();
				if(dict==null)  {
					throw new JvarkitException.DictionaryMissing("No Dictionary in input");
					}
				final IntervalParser intervalParser = new IntervalParser(dict);
				
				final Interval userInterval;
				iter= null;
				if(this.regionStr==null || this.regionStr.isEmpty()) {
					LOG.warn("You're currently scanning the whole BAM ???!!!");
					iter = sfr.iterator();
					userInterval = null;
					}
				else
					{
					userInterval = intervalParser.parse(this.regionStr);
					if(userInterval==null) {
						throw new JvarkitException.UserError("cannot parse interval "+this.regionStr);
					}
					iter = sfr.query(
							userInterval.getContig(),
							userInterval.getStart(),
							userInterval.getEnd(),
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
					if(this.printflag){
						for(final SAMFlag flg: SAMFlag.values()) {
							sql.append(flg.name()).append(",");
						}
					}
					sql.append("rname,pos,mapq,cigar,rnext,pnext,tlen,sequence,qualities,group_id,samfile_id) select ");

					
					
				    sql.append(quote(rec.getReadName())).append(",");
				    sql.append(rec.getFlags()).append(",");
				    
				    if(this.printflag){
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
					
					if(this.printcigar && !rec.getReadUnmappedFlag() && rec.getCigar()!=null) {
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
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
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
