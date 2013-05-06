package fr.inserm.umr1087.jvarkit.tools.cmpbams;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import com.sleepycat.bind.tuple.StringBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

public class CompareBams {
	private static final Logger LOG=Logger.getLogger(CompareBams.class.getName());
	private static final String DATABASENAME="read2pos";
	private File dbHome;
	private Environment environment=null;
	private Database database=null;
	private Transaction txn;
	private int countBams;
	private boolean useSamFlag=false;
	
	
	private static class Match
		{
		byte tid;
		int pos;
		int flag=0;

		@Override
		public int hashCode()
			{
			int result = 1;
			result = 31 * result + pos;
			result = 31 * result + tid;
			//NO result = 31 * result + flag;
			return result;
			}
		@Override
		public boolean equals(Object obj)
			{
			if (this == obj) { return true; }
			if (obj == null) { return false; }
			Match other = (Match) obj;
			if (tid != other.tid) { return false; }
			if(tid==-1) return true;
			if (pos != other.pos) { return false; }
			//NO if (flag != other.flag) { return false; }
			return true;
			}
		}
	
	/** fileid to matches */
	private List<Set<Match>> decode(final DatabaseEntry data)
		throws IOException
		{
		TupleInput in=new TupleInput(data.getData());
		List<Set<Match>> L=new ArrayList<Set<Match>>(this.countBams);
		for(int i=0;i< this.countBams;++i)
			{
			byte nmatches=in.readByte();
			Set<Match> set=new HashSet<Match>(nmatches);
			L.add(set);
			for(int j=0;j< nmatches;++j)
				{
				Match m=new Match();
				m.tid=in.readByte();
				m.pos=in.readInt();
				m.flag=this.useSamFlag?in.readInt():0;
				set.add(m);
				}
			}
		
		return L;
		}
	

	
	private void encode(DatabaseEntry data,final List<Set<Match>> L)
		{
		TupleOutput out=new TupleOutput();
		for(Set<Match> set:L)
			{
			out.writeByte((byte)Math.min(Byte.MAX_VALUE, set.size()));
			int count=0;
			for(Match m:set)
				{
				if(++count>=Byte.MAX_VALUE) break;
				out.writeByte(m.tid);
				out.writeInt(m.pos);
				if(this.useSamFlag) out.writeInt(m.flag);
				}
			}
		data.setData(out.getBufferBytes(),out.getBufferOffset(),out.getBufferLength());
		}
	
	private void print(final Set<Match> set,final SAMSequenceDictionary dict)
		{
		boolean first=true;
		for(Match m:set)
			{
			if(!first)System.out.print(',');
			first=false;
			if(m.tid<0){ System.out.print("unmapped"); continue;}
			SAMSequenceRecord ssr=(dict==null?null:dict.getSequence(m.tid));
			String seqName=(ssr==null?null:ssr.getSequenceName());
			if(seqName==null) seqName="tid"+m.tid;
			System.out.print(String.valueOf(seqName+":"+(m.pos)));
			if(this.useSamFlag) System.out.print("="+m.flag);
			}
		if(first) System.out.print("(empty)");
		}
	
	
	private void run(String[] args)
		throws Exception
		{
		int optind=0;
		while(optind< args.length)
			{
			if(args[optind].equals("-h") ||
			   args[optind].equals("-help") ||
			   args[optind].equals("--help"))
				{
				System.err.println("Pierre Lindenbaum PhD. 2013");
				System.err.println("Options:");
				System.err.println(" -h help; This screen.");
				System.err.println(" -d <berkeleydb-dir>.");
				System.err.println(" -F use samFlag");
				return;
				}
			else if(args[optind].equals("-d") && optind+1< args.length)
				{
				this.dbHome=new File(args[++optind]);
				}
			else if(args[optind].equals("-F") )
				{
				this.useSamFlag=true;
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unknown option "+args[optind]);
				return;
				}
			else 
				{
				break;
				}
			++optind;
			}
		if(this.dbHome==null)
			{
			System.err.println("db-home undefined");
			return;
			}
		this.countBams=args.length-optind;
		if(this.countBams<2)
			{
			System.err.println("Need more bams please");
			return;
			}
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();

		EnvironmentConfig envConfig= new EnvironmentConfig();
		envConfig.setAllowCreate(true);
		envConfig.setReadOnly(false);
		envConfig.setConfigParam(EnvironmentConfig.LOG_FILE_MAX,"250000000");
		envConfig.setTransactional(true);
		this.environment= new Environment(this.dbHome, envConfig);
		
		this.txn=this.environment.beginTransaction(null, null);
		DatabaseConfig cfg= new DatabaseConfig();
		cfg.setAllowCreate(true);
		cfg.setReadOnly(false);
		cfg.setTransactional(true);
		cfg.setKeyPrefixing(true);
		this.database= this.environment.openDatabase(txn,DATABASENAME,cfg);
		
		List<SAMSequenceDictionary> sequenceDictionaries=new ArrayList<SAMSequenceDictionary>(this.countBams);
		
		for(int currentSamFileIndex=0;
				currentSamFileIndex<this.countBams;
				currentSamFileIndex++ )
			{
			long nReads=0L;
			File samFile=new File(args[optind+currentSamFileIndex]);
			LOG.info("Opening "+samFile);
			SAMFileReader samFileReader=new SAMFileReader(samFile);
			samFileReader.setValidationStringency(ValidationStringency.SILENT);
			SAMSequenceDictionary dict=samFileReader.getFileHeader().getSequenceDictionary();
			if(dict.getSequences().size()>=Byte.MAX_VALUE)
				{
				System.err.println("Too many Ref Sequences ("+ dict.getSequences().size() +") . Limited to "+Byte.MAX_VALUE);
				return;
				}
			sequenceDictionaries.add(dict);
			
			for(Iterator<SAMRecord> iter=samFileReader.iterator();
					iter.hasNext(); )
				{
				if(nReads++%10000000==0) LOG.info("in "+samFile+" count:"+nReads);
				SAMRecord rec=iter.next();
				StringBinding.stringToEntry(rec.getReadName(), key);
				
				List<Set<Match>> matches=null;
				if(this.database.get(this.txn, key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
					{
					matches=decode(data); 
					}
				else
					{
					matches=new ArrayList<Set<Match>>( this.countBams);
					for(int i=0;i< this.countBams;++i) matches.add(new HashSet<Match>());
					}
				
				Match match=new Match();
				match.tid=(byte)rec.getReferenceIndex().intValue();
				match.pos=rec.getAlignmentStart();
				match.flag=rec.getFlags();
				matches.get(currentSamFileIndex).add(match);
				encode(data,matches);
				if(this.database.put(this.txn, key, data)!=OperationStatus.SUCCESS)
					{
					System.err.println("BDB error.");
					System.exit(-1);
					}
				}
			samFileReader.close();
			LOG.info("Close "+samFile);
			}
		LOG.info("Writing results....");
		//compute the differences for each read
		System.out.print("#READ-Name\t");
		for(int x=0;x<this.countBams;++x)
			{
			for(int y=x+1;y<this.countBams;++y)
				{
				if(!(x==0 && y==1)) System.out.print("|");
				System.out.print(args[optind+x]);
				System.out.print(" ");
				System.out.print(args[optind+y]);
				}
			}
		for(int x=0;x<this.countBams;++x)
			{
			System.out.print("\t"+args[optind+x]);
			}
		System.out.println();
		
		key=new DatabaseEntry();
		Cursor c=this.database.openCursor(txn,null);;
		while(c.getNext(key, data,LockMode.DEFAULT)==OperationStatus.SUCCESS)
			{
			System.out.print(StringBinding.entryToString(key));
			System.out.print("\t");
			
			List<Set<Match>> matches=decode(data);

			for(int x=0;x<this.countBams;++x)
				{
				Set<Match> first=matches.get(x);
				for(int y=x+1;y<this.countBams;++y)
					{
					if(!(x==0 && y==1)) System.out.print("|");
					Set<Match> second=matches.get(y);
					if(first.equals(second))
						{
						System.out.print("EQ");
						}
					else
						{
						System.out.print("NE");
						}
					}
				}

			for(int x=0;x<this.countBams;++x)
				{
				System.out.print("\t");
				print(matches.get(x),sequenceDictionaries.get(x));
				}
			
			System.out.println();
			}
		c.close();
		this.database.close();
		this.environment.removeDatabase(txn, DATABASENAME);
		this.txn.commit();
		this.environment.cleanLog();
		this.environment.close();
		}
		
	public static void main(String[] args) throws Exception
		{
		new CompareBams().run(args);
		}
}
