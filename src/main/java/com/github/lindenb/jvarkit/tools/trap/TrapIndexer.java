package com.github.lindenb.jvarkit.tools.trap;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Example:

```
java -jar dist/trapindexer.jar  -o chr22.dat  chr22.TraPv2.txt.gz
```

## See also

* VcfTrap

END_DOC
 */
@Program(name="trapindexer",
	description="Convert text data to binary format for the trap DATABASE database http://trap-score.org/. Those data can be used by the tool `vcftrap`.",
	keywords= {"trap"}
)
public class TrapIndexer extends Launcher{
private static final Logger LOG = Logger.build(TrapIndexer.class).make();
@Parameter(names= {"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
private File outfilename=null;

static final int SCORE_STRLEN=5;/* 0.123 */
private static final int SCORE_SIZEOF=SCORE_STRLEN - 2 /* remove leading '0.' byte[0]==1 is score=1f*/;
private static final int ENSG_STRLEN=15;
static final int RECORD_SIZOF= Integer.BYTES /* pos */ + 2/*ref and alt CANNOT be in one byte. Saw some data REF=M, ALT=N */+ Integer.BYTES /* ensgene-id */+ SCORE_SIZEOF /*score */;
static final byte MAGIC[]= "TRAP.1.1".getBytes();

public static TrapRecord decode(final String contig,byte array[]) {
	if(array.length!=RECORD_SIZOF) throw new IllegalStateException("byte.length "+array.length+"!="+RECORD_SIZOF);
	try {
		final DataInputStream dis = new DataInputStream(new ByteArrayInputStream(array));
		final int pos = dis.readInt();
		if(pos<0) throw new IOException("pos<0 : "+pos);
		final byte ref = dis.readByte();
		final byte alt = dis.readByte();
		int ensgId = dis.readInt();
		final String ensg = String.format("ENSG%0"+(ENSG_STRLEN-4)+"d",ensgId);
		byte score_bytes[]=new byte[SCORE_SIZEOF];
		dis.readFully(score_bytes);
		final float score ;
		if( score_bytes[0] == (byte)1)
			{	
			score=1.0f;
			}
		else
			{
			score = Float.parseFloat("0."+new String(score_bytes));
			}
		return new TrapRecord() {
			@Override
			public int getStart() { return pos; }
			@Override
			public int getEnd() { return pos; }
			@Override
			public String getContig() { return contig;}
			@Override
			public String getChr() { return getContig(); }
			@Override
			public float getScore() { return score; }
			@Override
			public char getRef() {return (char)ref;}
			@Override
			public String getGene() {return ensg;}
			@Override
			public char getAlt() { return  (char)alt; }
			@Override
			public String toString() {
				return contig+":"+pos+":"+(char)ref+"/"+(char)alt+" "+ensg+" "+score;
				}
		};
	} catch (IOException e) {
		throw new RuntimeIOException(e);
		}
	}


@Override
public int doWork(final List<String> args) {
	
	final Pattern TAB=Pattern.compile("[\t]");
	OutputStream fos = null;
	BufferedReader r = null;
	try
		{
		final File inputFile = new File(oneAndOnlyOneFile(args));
		IOUtil.assertFileIsReadable(inputFile);
	
		if(!inputFile.getName().startsWith("chr"))
			{
			LOG.error("File doesn't starts with chr :"+inputFile);
			return -1;
			}
		int dot = inputFile.getName().indexOf(".");
		if(dot==-1) {
			LOG.error("Cannot find dot in "+inputFile.getName());
			return -1;
			}

	final String contig = inputFile.getName().substring(3,dot);
	
	
	if(args.isEmpty()) {
		LOG.error("No input file defined");
		return -1;
		}
	if(this.outfilename==null)
		{
		fos = stdout();
		}
	else
		{
		if(!this.outfilename.getName().endsWith(".dat"))
			{
			LOG.error("output filename doesn't end with *.dat :"+outfilename);
			return -1;
			}
		if(!this.outfilename.getName().contains(contig+"."))
			{
			LOG.error("output should contain with \""+contig+".\"  :"+outfilename);
			return -1;
			}
		LOG.info("Opening "+this.outfilename);
		// outputstream;
		fos = new FileOutputStream(this.outfilename);
		}
	
	fos.write(MAGIC);
	
	final byte score_str[]=new byte[SCORE_SIZEOF];
	long nRecords=0L;
	int prev_pos=0;
	LOG.info("Reading "+inputFile+" contig:("+contig+")");
	r=IOUtils.openFileForBufferedReading(inputFile);
	String line;
	while((line=r.readLine())!=null)
		{
		if(StringUtil.isBlank(line)) {
			LOG.error("Blank line in "+inputFile+" ??");
			return -1;
			}
		nRecords++;
		if(nRecords%1000000==0)
			{
			LOG.info("indexing "+nRecords+" last \""+line+"\"");
			}
		final String tokens[] = TAB.split(line);
		
		if(tokens.length<5) throw new JvarkitException.TokenErrors(tokens);
		if(StringUtil.isBlank(tokens[0])) throw new JvarkitException.FileFormatError("empty pos in "+line);
		if(StringUtil.isBlank(tokens[1])) throw new JvarkitException.FileFormatError("empty ref in "+line);
		if(StringUtil.isBlank(tokens[2])) throw new JvarkitException.FileFormatError("empty alt in "+line);
		if(tokens[1]==tokens[2]) throw new JvarkitException.FileFormatError("ref==alt in "+line);
		if(tokens[1].length()!=1) throw new JvarkitException.FileFormatError("bad ref "+line);
		if(tokens[2].length()!=1) throw new JvarkitException.FileFormatError("bad alt "+line);
		if(tokens[3].length()>ENSG_STRLEN) throw new JvarkitException.FileFormatError("strlen(ensGene)> "+ENSG_STRLEN+" in "+line);
		if(!tokens[3].startsWith("ENSG")) throw new JvarkitException.FileFormatError("ensGene doesn't start with ENSG  in "+line);
		if(tokens[4].length()>SCORE_STRLEN) throw new JvarkitException.FileFormatError("strlen(score)> "+SCORE_STRLEN+" in "+line);
		final float f = Float.parseFloat(tokens[4]);
		if(f<0)  throw new JvarkitException.FileFormatError("bad score "+tokens[4]+"->"+f+"<0");
		if(f>1.0f)  throw new JvarkitException.FileFormatError("bad score "+tokens[4]+"->"+f+">1.0");
		
		final int pos = Integer.parseInt(tokens[0]);
		
		if(prev_pos>pos )
			{
			LOG.error("input is not sorted "+prev_pos +" before "+pos);
			return -1;
			}
		prev_pos=pos;
		
		
		final ByteArrayOutputStream baos = new ByteArrayOutputStream(RECORD_SIZOF);
		final DataOutputStream daos = new DataOutputStream(baos);
		daos.writeInt(pos);
		daos.writeByte(tokens[1].charAt(0));
		daos.writeByte(tokens[2].charAt(0));
		daos.writeInt(Integer.parseInt(tokens[3].substring(4)));//after 'ENSG'
		
		Arrays.fill(score_str,(byte)'0');
		if(tokens[4]. equals("1") || tokens[4].equals("1.0"))
			{
			score_str[0]=(byte)1;//not a char !! just value 1
			}
		else if(tokens[4].equals("0"))
			{
			//0 already set
			}
		else if(tokens[4].startsWith("0.")) {
			final byte b_array[] = tokens[4].substring(2).getBytes();
			if(b_array.length>SCORE_SIZEOF) throw new JvarkitException.FileFormatError("a(score)> "+SCORE_STRLEN+" in "+line);
			System.arraycopy(b_array,0,score_str,0,b_array.length);
			}
		else
			{
			LOG.error("Bad Score: "+tokens[4]);
			return -1;
			}
		
		daos.write(score_str, 0, score_str.length);
		
		
		
		final byte record_bytes[] = baos.toByteArray();
		if(record_bytes.length!=RECORD_SIZOF)
			{
			LOG.error("Bad record size: "+record_bytes.length);
			return -1;
			}
		
		fos.write(record_bytes);
		}

	r.close();

	fos.flush();
	fos.close();
	LOG.info("Done "+nRecords+" lines. File size:"+(RECORD_SIZOF*nRecords));
	return 0;
	}
catch(final Exception err)
	{
	LOG.error(err);
	return -1;
	}
finally
	{
	CloserUtil.close(fos);
	CloserUtil.close(r);
	}
}
	
public static void main(String[] args) {
	new TrapIndexer().instanceMainWithExit(args);
	}
}
