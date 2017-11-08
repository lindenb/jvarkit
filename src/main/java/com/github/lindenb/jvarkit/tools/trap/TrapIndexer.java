package com.github.lindenb.jvarkit.tools.trap;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

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

static final int SCORE_STRLEN=5;
static final int ENSG_STRLEN=15;
static final int RECORD_SIZOF= Integer.BYTES /* pos */ + 1/*ref and alt in one byte*/+ Integer.BYTES /* ensgene-id */+ SCORE_STRLEN /*score */;
static final byte MAGIC[]= "TRAP.1.0".getBytes();

public static TrapRecord decode(final String contig,byte array[]) {
	if(array.length!=RECORD_SIZOF) throw new IllegalStateException("byte.length "+array.length+"!="+RECORD_SIZOF);
	try {
		final DataInputStream dis = new DataInputStream(new ByteArrayInputStream(array));
		final int pos = dis.readInt();
		if(pos<0) throw new IOException("pos<0 : "+pos);
		final byte opcode = dis.readByte();
		final byte bases[] = opcode2bases(opcode);
		int ensgId = dis.readInt();
		final String ensg = String.format("ENSG%0"+(ENSG_STRLEN-4)+"d",ensgId);
		byte score_bytes[]=new byte[SCORE_STRLEN];
		dis.readFully(score_bytes);
		final float score = Float.parseFloat(new String(score_bytes));
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
			public char getRef() {return (char)bases[0];}
			@Override
			public String getGene() {return ensg;}
			@Override
			public char getAlt() { return  (char)bases[1]; }
			@Override
			public String toString() {
				return contig+":"+pos+":"+(char)bases[0]+"/"+(char)bases[1]+" "+ensg+" "+score;
				}
		};
	} catch (IOException e) {
		throw new RuntimeIOException(e);
		}
	}

private static byte[] opcode2bases(final byte opcode)
	{
	switch(opcode) {
		case 1 : return new byte[]{(byte)'A',(byte)'A'};
		case 2 : return new byte[]{(byte)'A',(byte)'C'};
		case 3 : return new byte[]{(byte)'A',(byte)'G'};
		case 4 : return new byte[]{(byte)'A',(byte)'T'};
		case 5 : return new byte[]{(byte)'A',(byte)'N'};
		case 6 : return new byte[]{(byte)'C',(byte)'A'};
		case 7 : return new byte[]{(byte)'C',(byte)'C'};
		case 8 : return new byte[]{(byte)'C',(byte)'G'};
		case 9 : return new byte[]{(byte)'C',(byte)'T'};
		case 10 : return new byte[]{(byte)'C',(byte)'N'};
		case 11 : return new byte[]{(byte)'G',(byte)'A'};
		case 12 : return new byte[]{(byte)'G',(byte)'C'};
		case 13 : return new byte[]{(byte)'G',(byte)'G'};
		case 14 : return new byte[]{(byte)'G',(byte)'T'};
		case 15 : return new byte[]{(byte)'G',(byte)'N'};
		case 16 : return new byte[]{(byte)'T',(byte)'A'};
		case 17 : return new byte[]{(byte)'T',(byte)'C'};
		case 18 : return new byte[]{(byte)'T',(byte)'G'};
		case 19 : return new byte[]{(byte)'T',(byte)'T'};
		case 20 : return new byte[]{(byte)'T',(byte)'N'};
		case 21 : return new byte[]{(byte)'N',(byte)'A'};
		case 22 : return new byte[]{(byte)'N',(byte)'C'};
		case 23 : return new byte[]{(byte)'N',(byte)'G'};
		case 24 : return new byte[]{(byte)'N',(byte)'T'};
		case 25 : return new byte[]{(byte)'N',(byte)'N'};
		}
	throw new IllegalArgumentException("bad opcode ("+(int)opcode+")");
	}
private static byte bases2opcode(final char ref,final char alt) {
	switch(ref) {
		case 'A': switch(alt) {	case 'A': return 1;	case 'C': return 2;	case 'G': return 3;	case 'T': return 4;	case 'N': return 5; } break;
		case 'C': switch(alt) {	case 'A': return 6;	case 'C': return 7;	case 'G': return 8;	case 'T': return 9;	case 'N': return 10; } break;
		case 'G': switch(alt) {	case 'A': return 11;	case 'C': return 12;	case 'G': return 13;	case 'T': return 14;	case 'N': return 15; } break;
		case 'T': switch(alt) {	case 'A': return 16;	case 'C': return 17;	case 'G': return 18;	case 'T': return 19;	case 'N': return 20; } break;
		case 'N': switch(alt) {	case 'A': return 21;	case 'C': return 22;	case 'G': return 23;	case 'T': return 24;	case 'N': return 25; } break;
		}
	throw new IllegalArgumentException("bad ref/alt ("+ref+"/"+alt+")");
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
	
	final byte score_str[]=new byte[SCORE_STRLEN];
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
		if(f<0 || f>1.f)  throw new JvarkitException.FileFormatError("bad score "+tokens[4]);
		
		final ByteArrayOutputStream baos = new ByteArrayOutputStream(RECORD_SIZOF);
		final DataOutputStream daos = new DataOutputStream(baos);
		final int pos = Integer.parseInt(tokens[0]);
		
		if(prev_pos>pos )
			{
			LOG.error("input is not sorted "+prev_pos +" before "+pos);
			return -1;
			}
		prev_pos=pos;
		
		daos.writeInt(pos);
		daos.writeByte(bases2opcode(tokens[1].charAt(0),tokens[2].charAt(0)));
		daos.writeInt(Integer.parseInt(tokens[3].substring(4)));//after 'ENSG'
		
		
		Arrays.fill(score_str,(byte)'0');
		if(tokens[4]. equals("1")) tokens[4]="1.0";/* else padding with 0 gives large number '10000' */
		final byte b_array[] = tokens[4].getBytes();
		if(b_array.length>SCORE_STRLEN) throw new JvarkitException.FileFormatError("a(score)> "+SCORE_STRLEN+" in "+line);
		System.arraycopy(b_array,0,score_str,0,b_array.length);
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
