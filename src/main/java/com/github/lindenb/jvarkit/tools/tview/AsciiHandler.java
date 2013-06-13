package com.github.lindenb.jvarkit.tools.tview;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.Writer;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;

public class AsciiHandler extends AbstractHandler
	{
	private boolean printReadName=false;
	private PrintWriter out;
	public AsciiHandler()
		{
		this(System.out);
		}
	
	public AsciiHandler(OutputStream out)
		{
		this.out=new PrintWriter(out);
		}
	public AsciiHandler(Writer out)
		{
		this.out=new PrintWriter(out);
		}
	
	public void endRow() {	this.out.println(); }
	public void whiteSpace() { this.out.print(' ');}
	public void deletion() { this.out.print('*');}
	public void base(
			SAMRecord rec,int readPos,
			IndexedFastaSequenceFile ref,int refPos
			)
		{
		
		char c1;
		if(printReadName)
			{
			c1=(rec.getReadNegativeStrandFlag()?',':'.');
			String readName=rec.getReadName();
			if(readPos<readName.length())
				{
				c1=readName.charAt(readPos);
				c1=(rec.getReadNegativeStrandFlag()?Character.toLowerCase(c1):Character.toUpperCase(c1));
				}
			}
		else
			{
			c1=(char)rec.getReadBases()[readPos];
			Character c2=getReferenceBaseAt(ref,rec.getReferenceName(),refPos);
			if(c2!=null)
				{
				if(c1==c2) c1=(rec.getReadNegativeStrandFlag()?',':'.');
				}
			}
		this.out.print(c1);
		}

	@Override
	public void reference(IndexedFastaSequenceFile ref, String seqName,int refPos)
		{
		if(ref==null)
			{
			out.print('N');
			return;
			}
		if(refPos==-1)
			{
			out.print('*');
			return;
			}
		Character c=getReferenceBaseAt(ref,seqName,refPos);
		if(c==null)
			{
			out.print('#');
			return;
			}
		out.print(c);
		}
	
	@Override
	public void endReferenceSeq() {
		this.endRow();
		}
	
	@Override
	public void endDocument()
		{
		this.out.println();
		this.out.flush();
		}
	}
