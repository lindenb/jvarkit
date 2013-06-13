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
	private StringBuilder ref1=null;//new StringBuilder(); 
	private StringBuilder ref2=null;//new StringBuilder(); 
	@Override
	public void beginReferenceSeq() {
		ref2=new StringBuilder();
		ref1=new StringBuilder();
		}
	@Override
	public void reference(IndexedFastaSequenceFile ref, String seqName,int refPos)
		{
		while(ref2.length()< ref1.length())
			{
			ref2.append(' ');
			}
		if(ref==null)
			{
			ref1.append('N');
			}
		else if(refPos==-1)
			{
			ref1.append('*');
			}
		else
			{
			Character c=getReferenceBaseAt(ref,seqName,refPos);
			if(c==null)
				{
				ref1.append('#');
				}
			else
				{
				ref1.append(c);
				}
			}
		
		if(refPos!=-1 && refPos%10==0)
			{
			ref2.append(String.valueOf(refPos));
			}
		}
	
	
	@Override
	public void endReferenceSeq()
		{
		this.out.print(ref2);
		this.endRow();
		this.out.print(ref1);
		this.endRow();
		ref1=null;
		ref2=null;
		}
	
	@Override
	public void endDocument()
		{
		this.out.println();
		this.out.flush();
		}
	}
