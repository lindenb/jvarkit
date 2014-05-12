package com.github.lindenb.jvarkit.util.ucsc;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.samtools.util.CloserUtil;

/**
 * PSL format align.
 * see http://genome.ucsc.edu/FAQ/FAQformat.html#format2
 */
public class PslAlign
	{
	private int matches=0;// - Number of bases that match that aren't repeats
    private int misMatches=0;// - Number of bases that don't match
    private int repMatches=0;// - Number of bases that match but are part of repeats
    private int nCount=0;// - Number of 'N' bases
    private int qNumInsert=0;// - Number of inserts in query
    private int qBaseInsert=0;// - Number of bases inserted in query
    private int tNumInsert=0;// - Number of inserts in target
    private int tBaseInsert=0;// - Number of bases inserted in target
    private char  strand='+';// -9 '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
    private String qName=null;// 9- Query sequence name
    private int qSize=0;// -10 Query sequence size
    private int qStart=0;//-11 Alignment start position in query
    private int qEnd=0;// -12 Alignment end position in query
    private String tName=null;//-13 Target sequence name
    private int tSize=0;// -14 Target sequence size
    private int tStart=0;// -15 Alignment start position in target
    private int tEnd=0;// -16 Alignment end position in target
    //int blockCount;//17 - Number of blocks in the alignment (a block contains no gaps)
    //int blockSizes[];//18 - Comma-separated list of sizes of each block
    //int qStarts[];//19 - Comma-separated list of starting positions of each block in query
    //int tStarts[];//20 - Comma-separated list of starting positions of each block in target 
    private List<Block> blocks=new ArrayList<Block>();
    
    public class Block
    	{
    	private int qStart;
    	private int tStart;
    	private int len;
    	
    	public int getQStart()
    		{
			return qStart;
			}
    	public int getTStart()
			{
			return tStart;
			}
    	
    	/** Alignment end position in block 0-based */
    	public int getTargetEnd()
			{
			return getTStart()+size();
			}
    	
    	public int size()
    		{
    		return this.len;
    		}
    	
    	@Override
    	public String toString() {
    		return qStart+":"+tStart+":"+len;
    		}
    	}
    
    private static final Pattern COMMA=Pattern.compile("[,]");
   
    public PslAlign()
    	{
    	
    	}

    public PslAlign(String tokens[])
    	{
    	strand=tokens[8].charAt(0);
    	qName=tokens[9];
    	qStart=Integer.parseInt(tokens[11]);
    	qEnd=Integer.parseInt(tokens[12]);
    	tName=tokens[13];
    	tStart=Integer.parseInt(tokens[15]);
    	tEnd=Integer.parseInt(tokens[16]);
    	int blockCount=Integer.parseInt(tokens[17]);
    	for(int i=0;i< blockCount;++i)
    		{
        	this.blocks.add(new Block());
    		}
    	
    	String ss[]=COMMA.split(tokens[18]);
    	for(int i=0;i< blockCount;++i) this.blocks.get(i).len=Integer.parseInt(ss[i]);
    	ss=COMMA.split(tokens[19]);
    	for(int i=0;i< blockCount;++i) this.blocks.get(i).qStart=Integer.parseInt(ss[i]);
    	ss=COMMA.split(tokens[20]);
    	for(int i=0;i< blockCount;++i) this.blocks.get(i).tStart=Integer.parseInt(ss[i]);
    	}
    
    public void addBlock(int qStart,int tStart,int len)
    	{
    	Block b=new Block();
    	b.qStart=qStart;
    	b.tStart=tStart;
    	b.len=len;
    	this.blocks.add(b);
    	}
    
    
    public String toString()
    	{
    	StringBuilder b=new StringBuilder();
    	b.append(getMatches()).append('\t');
    	b.append(getMisMatches()).append('\t');
    	b.append(getRepMatches()).append('\t');
    	b.append(getNCount()).append('\t');
    	b.append(getQNumInsert()).append('\t');
    	b.append(getQBaseInsert()).append('\t');
    	b.append(getTNumInsert()).append('\t');
    	b.append(getTBaseInsert()).append('\t');
    	b.append(getStrand()).append('\t');
    	b.append(getQName()).append('\t');
    	b.append(getQSize()).append('\t');
    	b.append(getQStart()).append('\t');
    	b.append(getQEnd()).append('\t');
    	b.append(getTName()).append('\t');
    	b.append(getTSize()).append('\t');
    	b.append(getTStart()).append('\t');
    	b.append(getTEnd()).append('\t');
    	b.append(getBlockCount()).append('\t');
    	
    	//
    	for(Block ck:this.blocks)
    		{
    		b.append(ck.len);
    		b.append(',');
    		}
    	b.append('\t');
    	
    	//
    	for(Block ck:this.blocks)
			{
			
			b.append(ck.qStart);
			b.append(',');
			}
		b.append('\t');
		
    	//
    	for(Block ck:this.blocks)
			{
			b.append(ck.tStart);
			b.append(',');
			}

		
    	return b.toString();
    	}
    
    @Deprecated
    public String getQueryName()
    	{
    	return qName;
    	}
    @Deprecated
    public String getTargetName()
    	{
		return tName;
		}
    @Deprecated
    public int getTargetStart()
    	{
		return tStart;
		}
    @Deprecated
    public int getTargetEnd()
    	{
		return tEnd;
		}
    
    @Deprecated
    /* Alignment start position in query 0-based */
	public int getQueryStart()
		{
		return qStart;
		}
    @Deprecated
	/* Alignment end position in query 0-based */
	public int getQueryEnd()
		{
		return qEnd;
		}
	
	public int getBlockCount()
		{
		return this.blocks.size();
		}
    
    public List<Block> getBlocks()
    	{
    	return this.blocks;
    	}
    
    public Block getBlock(int i)
		{
		return this.blocks.get(i);
		}
    
    
    
    
    public int getMatches() {
		return matches;
	}

	public void setMatches(int matches) {
		this.matches = matches;
	}

	public int getMisMatches() {
		return misMatches;
	}

	public void setMisMatches(int misMatches) {
		this.misMatches = misMatches;
	}

	public int getRepMatches() {
		return repMatches;
	}

	public void setRepMatches(int repMatches) {
		this.repMatches = repMatches;
	}

	public int getNCount() {
		return nCount;
	}

	public void setNCount(int nCount) {
		this.nCount = nCount;
	}

	public int getQNumInsert() {
		return qNumInsert;
	}

	public void setQNumInsert(int qNumInsert) {
		this.qNumInsert = qNumInsert;
	}

	public int getQBaseInsert() {
		return qBaseInsert;
	}

	public void setQBaseInsert(int qBaseInsert) {
		this.qBaseInsert = qBaseInsert;
	}

	public int getTNumInsert() {
		return tNumInsert;
	}

	public void setTNumInsert(int tNumInsert) {
		this.tNumInsert = tNumInsert;
	}

	public int getTBaseInsert() {
		return tBaseInsert;
	}

	public void setTBaseInsert(int tBaseInsert) {
		this.tBaseInsert = tBaseInsert;
	}

	public char getStrand() {
		return strand;
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}

	public String getQName() {
		return qName;
	}

	public void setQName(String qName) {
		this.qName = qName;
	}

	public int getQSize() {
		return qSize;
	}

	public void setQSize(int qSize) {
		this.qSize = qSize;
	}

	public int getQStart() {
		return qStart;
	}

	public void setQStart(int qStart) {
		this.qStart = qStart;
	}

	public int getQEnd() {
		return qEnd;
	}

	public void setQEnd(int qEnd) {
		this.qEnd = qEnd;
	}

	public String getTName() {
		return tName;
	}

	public void setTName(String tName) {
		this.tName = tName;
	}

	public int getTSize() {
		return tSize;
	}

	public void setTSize(int tSize) {
		this.tSize = tSize;
	}

	public int getTStart() {
		return tStart;
	}

	public void setTStart(int tStart) {
		this.tStart = tStart;
	}

	public int getTEnd() {
		return tEnd;
	}

	public void setTEnd(int tEnd) {
		this.tEnd = tEnd;
	}

	
    
    public static Iterator<PslAlign> iterator(LineIterator iter)
    	{
    	return new MyIter(iter);
    	}
    
    private static class MyIter implements Iterator<PslAlign>
    	{
    	private Pattern tab=Pattern.compile("[\t]");
    	private LineIterator delegate=null;
    	MyIter(LineIterator delegate)
    		{
    		this.delegate=delegate;
    		}
    	@Override
    	public boolean hasNext()
    		{
    		if(delegate==null) return false;
    		while(delegate.hasNext())
    			{
    			String line=delegate.peek();
    			if(line.isEmpty() || !Character.isDigit(line.charAt(0)))//PSL header, empty lines...
    				{
    				delegate.next();
    				continue;
    				}
    			return true;
    			}
    		CloserUtil.close(delegate);
    		delegate=null;
    		return false;
    		}
    	@Override
    	public PslAlign next()
    		{
    		if(!hasNext()) throw new IllegalStateException();
    		return new PslAlign(this.tab.split(delegate.next()));
    		}
    	
    	@Override
    	public void remove() {
    		throw new UnsupportedOperationException();
    		}
    	
    	@SuppressWarnings("unused")
		public void close()
    		{
    		CloserUtil.close(delegate);
    		}
    	}
    
	}
