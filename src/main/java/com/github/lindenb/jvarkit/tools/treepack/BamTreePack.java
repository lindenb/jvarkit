package com.github.lindenb.jvarkit.tools.treepack;

import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;


public class BamTreePack extends AbstractBamTreePack
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BamTreePack.class);

	/** Group by SAMReadGroupRecord */
	private class GroupIdNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(GroupIdNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(SAMRecord seq)
				{
				SAMReadGroupRecord g=seq.getReadGroup();
				String gid=(g==null?"*":g.getReadGroupId());
				if(gid==null || gid.isEmpty()) gid="*";
				this.get(gid).watch(seq);
				}
			}
		@Override
		public String getName() {
			return "group";
			}
		@Override
		public String getDescription()
			{
			return "SAMReadGroupRecord.groupId)";
			}
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	
	
	/** Group by Sample */
	private class SampleNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(SampleNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(SAMRecord seq)
				{
				SAMReadGroupRecord g=seq.getReadGroup();
				String sample=(g==null?"*":g.getSample());
				if(sample==null || sample.isEmpty()) sample="*";
				this.get(sample).watch(seq);
				}
			}
		
		@Override
		public String getDescription()
			{
			return "Sample name";
			}
		
		@Override
		public String getName() {
			return "sample";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	
	/** Group by InsertSize */
	private class InsertSizeNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(InsertSizeNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(SAMRecord rec)
				{
				String s="*";
				int v=Math.abs(rec.getInferredInsertSize());
				if(v>0) s=intervalToString(v, 10);
					
				this.get(s).watch(rec);
				}
			}
		
		@Override
		public String getDescription()
			{
			return "Insert-Size";
			}
		
		@Override
		public String getName() {
			return "insertsize";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	/** Group by Chromosome */
	private class ChromosomeNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(ChromosomeNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(SAMRecord seq)
				{
				String chr=seq.getReferenceName();
				if(chr==null) chr=SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
				this.get(chr).watch(seq);
				}
			}
		@Override
		public String getDescription()
			{
			return "Chromosome";
			}
		@Override
		public String getName()
			{
			return "chrom";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	/** Group by MAPQ */
	private class MappingQualityFactory
		extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(MappingQualityFactory.this,parent,label);
				}
			
			@Override
			public void watch(SAMRecord seq)
				{
				String qualStr="*";
				if(!seq.getReadUnmappedFlag())
					{
					int qual=seq.getMappingQuality();
					if(qual==255)
						{
						qualStr="?";
						}
					else if(qual==0)
						{
						qualStr="*";
						}
					else
						{
						qual=((int)(qual/10.0));
						qualStr="["+(int)(qual)*10+"-"+((qual+1)*10)+"[";
						}
					}
				this.get(qualStr).watch(seq);
				}
			}
		@Override
		public String getDescription()
			{
			return "Mapping quality window.size=10";
			}
		@Override
		public String getName() {
			return "mapq";
			}

		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}

	/** Group by MAPQ */
	private class PropertyPairedFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(PropertyPairedFactory.this,parent,label);
				}
			
			@Override
			public void watch(SAMRecord seq)
				{
				String s="not-paired";
				if(seq.getReadPairedFlag())
					{
					s=(seq.getReadPairedFlag()?"properly-paired":"not-property-paired");
					}
				this.get(s).watch(seq);
				}
			}
		@Override
		public String getDescription()
			{
			return "Properly paired reads";
			}
		@Override
		public String getName() {
			return "properly-paired";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	/** Group by NotSameChromNodeFactory */
	private class NotSameChromNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(NotSameChromNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(SAMRecord rec)
				{
				String ref1=null;
				String ref2=null;
				if(rec.getReadPairedFlag() && 
					!rec.getReadUnmappedFlag() &&
					!rec.getMateUnmappedFlag() &&
					(ref1=rec.getReferenceName())!=(ref2=rec.getMateReferenceName())
					)
					{
					String s=ref1.compareTo(ref2)<0?
							ref1+"/"+ref2:
							ref2+"/"+ref1
							;
					this.get(s).watch(rec);
					}
				
				}
			}
		
		@Override
		public String getDescription()
			{
			return "read chromosome != mate chromosome";
			}
		
		@Override
		public String getName() {
			return "splitchrom";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}

	
	private class SamFlagNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(SamFlagNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(SAMRecord rec)
				{
				String s=String.valueOf(rec.getFlags());
				this.get(s).watch(rec);	
				}
			}
		
		@Override
		public String getDescription()
			{
			return "sam flag";
			}
		
		@Override
		public String getName() {
			return "samflag";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	

	
	
	private BamTreePack()
		{
		
		}
	
	
	 private void scan(SamReader sfr)
		 {
		 SAMRecordIterator iter=sfr.iterator();
		 SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(sfr.getFileHeader().getSequenceDictionary());
		 while(iter.hasNext())
			 {
			 SAMRecord rec=iter.next();
			 progress.watch(rec);
			 root.watch(rec);
			 }
		 progress.finish();
		 }
	  
	@Override
	public void printOptions(PrintStream out)
		{
		super.printOptions(out);
		}
	
	
	private  List<NodeFactory> _factories=null;
	@Override
	protected List<NodeFactory> getAllAvailableFactories()
		{
		if(_factories==null)
			{
			_factories=new ArrayList<NodeFactory>();
			_factories.add(new ChromosomeNodeFactory());
			_factories.add(new SampleNodeFactory());
			_factories.add(new PropertyPairedFactory());
			_factories.add(new MappingQualityFactory());
			_factories.add(new GroupIdNodeFactory());
			_factories.add(new InsertSizeNodeFactory());
			_factories.add(new NotSameChromNodeFactory());
			_factories.add(new SamFlagNodeFactory());
			}
		return _factories;
		}
	

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBamTreePack.AbstractBamTreePackCommand
		{    
		@Override
		public Collection<Throwable> call() throws Exception
			{
			if(super.nodeFactoryChain.isEmpty())
				{
				LOG.error("no path defined");
				return -1;
				}
			List<String> args= getInputFiles();
			SamReader sfr=null;
			try
				{
				if(args.isEmpty())
					{
					LOG.info("Reading stdin");
					sfr= openSamReader(null);
					scan(sfr);
					CloserUtil.close(sfr);
					LOG.info("Done stdin");
					}
				else
					{
					for(String filename:args)
						{
						InputStream in=null;
						LOG.info("Reading "+filename);
						sfr= openSamReader(filename);
						scan(sfr);
						LOG.info("Done "+filename);
						CloserUtil.close(sfr);
						CloserUtil.close(in);
						
						}
					}
				
				this.layout();
				this.svg();
				return RETURN_OK;
				}
			catch (Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(sfr);
				}
			}
		}
	
	public static void main(String[] args)
		{
		new BamTreePack().instanceMainWithExit(args);
		}

	}
