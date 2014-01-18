package com.github.lindenb.jvarkit.tools.treepack;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloserUtil;


public class BamTreePack extends AbstractTreePackCommandLine<SAMRecord>
	{
	
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
	
	
	
	
	

	
	
	private BamTreePack()
		{
		
		}
	
	
	 private void scan(SAMFileReader sfr)
		 {
		 sfr.setValidationStringency(ValidationStringency.LENIENT);
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
	
	
	
	@Override
	protected List<NodeFactory> getAllAvailableFactories() {
		List<NodeFactory> L=new ArrayList<NodeFactory>();
		L.add(new ChromosomeNodeFactory());
		L.add(new SampleNodeFactory());
		L.add(new PropertyPairedFactory());
		L.add(new MappingQualityFactory());
		L.add(new GroupIdNodeFactory());
		return L;
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()))!=-1)
			{
			switch(c)
				{
				default: 
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(super.nodeFactoryChain.isEmpty())
			{
			error("no path defined");
			return -1;
			}
		
		SAMFileReader sfr=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading stdin");
				sfr=new SAMFileReader(System.in);
				scan(sfr);
				CloserUtil.close(sfr);
				}
			else
				{
				for(int optind =opt.getOptInd(); optind < args.length;++optind)
					{
					File f=new File(args[optind]);
					info("Reading stdin");
					sfr=new SAMFileReader(f);
					scan(sfr);
					CloserUtil.close(sfr);
					}
				}
			this.layout();
			this.svg();
			return 0;
			}
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
			}
		}
	
	public static void main(String[] args)
		{
		new BamTreePack().instanceMainWithExit(args);

		}

	}
