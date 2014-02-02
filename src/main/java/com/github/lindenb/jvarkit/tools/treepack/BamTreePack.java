package com.github.lindenb.jvarkit.tools.treepack;

import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.io.IOUtils;
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
	
	
	 private void scan(SAMFileReader sfr)
		 {
		 sfr.setValidationStringency(ValidationStringency.SILENT);
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
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/BamTreePack";
		}
	@Override
	public String getProgramDescription()
		{
		return "Create a TreeMap from one or more SAM/BAM file. Ouput is a SVG file.";
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
					switch(handleOtherOptions(c, opt, null))
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
				info("Done stdin");
				}
			else
				{
				for(int optind =opt.getOptInd(); optind < args.length;++optind)
					{
					InputStream in=null;
					String filename= args[optind];
					info("Reading "+filename);
					if(IOUtils.isRemoteURI(filename))
						{
						in=IOUtils.openURIForReading(filename);
						sfr=new SAMFileReader(in);
						}
					else
						{
						sfr=new SAMFileReader(new File(filename));
						}
					scan(sfr);
					info("Done "+filename);
					CloserUtil.close(sfr);
					CloserUtil.close(in);
					
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
