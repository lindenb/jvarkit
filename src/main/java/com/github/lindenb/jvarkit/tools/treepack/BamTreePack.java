package com.github.lindenb.jvarkit.tools.treepack;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloserUtil;


public class BamTreePack extends AbstractTreePackCommandLine<SAMRecord>
	{
	

	
	
	
	/** Group by Sample */
	private class SampleNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent)
				{
				super(SampleNodeFactory.this,parent);
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
		public String getName() {
			return "sample";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent)
			{
			return new MyNode(parent);
			}
		}
	
	/** Group by Chromosome */
	private class ChromosomeNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent)
				{
				super(ChromosomeNodeFactory.this,parent);
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
		public String getName()
			{
			return "chrom";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent)
			{
			return new MyNode(parent);
			}
		}
	
	/** Group by MAPQ */
	private class MappingQualityFactory
		extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent)
				{
				super(MappingQualityFactory.this,parent);
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
					else
						{
						qual=((int)(qual/10.0))*10;
						qualStr=String.valueOf(qual);
						}
					}
				this.get(qualStr).watch(seq);
				}
			}
		@Override
		public String getName() {
			return "mapq";
			}

		public BranchNode createBranch(BranchNode parent)
			{
			return new MyNode(parent);
			}
		}

	/** Group by MAPQ */
	private class PropertyPairedFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent)
				{
				super(PropertyPairedFactory.this,parent);
				}
			
			@Override
			public void watch(SAMRecord seq)
				{
				String s="";
				if(seq.getReadPairedFlag())
					{
					s=(seq.getReadPairedFlag()?"properly-paired":"not-property-paired");
					}
				this.get(s).watch(seq);
				}
			}
		@Override
		public String getName() {
			return "properly-paired";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent)
			{
			return new MyNode(parent);
			}
		}
	
	
	
	
	

	
	
	private BamTreePack()
		{
		
		}
	
	
	 private void scan(SAMFileReader sfr)
		 {
		 sfr.setValidationStringency(ValidationStringency.LENIENT);
		 SAMRecordIterator iter=sfr.iterator();
		 while(iter.hasNext())
			 {
			 root.watch(iter.next());
			 }
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
		return L;
		}
	
	@Override
	public int doWork(String[] args)
		{
		String path=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "e:"))!=-1)
			{
			switch(c)
				{
				case 'e': path=opt.getOptArg();break;
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
		if(buildFactoryChain(path)!=0)
			{
			error("Cannot parse "+path);
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
