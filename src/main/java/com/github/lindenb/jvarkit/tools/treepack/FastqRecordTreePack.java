package com.github.lindenb.jvarkit.tools.treepack;

import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.CloserUtil;


public class FastqRecordTreePack extends AbstractFastqRecordTreePack
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(FastqRecordTreePack.class);

	private class LengthNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(LengthNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(FastqRecord rec)
				{
				String key=String.valueOf(rec.getReadString().length());
				this.get(key).watch(rec);
				}
			}
		@Override
		public String getName() {
			return "length";
			}
		@Override
		public String getDescription()
			{
			return "Read Length";
			}
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	
	
	/** Group by Sample */
	private class SubSequenceFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(SubSequenceFactory.this,parent,label);
				}
			
			@Override
			public void watch(FastqRecord seq)
				{
				for(int i=0;i+ SubSequenceFactory.this.len <= seq.getReadString().length();++i)
					{
					String key=seq.getReadString().substring(i,i+SubSequenceFactory.this.len);
					this.get(key).watch(seq);
					}
				}
			}
		
		@Override
		public String getDescription()
			{
			return "DNA words. Length:"+len+" bases";
			}
		
		@Override
		public String getName() {
			return "word"+len;
			}
		
		private int len;
		SubSequenceFactory(int len)
			{
			this.len=len;
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	
	/** Group by InsertSize */
	private class QualityNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(QualityNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(FastqRecord rec)
				{
				String qual=rec.getBaseQualityString();
				int x=0;
				for(int i=0;i< qual.length();++i)
					{
					x+=SAMUtils.fastqToPhred(qual.charAt(i));
					}
				String s= intervalToString((int)(x/qual.length()),5);
				this.get(s).watch(rec);
				}
			}
		
		@Override
		public String getDescription()
			{
			return "Mean quality (phred>=33)";
			}
		
		@Override
		public String getName() {
			return "qual";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	
	

	
	
	private FastqRecordTreePack()
		{
		
		}
	
	
	 private void scan(FastqReader iter)
		 {
		 while(iter.hasNext())
			 {
			 FastqRecord rec=iter.next();
			 root.watch(rec);
			 }
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
			for(int i=1;i<=5;++i) _factories.add(new SubSequenceFactory(i));
			_factories.add(new LengthNodeFactory());
			_factories.add(new QualityNodeFactory());
			}
		return _factories;
		}


	@Override
	public Collection<Throwable> call() throws Exception {
		if(super.listFactories) 
			{
			super.printAvailableFactories();
			return RETURN_OK;
			}
		setDimension(super.dimensionStr);
		buildFactoryChain(super.chainExpression);
		if(super.nodeFactoryChain.isEmpty())
			{
			return wrapException("no path defined");
			}
		
		FastqReader fqr=null;
		final List<String> args = getInputFiles();
		try
			{
			if(args.isEmpty())
				{
				LOG.info("Reading stdin");
				fqr=new FourLinesFastqReader(stdin());
				scan(fqr);
				CloserUtil.close(fqr);
				LOG.info("Done stdin");
				}
			else
				{
				for(final String filename:args)
					{
					InputStream in=null;
					LOG.info("Reading "+filename);
					if(IOUtils.isRemoteURI(filename))
						{
						in=IOUtils.openURIForReading(filename);
						fqr=new FourLinesFastqReader(in);
						}
					else
						{
						fqr=new FourLinesFastqReader(new File(filename));
						}
					scan(fqr);
					LOG.info("Done "+filename);
					CloserUtil.close(fqr);
					CloserUtil.close(in);
					}
				}
			
			this.layout();
			this.svg();
			return RETURN_OK;
			}
		catch (Exception e)
			{
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(fqr);
			}
		}
	
	public static void main(String[] args)
		{
		new FastqRecordTreePack().instanceMainWithExit(args);
		}

	}
