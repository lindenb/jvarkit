package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

@Deprecated /* use picard/RevertSam http://broadinstitute.github.io/picard/command-line-overview.html#RevertSam */
public class Biostar106668 extends AbstractBiostar106668
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar106668.class);

	
	
	
		@Override
		public Collection<Throwable> call() throws Exception
			{
			
			final List<String> args= getInputFiles();
			SamReader samReader=null;
			SAMRecordIterator iter=null;
			SAMFileWriter samWriter=null;
			long nConvert=0;
			try
				{
				SamFileReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
				if(args.isEmpty())
					{
					LOG.info("Reading sfomr stdin");
					samReader=SamFileReaderFactory.mewInstance().openStdin();
					}
				else if(args.size()==1)
					{
					File filename=new File(args.get(0));
					LOG.info("Reading from "+filename);
					samReader=SamFileReaderFactory.mewInstance().open(filename);
					}
				else
					{
					return wrapException("Illegal number of arguments.");
					}
				SAMFileHeader header=samReader.getFileHeader();
				header.addComment(getProgramCommandLine());
							
				SAMFileWriterFactory sfw=new SAMFileWriterFactory();
				sfw.setCreateIndex(false);
				sfw.setCreateMd5File(false);
				if(getOutputFile()!=null)
					{
					samWriter=sfw.makeSAMOrBAMWriter(header, true, getOutputFile());
					}
				else if(super.binary)
					{
					samWriter=sfw.makeBAMWriter(header,true, stdout());
					}
				else
					{
					samWriter=sfw.makeSAMWriter(header,true,stdout());
					}
				
				iter=samReader.iterator();
				SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
				while(iter.hasNext())
					{
					SAMRecord rec=progress.watch(iter.next());
					if(rec.getDuplicateReadFlag())
						{
						rec.setDuplicateReadFlag(false);
						++nConvert;
						}
					samWriter.addAlignment(rec);
					}
				progress.finish();
				LOG.info("Convert :"+nConvert);
				return Collections.emptyList();
				}
			catch (Exception e)
				{
				LOG.error(e);
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(samReader);
				CloserUtil.close(samWriter);
				}
			
			}
		
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException
		{
		new Biostar106668().instanceMainWithExit(args);
		}
		

	}
