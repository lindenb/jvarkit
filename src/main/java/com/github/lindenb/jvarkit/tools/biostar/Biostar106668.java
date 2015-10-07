package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

@Deprecated /* use picard/RevertSam http://broadinstitute.github.io/picard/command-line-overview.html#RevertSam */
public class Biostar106668 extends AbstractBiostar106668
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar106668.class);
	
	@Override
	public Command createCommand()
		{
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBiostar106668.AbstractBiostar106668Command
		{
		@Override
		protected Collection<Throwable> call(final String inputName) throws Exception {
			
			SamReader samReader=null;
			SAMRecordIterator iter=null;
			SAMFileWriter samWriter=null;
			long nConvert=0;
			try
				{
				samReader =  openSamReader(inputName);
				SAMFileHeader header=samReader.getFileHeader();
				samWriter= openSAMFileWriter(header, true);
				
				iter=samReader.iterator();
				SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samReader.getFileHeader().getSequenceDictionary());
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
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(samReader);
				CloserUtil.close(samWriter);
				}

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
