package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMTextHeaderCodec;
import net.sf.samtools.util.BufferedLineReader;

/** utility to load a SAMSequenceDictionary */
public class SAMSequenceDictionaryFactory 
{
public SAMSequenceDictionaryFactory()
	{
	}
/** get a SAMSequenceDictionary from a file. try to find the right file if sequence is a fasta file */
public  SAMSequenceDictionary load(File file) throws IOException
	{
	final File dictionary = findSequenceDictionary(file); 
	if(dictionary==null) throw new FileNotFoundException("Cannot find dict file for "+file+". Was the reference sequence indexed with Picard ?");
	IoUtil.assertFileIsReadable(dictionary);
	BufferedLineReader blr=null;
    try {
        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        blr=new BufferedLineReader(new FileInputStream(dictionary));
        final SAMFileHeader header = codec.decode(blr, dictionary.toString());
        if (header.getSequenceDictionary() != null && header.getSequenceDictionary().size() > 0)
        	{
            return header.getSequenceDictionary();
        	}
        else
        	{
        	 throw new PicardException("no sequence dictionary file in  " + dictionary);
        	}
    	}
    catch (Exception e) {
        throw new PicardException("Could not open sequence dictionary file: " + dictionary, e);
    }
    finally
    	{
    	if(blr!=null) blr.close();
    	}
	}

/** copied from picard: AbstractFastaSequenceFile . Find the dict or return null */
public File findSequenceDictionary(final File file) {
	if(!file.exists()) return null;
	if(!file.isFile()) return null;
	String filename=file.getName();
    for (final String extension : ReferenceSequenceFileFactory.FASTA_EXTENSIONS)
    	{
        if (!filename.endsWith(extension)) continue;
        File dictFile=new File(file.getParentFile(),filename+".dict");
        if(dictFile.exists()) return dictFile;
        dictFile=new File(file.getParentFile(),filename.substring(0, filename.lastIndexOf(extension))+".dict");
        if(dictFile.exists()) return dictFile;
    	}
    //the index was given
    if(filename.endsWith(".dict")) return file;
    return null;
	}

}
