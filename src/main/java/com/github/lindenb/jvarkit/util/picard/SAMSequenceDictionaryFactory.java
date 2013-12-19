package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
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
	if(dictionary==null)
		{
		final File dict2=findSequenceFaidx(file);
		if(dict2==null) throw new FileNotFoundException("Cannot find dict file for "+file+". Was the reference sequence indexed with Picard ?");
		return loadFaidxSequenceDict(dict2);
		}
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


/**  find fasta + faidx */
public File findSequenceFaidx(final File file)
	{
	final String faidx_extension=".fai";
	if(!file.exists()) return null;
	if(!file.isFile()) return null;
	String filename=file.getName();
    for (final String extension : ReferenceSequenceFileFactory.FASTA_EXTENSIONS)
    	{
        if (!filename.endsWith(extension)) continue;
        File dictFile=new File(file.getParentFile(),filename+faidx_extension);
        if(dictFile.exists()) return dictFile;
    	}
    //the index was given
    if(filename.endsWith(faidx_extension)) return file;
    return null;
	}

private SAMSequenceDictionary loadFaidxSequenceDict(File faidx)
	throws IOException
	{
	Set<String> seen=new HashSet<String>();
	Pattern tab=Pattern.compile("[\t]");
	SAMSequenceDictionary dict=null;
	BufferedLineReader blr=null;
    try {
        blr=new BufferedLineReader(new FileInputStream(faidx));
		List<SAMSequenceRecord> L=new ArrayList<SAMSequenceRecord>();
		String line;
		while((line=blr.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			String tokens[]=tab.split(line);
			if(tokens.length<5 || tokens[0].isEmpty() || tokens[1].isEmpty()) 
				{
				throw new PicardException("Bad line in  file: " + faidx+ " "+line.replaceAll("\t", "\\t"));
				}
			
			if(seen.contains(tokens[0])) 
				{
				throw new PicardException("Duplicate sequence name in  file: " + faidx+ " "+line.replaceAll("\t", "\\t"));
				}
			seen.add(tokens[0]);
			SAMSequenceRecord ssr=new SAMSequenceRecord(tokens[0], Integer.parseInt(tokens[1]));
			L.add(ssr);
			}
		dict=new SAMSequenceDictionary(L);
		}
	catch(Exception err)
		{
		throw new PicardException("Could not read sequence dictionary file: " + faidx, err);
		}
	finally
		{
		if(blr!=null) blr.close();
		}
	 if(dict!= null && dict.size() > 0)
	  	{
	    return dict;
	  	}
	  else
	  	{
	  	throw new PicardException("no sequence dictionary file in  " + faidx);
	  	}
	}
}
