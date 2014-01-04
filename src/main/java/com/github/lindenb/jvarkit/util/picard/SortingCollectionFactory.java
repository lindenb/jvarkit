package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import net.sf.picard.PicardException;
import net.sf.samtools.SAMFileWriterImpl;
import net.sf.samtools.util.SortingCollection;
import net.sf.samtools.util.SortingCollection.Codec;

public class SortingCollectionFactory<T>
	{
	private Comparator<T> comparator=null;
	private Class<T> componentType=null;
	private Codec<T> codec=null;
	private int maxRecordsInRAM = SAMFileWriterImpl.getDefaultMaxRecordsInRam();;
	private Collection<File> tmpDirs = Collections.emptyList();
	
	public void setComponentType(Class<T> componentType)
		{
		this.componentType = componentType;
		}
	
	public Class<T> getComponentType()
		{
		return componentType;
		}
	
	public void setComparator(Comparator<T> comparator)
		{
		this.comparator = comparator;
		}
	public Comparator<T> getComparator()
		{
		return comparator;
		}
	
	public void setCodec(Codec<T> codec)
		{
		this.codec = codec;
		}
	public Codec<T> getCodec()
		{
		return codec;
		}
	public void setMaxRecordsInRAM(int maxRecordsInRAM)
		{
		this.maxRecordsInRAM = maxRecordsInRAM;
		}
	
	public int getMaxRecordsInRAM()
		{
		return maxRecordsInRAM;
		}
	
	public void setTmpDirs(Collection<File> tmpDirs)
		{
		this.tmpDirs = tmpDirs;
		}
	
	public Collection<File> getTmpDirs()
		{
		return tmpDirs;
		}
	
	public SortingCollection<T> make()
		{
		if(componentType==null) throw new PicardException("componentType is undefined ");
		if(codec==null) throw new PicardException("codec is undefined ");
		if(comparator==null) throw new PicardException("comparator is undefined ");
		List<File> dirs=new ArrayList<File>(getTmpDirs());
		if(dirs.isEmpty())
			{
			dirs.add(new File(System.getProperty("java.io.tmpdir")));
			}
		return SortingCollection.newInstance(
				componentType,
				codec,
				comparator,
				maxRecordsInRAM,
				tmpDirs
				);
		}
	}
