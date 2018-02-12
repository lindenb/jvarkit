package com.github.lindenb.jvarkit.tools.tests;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;

public class TestUtils {

protected static class ParamCombiner
	{
	private List<List<Object>> data = new ArrayList<>();
	
	public ParamCombiner() {
		
	}
	
	public ParamCombiner initList(Object list[]) {
		this.data.clear();
		for(final Object o:list)
			{
			List<Object> row = new ArrayList<>(1);
			row.add(o);
			this.data.add(row);
			}
		return this;
		}
	
	ParamCombiner append(Object right)
		{
		for(int i=0;i< this.data.size();i++)
			{
			this.data.get(i).add(right);
			}
		return this;
		}
	
	public ParamCombiner when(Function<List<Object>, Object> fun)
		{
		for(int i=0;i< this.data.size();i++)
			{
			this.data.get(i).add(fun.apply(this.data.get(i)));
			}
		return this;
		}
	
	public ParamCombiner filter(Predicate<List<Object>> pred)
		{
		this.data.removeIf(pred);
		return this;
		}
	
	public ParamCombiner product(Object...right)
		{
		final List<List<Object>> data2 = new ArrayList<>();
		
		for(int i=0;i< this.data.size();i++)
			{
			final List<Object> datarow = this.data.get(i);
			for(int j=0;j< right.length;j++)
				{
				final List<Object> row2 = new ArrayList<>(datarow);
				row2.add(right[j]);
				data2.add(row2);
				}
			}
		this.data.clear();
		this.data.addAll(data2);
		return this;
		}
	
	public Object[][] build() {
		final Object[][] table = new Object[this.data.size()][];
		for(int i=0;i< this.data.size();i++)
			{
			table[i] = this.data.get(i).toArray(new Object[this.data.get(i).size()]);
			}
		return table;
		}
	}
	
private List<File> _collectFiles(final File dir, FilenameFilter filter) {
	List<File> list = new ArrayList<>();
	if(dir==null || !dir.isDirectory()) return list;
	File child[] = dir.listFiles(filter);
	if(child==null || child.length==0) return list;
	for(File c : child) {
		if(c==null) continue;
		if(c.isDirectory())
			{
			list.addAll(_collectFiles(c,filter));
			continue;
			}
		if(!c.canRead() ) continue;
		if(!c.isFile()) continue;
		
		if(filter.accept(dir, c.getName())) {
			list.add(c);
			}
		}
	return list;
	}

protected Object[] _collectAllFiles( FilenameFilter filter) {
	return _collectFiles(new File("./src/test/resources/"),filter).
			stream().
			map(F->(Object)F.getPath()).
			toArray(i->new Object[i]);
		}

protected Object[] collectAllVcfs() {
	return _collectAllFiles((D,N)->N.endsWith(".vcf") || N.endsWith(".vcf.gz"));
	}

protected Object[] collectAllSamOrBam() {
	return _collectAllFiles((D,N)->N.endsWith(".sam") || N.endsWith(".bam"));
	}

protected Object[] collectAllFasta() {
	return _collectAllFiles((D,N)->N.endsWith(".fa") || N.endsWith(".fasta"));
	}

}
