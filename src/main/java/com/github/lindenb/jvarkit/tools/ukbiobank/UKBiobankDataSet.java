package com.github.lindenb.jvarkit.tools.ukbiobank;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;

import htsjdk.samtools.util.StringUtil;

public class UKBiobankDataSet {
private final Map<String,List<Coding>> codings = new HashMap<>();
private final Map<String,DataDict> column2datadict = new HashMap<>();
class Coding   {
	private final Map<String,String> hash;
	Coding(Map<String,String> hash) {
		this.hash = hash;
		}
	String getCodingName() {
		return this.hash.get("coding_name");
		}
	String getCode() {
		return this.hash.get("code");
		}
	String getMeaning() {
		return this.hash.get("meaning");
		}
	String getParentCode() {
		return this.hash.get("parent_code");
		}
	@Override
	public int hashCode() {
		return getCodingName().hashCode()*31 + getCode().hashCode();
		}
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof Coding)) return true;
		final Coding o=Coding.class.cast(obj);
		return o.getCodingName().equals(this.getCodingName()) &&
				o.getCode().equals(this.getCode());
		}
	public Set<Coding> _getAllDescendants(final Set<Coding> set) {
		if(set.contains(this)) return set;
		set.add(this);
		if(this.getCode().equals(this.getParentCode())) return set;
		if(this.getParentCode().isEmpty()) return set;
		for(Coding dd : codings.get(this.getCodingName())) {
			if(dd.getCode().equals(this.getParentCode())) {
				dd._getAllDescendants(set);
				}
			}
		return set;
		}
	public Set<Coding> getAllDescendants() {
		return _getAllDescendants(new HashSet<>());
		}
	}

class DataDict   {
	private final Map<String,String> hash;
	DataDict(Map<String,String> hash) {
		this.hash = hash;
		}
	String getEntity() {
		return this.hash.get("entity");
		}
	String getName() {
		return this.hash.get("name");
		}
	String getCodingName() {
		return this.hash.get("coding_name");
		}
	List<Coding> getCodings() {
		final String c = getCodingName();
		if(StringUtil.isBlank(c)) return Collections.emptyList();
		return UKBiobankDataSet.this.codings.getOrDefault(c,Collections.emptyList());
		}
	
	String getTitle() {
		return this.hash.get("title");
		}
	String getLinkOut() {
		return this.hash.get("linkout");
		}
	String getFolderPath() {
		return this.hash.get("folder_path");
		}
	String getReferencedEntityField() {
		return this.hash.get("referenced_entity_field");
		}
	String getRelationShip() {
		return this.hash.get("relationship");
		}
	String getUnits() {
		return this.hash.get("units");
		}
	
	private String translate(String s) {
		for(Coding  coding : getCodings()) {
			if(coding.getCode().equals(s)) {
				return coding.getMeaning();
				}
			}
		return s;
		}
	
	String format1(String s) {
		if(StringUtil.isBlank(s)) return "";
		s = translate(s);
		if(!StringUtil.isBlank(getUnits())) {
			s+="("+getUnits()+")";
			}
		return s;
		}
	
	List<String> getFolderPaths() {
		String s = getFolderPath();
		List<String> L = new ArrayList<>();
		for(;;) {
			int x = s.indexOf(" > ");
			if(x==-1) {
				L.add(s);
				break;
				}
			else
				{
				L.add(s.substring(0,x));
				s=s.substring(x+3);
				}
			}
		return L;
		}
	boolean isMultiSelect() {
		return this.hash.get("is_multi_select").equals("yes");
		}
	boolean isSparseCoding() {
		return this.hash.get("is_sparse_coding").equals("yes");
		}
	String getColumnName() {
		return getEntity() + "." + getName();
		}
	@Override
	public int hashCode() {
		return getColumnName().hashCode();
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof DataDict)) return true;
		return DataDict.class.cast(obj).getColumnName().equals(this.getColumnName());
		}
	}

private UKBiobankDataSet() {
	}

/** return dataset by column of null if not found */
public DataDict getDataDictByColumn(final String col) {
	return this.column2datadict.get(col);
	}

/** return codings or empty list  */
public List<Coding> findCodingsByName(final String coding_name) {
	return this.codings.getOrDefault(coding_name, Collections.emptyList());
	}

public static UKBiobankDataSet load(final String base) throws IOException {
	final UKBiobankDataSet ds= new UKBiobankDataSet();
	final CharSplitter splitter = CharSplitter.TAB;
	Path path = Paths.get(base+".dataset.codings.tsv");
	try(BufferedReader br=IOUtils.openPathForBufferedReading(path)) {
		String line = br.readLine();
		final FileHeader header = new FileHeader(splitter.split(line));
		while((line=br.readLine())!=null)  {
			final Coding coding = ds.new Coding(header.toMap(splitter.split(line)));
			List<Coding> L = ds.codings.get(coding.getCodingName());
			if(L==null) {
				L=new ArrayList<>();
				ds.codings.put(coding.getCodingName(),L);
				}
			L.add(coding);
			}
		}
	
	path = Paths.get(base+".data_dictionary.tsv");
	try(BufferedReader br=IOUtils.openPathForBufferedReading(path)) {
		String line = br.readLine();
		final FileHeader header = new FileHeader(splitter.split(line));
		while((line=br.readLine())!=null)  {
			final DataDict datadict = ds.new DataDict(header.toMap(splitter.split(line)));
			ds.column2datadict.put(datadict.getColumnName(),datadict);
			}
		}
	return ds;
	}
}
