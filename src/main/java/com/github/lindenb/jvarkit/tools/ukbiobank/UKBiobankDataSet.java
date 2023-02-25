package com.github.lindenb.jvarkit.tools.ukbiobank;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

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
	}

private UKBiobankDataSet() {
	}

public DataDict getDataDictByColumn(final String col) {
	return this.column2datadict.get(col);
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
