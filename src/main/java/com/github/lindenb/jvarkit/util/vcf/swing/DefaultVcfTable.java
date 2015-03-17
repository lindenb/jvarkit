package com.github.lindenb.jvarkit.util.vcf.swing;

import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import java.util.Collections;
import java.util.List;


public class DefaultVcfTable extends AbstractVcfTable
	{
	private static final long serialVersionUID = 1L;
	private VCFHeader header=null;
	private AbstractVCFCodec codec=null;
	
	public DefaultVcfTable()
		{
		this(null, null);
		}
	
	public DefaultVcfTable(VCFHeader header,AbstractVCFCodec codec)
		{
		this.header=header;
		this.codec=codec;
		}
	
	public void reset(VCFHeader header,AbstractVCFCodec codec)
		{
		super.rows.clear();
		this.header=header;
		this.codec=codec;
		fireTableStructureChanged();
		}
	
	public VCFHeader getVCFHeader()
		{
		return header;
		}
	
	@Override
	public List<String> getSamples()
		{
		if(getVCFHeader()==null) return Collections.emptyList();
		return getVCFHeader().getSampleNamesInOrder();
		}
	
	@Override
	public List<String> getRows()
		{
		if(getVCFHeader()==null || getCodec()==null) return Collections.emptyList();
		return super.getRows();
		}
	
	@Override
	public AbstractVCFCodec getCodec() {
		return this.codec;
		}

	}
