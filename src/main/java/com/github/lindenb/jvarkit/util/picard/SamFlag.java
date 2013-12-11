package com.github.lindenb.jvarkit.util.picard;

import java.util.HashSet;
import java.util.Set;

public enum SamFlag
	{
		READ_PAIRED(0x1),
		READ_MAPPED_IN_PROPER_PAIR(0x2),
		READ_UNMAPPED(0x4),
		MATE_UNMAPPED(0x8),
		READ_REVERSE_STRAND(0x10),
		MATE_REVERSE_STRAND(0x20),
		FIRST_IN_PAIR(0x40),
		SECOND_IN_PAIR(0x80),
		NOT_PRIMARY_ALIGNMENT(0x100),
		READ_FAILS_VENDOR_QUALITY_CHECK(0x200),
		READ_IS_DUPLICATE(0x400),
		SUPPLEMENTARY_ALIGNMENT(0x800)
		;
	private int flag;
	SamFlag(int flag)
		{
		this.flag=flag;
		}
	
	public int getFlag()
		{
		return flag;
		}

	public String getLabel()
		{
		return name().toLowerCase().replace('_', ' ');
		}
	
	public static SamFlag valueOf(int flg)
		{
		for(SamFlag f:values())
			{
			if(flg==f.getFlag()) return f;
			}
		return null;
		}
	
	public boolean isSet(int flag)
		{
		return (getFlag() & flag)!=0; 
		}
	
	public boolean isUnset(int flag)
		{
		return !isSet(flag);
		}
	
	public static Set<SamFlag> getFlags(int flag)
		{
		Set<SamFlag> set=new HashSet<SamFlag>();
		for(SamFlag f:values())
			{
			if(f.isSet(flag)) set.add(f);
			}
		return set;
		}
	
	}
