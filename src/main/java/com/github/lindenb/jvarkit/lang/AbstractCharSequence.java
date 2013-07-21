package com.github.lindenb.jvarkit.lang;

public abstract class AbstractCharSequence
	implements CharSequence
	{
	@Override
	public int hashCode()
		{
		return getString().hashCode();
		}
	
	public String getString()
		{
		StringBuilder b=new StringBuilder(length());
		for(int i=0;i< length();++i) b.append(charAt(i));
		return b.toString();
		}
	
	@Override
	public String toString()
		{
		return getString();
		}
	@Override
	public CharSequence subSequence(int start, int end)
		{
		return new SubSequence(this,start,end);
		}
	}
