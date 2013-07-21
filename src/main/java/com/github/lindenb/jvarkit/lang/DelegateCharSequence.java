package com.github.lindenb.jvarkit.lang;

public class DelegateCharSequence
	extends AbstractCharSequence
	{
	private CharSequence delegate;
	public DelegateCharSequence(CharSequence delegate)
		{
		this.delegate=delegate;
		}
	public CharSequence getDelegate()
		{
		return delegate;
		}
	@Override
	public char charAt(int index)
		{
		return getDelegate().charAt(index);
		}
	@Override
	public int length()
		{
		return getDelegate().length();
		}
	
	}
