package com.github.lindenb.jvarkit.lang;

public class SubSequence extends DelegateCharSequence
	{
	private int start;
	private int end;
	public SubSequence(CharSequence delegate,int start, int end)
		{
		super(delegate);
		this.start=start;
		this.end=end;
		}
	@Override
	public int length()
		{
		return end-start;
		}
	@Override
	public char charAt(int index)
		{
		return getDelegate().charAt(start+index);
		}
	}
