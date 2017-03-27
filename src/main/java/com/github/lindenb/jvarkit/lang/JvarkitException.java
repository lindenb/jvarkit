/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.lang;
@SuppressWarnings("serial")
public class JvarkitException   {

public static class ReferenceMissing extends Error
	{
	public ReferenceMissing(final String msg) {
		super("Reference missing : "+msg);
		}
	}
	
public static class DictionaryMissing extends Error
	{
	public DictionaryMissing(final String msg) {
		super("Dictionary missing : "+msg);
		}
	}
public static class FastaDictionaryMissing extends DictionaryMissing
{
public FastaDictionaryMissing(final String file) {
	super(file);
	}
}

	
/** exception thrown when the user do something wrong */
public static class UserError extends Error
	{
	public UserError(final String msg) {
		super("User Error : "+msg);
		}
	}

	
/** exception thrown when the user made an error on the command line */
public static class CommandLineError extends Error
	{
	public CommandLineError(final String msg) {
		super("There was an error in the command line arguments. Please check your parameters : "+msg);
		}
	}

/** programming error */
public static class ProgrammingError extends Error
	{
	public ProgrammingError(final String msg) {
		super("Hum... There is something wrong, must be a programming error...: "+msg);
		}
	}
/** should never happen */
public static class ShouldNeverHappen extends Error
	{
	public ShouldNeverHappen(final String msg) {
		super("Hum... There is something wrong. This should never happen."+msg);
		}
	}

/** should never happen */
public static class MixingApplesAndOranges extends ShouldNeverHappen
	{
	public MixingApplesAndOranges(final Object a,final Object b) {
		this((a==null?"null":a.getClass().toString())+
				" and "+
			 (b==null?"null":b.getClass().toString()));
		}
	public MixingApplesAndOranges(final String msg) {
		super("Mixing apples and oranges."+msg);
		}
	}

}
