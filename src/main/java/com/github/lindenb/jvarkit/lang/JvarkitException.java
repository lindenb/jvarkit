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

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.w3c.dom.Node;

import htsjdk.samtools.SAMSequenceDictionary;

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
	public FastaDictionaryMissing(final File file) {
		this(file.getPath());
		}
	public FastaDictionaryMissing(final String file) {
		super("A Sequence dictionary is missing for "+file+". A reference should have an associated .dict and .fai file. See http://gatkforums.broadinstitute.org/gatk/discussion/1601 .");
		}
	}

public static class VcfDictionaryMissing extends DictionaryMissing
{
public VcfDictionaryMissing(final File file) {
	this(file.getPath());
	}
public VcfDictionaryMissing(final String file) {
	super("A Sequence dictionary is missing for "+file+". A VCF should have a set of `##contig` in its header. See also https://broadinstitute.github.io/picard/command-line-overview.html#UpdateVcfSequenceDictionary ");
	}
}


public static class SampleMissing extends DictionaryMissing
{
public SampleMissing(final String file) {
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
/** exception thrown when we cannot convert a contig */
public static class ContigNotFound extends Error
	{
	public ContigNotFound(final String msg) {
		super(msg);
		}
	}
/** exception thrown when we cannot find a contig in a dict */
public static class ContigNotFoundInDictionary extends ContigNotFound
	{
	public ContigNotFoundInDictionary(final String contig,final SAMSequenceDictionary dict) {
		super("Cannot find contig \""+contig+"\" in dictionary:["+dict.getSequences().stream().map(SSR->SSR.getSequenceName()).collect(Collectors.joining(","))+"]");
		}
	}

/** exception two dicts are not the same */
public static class DictionariesAreNotTheSame extends Error
	{
	public DictionariesAreNotTheSame(final SAMSequenceDictionary dict1,final SAMSequenceDictionary dict2) {
		super("Two dictionaries are not the same:"+
				dict1.getSequences().stream().map(SSR->SSR.getSequenceName()).collect(Collectors.joining(","))
				+"\nand\n"+
				dict2.getSequences().stream().map(SSR->SSR.getSequenceName()).collect(Collectors.joining(","))
				);
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

public static class BerkeleyDbError  extends Error
	{	
	public BerkeleyDbError(final String msg) {
		super(msg);
		}
	}

public static class XmlDomError  extends Error
	{	
	private static String node2str(final Node node) {
		return "xml";
		}
	public XmlDomError(final String msg) {
		super(msg);
		}
	public XmlDomError(final Node node,final String msg) {
		this(node2str(node)+" : "+ msg);
		}
	}
public static class ScriptingError  extends Error
	{	
	public ScriptingError(final String error) {
		super(error);
		}
	public ScriptingError(final Exception error) {
		super(error);
		}
	}
public static class ScriptEngineNotFound  extends ScriptingError
	{	
	public ScriptEngineNotFound(final String engine) {
		super("ScriptEngineManager for \""+engine+"\". Do you use the SUN/Oracle JDK ?");
		}
	}
public static class JavaScriptEngineNotFound  extends ScriptEngineNotFound
	{
	public JavaScriptEngineNotFound() {
		super("javascript");
		}
	}


public static class FileFormatError  extends Error
	{	
	public FileFormatError(final String msg) {
		super(msg);
		}
	}
public static class TokenErrors  extends FileFormatError
	{	
	private static String vertical(final List<String> tokens) {
		if(tokens==null) return "(tokens:null)";
		StringBuilder sb=new StringBuilder("\n");
		for(int i=0;i<tokens.size();++i)
				sb.append(" * (").append(String.valueOf(i+1)).append(")\t\"").append(tokens.get(i)).append("\"\n");
		return sb.toString();
		}
	public TokenErrors(String tokens[]) {
		this(Arrays.asList(tokens));
		}

	public TokenErrors(List<String> tokens) {
		super(vertical(tokens));
		}

	public TokenErrors(final String msg,List<String> tokens) {
		super(msg+vertical(tokens));
		}
	public TokenErrors(final String msg,String tokens[]) {
		this(msg,Arrays.asList(tokens));
		}
	}

}
