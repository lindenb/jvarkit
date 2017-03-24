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
package com.github.lindenb.jvarkit.util.command;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.lang.JvarkitException;


public abstract class AbstractArgument<T,BUILDER extends AbstractArgument<T,BUILDER>>
	implements Argument<T>
	{
	private static int ID_GENERATOR=0;
	private String id="id"+(++ID_GENERATOR);
	protected String argName=null;
	protected String shortArg=null;
	protected String longArg=null;
	protected String shortDesc=null;
	protected String longDesc=null;
	protected String deprecatedMsg=null;
	/** this command should be hidden from the help menu */
	protected boolean hidden=false;
	protected final Set<String> categories=new HashSet<>();
	protected final Map<String,String> metaData=new LinkedHashMap<>();
	private T value=null;
	
	protected abstract BUILDER getThis();	
	public BUILDER shortOpt(final char c) { return shortOpt(String.valueOf(c));} 
	public BUILDER shortOpt(final String s) { this.shortArg = s; return getThis();} 
	public BUILDER longOpt(final String s) { this.longArg = s; return getThis();} 
	public BUILDER shortDesc(final String s) { this.shortDesc = s; return getThis();} 
	public BUILDER longDesc(final String s) { this.longDesc = s; return getThis();} 
	public BUILDER validate() {  return getThis();} 
	public BUILDER deprecated(final String msg) { this.deprecatedMsg=msg; return getThis();} 
	/** shortcut for marking the option as deprecated with a default message */
	public BUILDER deprecated() { return deprecated("This option is deprecated");}
	/** if true, the option will not be listed in the help message */
	public BUILDER hidden(boolean b) { this.hidden=b; return getThis();} 
	/** shortcut for hidden(false) */
	public BUILDER hidden() { return hidden(true);} 		
	public BUILDER categories(final String...s) { this.categories.addAll(Arrays.asList(s)); return getThis();} 
	public BUILDER category(final String s) { this.categories.add(s); return getThis();} 
	
	protected AbstractArgument() {
		this(null);
	}
	protected AbstractArgument(final T value) {
		this.value=null;
	}

	@Override
	public String getLongDesc() {
		return this.longDesc;
		}
	
	@Override
	public Map<String, String> getMetaData() {
		return this.metaData;
		}
	
	public String getDeprecationMsg() {
		return deprecatedMsg;
		}
	@Override
	public boolean isDeprecated() {
		return this.getDeprecationMsg()!=null;
		}
	@Override
	public boolean isHidden() {
		return hidden;
		}
	@Override
	public String getShortOpt() {
		return shortArg;
		}
	@Override
	public String getLongOpt() {
		return longArg;
		}
	
	/** get an unique identifier for this argument */
	protected String getId() {
		return id;
		}
	/** returns the hashcode the getId() */
	@Override
	public int hashCode() {
		return getId().hashCode();
		}
	
	/** two arguments are the same if they have the same getId() */
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj == null ) return false;
		if(!(obj instanceof AbstractArgument)) throw new JvarkitException.MixingApplesAndOranges(this,obj);
		return AbstractArgument.class.cast(obj).getId().equals(this.getId());
		}
	
	/** get the categories associated to this arguments.*/
	@Override
	public Set<String> getCategories() {
		return Collections.unmodifiableSet(this.categories);
		}
	
	/** returns a short description for this argument */
	public String getShortDesc() {
		return shortDesc;
		}
	
	public String getArgName() {
		return argName;
		}
	
	public String getLabel() {
		StringBuilder sb=new StringBuilder();
		if(hasShortOpt()) sb.append("-").append(getShortOpt());
		if(hasLongOpt() && hasShortOpt()) sb.append("|");
		if(hasLongOpt()) sb.append("--").append(getLongOpt());
		return sb.toString();
		}	
	
	protected static class ParsedLongArg
		{
		Map<String,String> metadata=Collections.emptyMap();
		String value=null;
		}
	
	protected Map<String,String> parseMetaData(final String s) {
		if(s==null || s.isEmpty()) return Collections.emptyMap();
		final Map<String,String> hash = new LinkedHashMap<>();
		for(final String str:s.split("[,]"))
			{
			final String key,value;
			int eq=str.indexOf("=");
			if(eq==-1) 
				{
				key=str;
				value="";
				}
			else
				{
				key=str.substring(0,eq);
				value=str.substring(eq+1);
				}
			if(hash.containsKey(key)) throw new JvarkitException.CommandLineError("already set "+key);
			hash.put(key, value);
			}
		return hash;
		}
	
	protected boolean isSupportingMetaData() 
	{
		return true;
	}
	
	protected boolean isSupportingEqualSignInArg() 
	{
		return true;
	}
	
	protected ParsedLongArg parseLongArg(final String optarg) {
		if(!hasLongOpt()) return null;
		
		if(!optarg.startsWith("--")) return null;
		if(optarg.equals("--")) return null;
		final ParsedLongArg pa= new ParsedLongArg();
		String s=optarg.substring(2);
		final int colon=s.indexOf(':');
		if(colon!=-1 && !isSupportingMetaData()) {
			throw new JvarkitException.CommandLineError(this.getLabel()+"Option doesn't support metadata "+optarg);
			}
		final int eq=s.indexOf('=');
		if(eq!=-1 && !isSupportingEqualSignInArg()) {
			throw new JvarkitException.CommandLineError(this.getLabel()+"Option doesn't support option "+optarg);
			}
		
		if(colon!=-1)
			{
			pa.metadata= parseMetaData(s.substring(colon+1,eq));
			pa.value = null;
			s=s.substring(0,colon);
			}
		else if(eq!=1)
			{
			pa.value = s.substring(eq+1);
			s=s.substring(0,eq);
			}
		if(!s.equals(getLongOpt())) return null;	
		return pa;
		}
	
	protected boolean isLongOpt( String s){
		return parseLongArg(s)!=null;
		}
	
	protected boolean isShortOpt( String s){
		if(!hasShortOpt()) return false;
		if(!s.startsWith("-")) return false;
		if(s.startsWith("--")) return false;
		s=s.substring(1);
		return s.equals(getShortOpt());
	}
	
	@Override
	public T get() {
		return this.value;
		}
	
	@Override
	public void set(final T t)
		{
		this.value=t;
		}	
	
	
	}
