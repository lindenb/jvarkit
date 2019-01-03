/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.PrintStream;

/**
 * used for knime support: each thread must have its own class to print to stdout
 * @author lindenb
 * usage:
 * <pre>
    PrintStream stdout=System.out;
	ThreadPrintStream tps=new ThreadPrintStream(System.out);
	System.setOut(tps);
	PrintStream out=new PrintStream("/tmp/x/y");
	MyMain.main(new String[]{});
	out.flush();
	out.close();
	System.setOut(stdout);
 * </pre>
 *
 */
public class ThreadPrintStream 
	extends PrintStream
	{
	private final static ThreadLocal<PrintStream> THREAD_TO_STDOUT = new ThreadLocal<PrintStream>();
	private PrintStream originalprintStream=null;
	public ThreadPrintStream(PrintStream originalprintStream)
		{
		super(originalprintStream);
		this.originalprintStream=originalprintStream;
		}
	
	public void set(PrintStream out)
		{
		if(THREAD_TO_STDOUT.get()!=null)
			 throw new IllegalStateException("PrintStream already opened for this thead");
		THREAD_TO_STDOUT.set(out);
		}
	
	public void dispose()
		{
		THREAD_TO_STDOUT.remove();
		}
	
	private PrintStream stream()
		{
		PrintStream out = THREAD_TO_STDOUT.get();
		return out==null?originalprintStream:out;
		}
	
	@Override
	public void println() {
		stream().println();
		}
	@Override
	public void flush() {
		stream().flush();
		}
	@Override
	public void print(String s) {
		stream().print(s);
		}
	
	@Override
	public void write(byte[] b) throws IOException {
		stream().write(b);
		}
	
	@Override
	public void write(byte[] buf, int off, int len) {
		stream().write(buf, off, len);
		}
	
	@Override
	public void write(int b) {
		stream().write(b);
		}
	
	@Override
	public boolean checkError() {
		return stream().checkError();
		}
	@Override
	public PrintStream append(char c) {
		return stream().append(c);
		}
	
	@Override
	public PrintStream append(final CharSequence csq) {
		return stream().append(csq);
		}
	@Override
	public PrintStream append(final CharSequence csq, int start, int end) {
		return stream().append(csq, start, end);
		}
	
	
	}
