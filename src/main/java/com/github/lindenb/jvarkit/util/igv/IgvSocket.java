package com.github.lindenb.jvarkit.util.igv;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Controlling IGV through a Port
 * 
 * http://www.broadinstitute.org/software/igv/PortCommands
 * @author lindenb
 *
 */
public class IgvSocket
	implements Closeable
	{
	public static final int DEFAULT_PORT=60151;
	public static final String DEFAULT_HOST="127.0.0.1";
	
	private String host=DEFAULT_HOST;
	private int port=DEFAULT_PORT;
	private Socket _socket=null;
	private PrintWriter _out=null;
	private BufferedReader _in=null;
	
	public IgvSocket()
		{
		
		}
	private IgvSocket(final IgvSocket src)
		{
		this.host = src.host;
		this.port = src.port;
		}
	
	public String getHost() {
		return host;
		}
	public int getPort() {
		return port;
		}
	public void setHost(final String host) {
		this.host = host;
		}
	public void setPort(int port) {
		this.port = port;
		}
	
	public Socket getSocket() throws IOException
		{
		if(_socket==null)
			{
			_socket=new Socket(getHost(), getPort());
			}
		return _socket;
		}
	
	public PrintWriter getWriter() throws IOException
		{
		if(_out==null)
			{
			_out= new PrintWriter(getSocket().getOutputStream(), true);
			}
		return _out;
		}

	public BufferedReader getReader() throws IOException
		{
		if(_in==null)
			{
			_in=  new BufferedReader(new InputStreamReader(getSocket().getInputStream()));
			}
		return _in;
		}
	
	
	@Override
	public void close() {
		CloserUtil.close(_in);
		_in=null;
		CloserUtil.close(_out);
		_out=null;
		CloserUtil.close(_socket);
		_socket=null;
		}
	
	public void show(final VariantContext ctx)
		{
		if( ctx ==null ) return ;
		show(ctx.getContig(),ctx.getStart());
		}	
	
	public void show(final String chrom,int chromStart)
		{
		if(chrom==null || chrom.trim().isEmpty() || chromStart<0) return;
		execute(Collections.singletonList("goto "+chrom+":"+chromStart));
		}
	
	public Runnable buildRunnable(final List<String> cmds) 
		{
		return new ShowRunner(this,cmds);
		}

	
	
	public void execute(final List<String> cmds) 
		{
		new Thread(buildRunnable(cmds)).start();
		}
	
	
	private static class ShowRunner implements Runnable
		{
		final IgvSocket copy;
		//String chrom;
		//int start;
		final List<String> commands;
		
		ShowRunner(final IgvSocket socket,final List<String> command) {
			this.commands = new ArrayList<>(command);
			this.copy =  new IgvSocket(socket);
			}
		@Override
		public void run() {

			PrintWriter out=null;
			try
				{
				out = this.copy.getWriter();
				this.copy.getReader();
				for(final String command:this.commands)
					{
					out.println(command);
					try{ Thread.sleep(5*1000);}
					catch(InterruptedException err2) {}
					 }
				}
			catch(final Exception err)
				{
				err.printStackTrace();
				}
			finally
				{
				CloserUtil.close(this.copy);
				}			
			}
		};
	

	}
