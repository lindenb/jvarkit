package com.github.lindenb.jvarkit.util.igv;

import htsjdk.samtools.util.CloserUtil;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;

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
	
	public String getHost() {
		return host;
		}
	public int getPort() {
		return port;
		}
	public void setHost(String host) {
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
	}
