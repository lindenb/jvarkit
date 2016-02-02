package com.github.lindenb.jvarkit.tools.server;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.handler.DefaultHandler;

public  class ProjectServer extends AbstractProjectServer {

	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(ProjectServer.class);
	
	private static class ProjectHandler extends DefaultHandler{
		final List<String> args;
		ProjectHandler(final List<String> args) {
			this.args=args;
		}
		@Override
		public void handle(String target, Request baseRequest, HttpServletRequest request, HttpServletResponse response)
				throws IOException, ServletException {
			 response.setContentType("text/plain; charset=utf-8");
			 PrintWriter pw=response.getWriter();
			 for(String s:this.args)
			 	{
				 pw.println(s);
			 	}
			 baseRequest.setHandled(true);
		}
	}
	
	@Override
	protected Handler createDefaultHandler(List<String> args) {
		return new ProjectHandler(args);
	}
	
	@Override
	public java.util.Collection<Throwable> call() throws Exception
	{
	try {
		super.createAndRunServer(super.getInputFiles());
	    return RETURN_OK;
		}
	catch (Exception e) {
		return wrapException(e);
		}
	
	}	
	
public static void main(String[] args) throws Exception{
    new ProjectServer().instanceMainWithExit(args);
	}

}
