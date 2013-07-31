package com.github.lindenb.jvarkit.tools.ws.rmi;

import java.rmi.Remote;
import java.rmi.RemoteException;
import java.util.List;

import com.github.lindenb.jvarkit.tools.ws.WSBam;
import com.github.lindenb.jvarkit.tools.ws.WSProject;
import com.github.lindenb.jvarkit.tools.ws.WSVcf;

public interface NGSService extends Remote
	{
	public static final String SERVICE_NAME="ngsservice";
	public  List<? extends WSProject> getProjects() throws RemoteException;
	public  WSProject getProjectById(String id) throws RemoteException;
	public  WSBam getBamById(String id) throws RemoteException;
	public  WSVcf getVcfById(String id) throws RemoteException;
	}
