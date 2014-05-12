package com.github.lindenb.jvarkit.tools.ws.rmi.server;

import java.io.File;
import java.io.IOException;
import java.rmi.RMISecurityManager;
import java.rmi.RemoteException;
import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;
import java.rmi.server.UnicastRemoteObject;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import htsjdk.samtools.cmdline.Option;
import htsjdk.samtools.cmdline.Usage;
import htsjdk.samtools.util.Log;

import com.github.lindenb.jvarkit.tools.ws.WSBam;
import com.github.lindenb.jvarkit.tools.ws.WSProject;
import com.github.lindenb.jvarkit.tools.ws.WSVcf;
import com.github.lindenb.jvarkit.tools.ws.impl.WSBamImpl;
import com.github.lindenb.jvarkit.tools.ws.impl.WSProjectImpl;
import com.github.lindenb.jvarkit.tools.ws.impl.WSReferenceImpl;
import com.github.lindenb.jvarkit.tools.ws.impl.WSVcfImpl;
import com.github.lindenb.jvarkit.tools.ws.rmi.NGSService;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

public class NGSServiceImpl extends AbstractCommandLineProgram
	implements NGSService
	{
	private static final Log LOG = Log.getInstance(NGSServiceImpl.class);
    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " RMI Server for NGSService.";

    @Option(shortName="CFG", doc="XML Config.",optional=false)
    public File CONFIG=null;

	
	
	private Map<String,WSReferenceImpl> id2reference=new HashMap<String,WSReferenceImpl>();
	private Map<String,WSProjectImpl> id2project=new HashMap<String,WSProjectImpl>();
	private Map<String,WSBamImpl> id2bam=new HashMap<String,WSBamImpl>();
	private Map<String,WSVcfImpl> id2vcf=new HashMap<String,WSVcfImpl>();
	
	@Override
	public List<? extends WSProject> getProjects()
		{
		return new ArrayList<WSProject>(id2project.values());
		}
	
	@Override
	public WSProject getProjectById(String id)
			throws RemoteException
		{
		return this.id2project.get(id);
		}
	
	@Override
	public WSBam getBamById(String id) throws RemoteException {
		return id2bam.get(id);
	}
	
	@Override
	public WSVcf getVcfById(String id) throws RemoteException {
		return id2vcf.get(id);
	}
	
	
	public String getBamHeader(String bamId)
		{
		return null;//TODO
		}
	
	
	private void loadConfig(File xmlFile) throws IOException,XMLStreamException,JAXBException
		{
		LOG.info("reading config "+xmlFile);
		JAXBContext ctx=JAXBContext.newInstance(WSReferenceImpl.class);
		Unmarshaller u=ctx.createUnmarshaller();
		XMLInputFactory xif=XMLInputFactory.newFactory();
		XMLEventReader r=xif.createXMLEventReader(new StreamSource(xmlFile));
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
			if(!evt.isStartElement() ) { r.next(); continue;}
			String localName=evt.asStartElement().getName().getLocalPart();
			if(localName.equals(WSReferenceImpl.xmlName))
				{
				WSReferenceImpl ref=u.unmarshal(r,WSReferenceImpl.class).getValue();
				id2reference.put(ref.getId(), ref);
				}
			else if(localName.equals(WSProjectImpl.xmlName))
				{
				WSProjectImpl p=u.unmarshal(r,WSProjectImpl.class).getValue();
				id2project.put(p.getId(), p);
				}
			else
				{
				r.next();
				}
			}
		r.close();
		}
	
	@Override
	protected int doWork()
		{
		try
			{
			if (System.getSecurityManager() == null)
				{
				System.setSecurityManager(new RMISecurityManager()
				);
				}
			NGSServiceImpl engine = new NGSServiceImpl();
			engine.loadConfig(CONFIG);
			/* exports the supplied remote object so that it can
			receive invocations of its remote methods from remote clients. */
			NGSService stub = (NGSService) UnicastRemoteObject.exportObject(
								engine,
								0 //TCP port
								);
			/* adds the name to the RMI registry running on the server */
			Registry registry = LocateRegistry.getRegistry();
			registry.rebind(NGSService.SERVICE_NAME, stub);
			LOG.info("Service bound.");
			return 0;
			}
		catch(SecurityException err)
			{
			LOG.error(err, "Error occured");
			return -1;
			}
		catch(Exception err)
			{
			err.printStackTrace();
			return -1;
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new NGSServiceImpl().instanceMain(args);
		}

}
