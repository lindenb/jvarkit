//
// This file was generated by the Eclipse Implementation of JAXB, v4.0.5 
// See https://eclipse-ee4j.github.io/jaxb-ri 
// Any modifications to this file will be lost upon recompilation of the source schema. 
//


package org.uniprot;

import java.util.ArrayList;
import java.util.List;
import jakarta.xml.bind.annotation.XmlAccessType;
import jakarta.xml.bind.annotation.XmlAccessorType;
import jakarta.xml.bind.annotation.XmlAttribute;
import jakarta.xml.bind.annotation.XmlType;


/**
 * Describes non-nuclear gene locations (organelles and plasmids).
 *             Equivalent to the flat file OG-line.
 * 
 * <p>Java class for geneLocationType complex type</p>.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.</p>
 * 
 * <pre>{@code
 * <complexType name="geneLocationType">
 *   <complexContent>
 *     <restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       <sequence>
 *         <element name="name" type="{http://uniprot.org/uniprot}statusType" maxOccurs="unbounded" minOccurs="0"/>
 *       </sequence>
 *       <attribute name="type" use="required">
 *         <simpleType>
 *           <restriction base="{http://www.w3.org/2001/XMLSchema}string">
 *             <enumeration value="apicoplast"/>
 *             <enumeration value="chloroplast"/>
 *             <enumeration value="organellar chromatophore"/>
 *             <enumeration value="cyanelle"/>
 *             <enumeration value="hydrogenosome"/>
 *             <enumeration value="mitochondrion"/>
 *             <enumeration value="non-photosynthetic plastid"/>
 *             <enumeration value="nucleomorph"/>
 *             <enumeration value="plasmid"/>
 *             <enumeration value="plastid"/>
 *           </restriction>
 *         </simpleType>
 *       </attribute>
 *       <attribute name="evidence" type="{http://uniprot.org/uniprot}intListType" />
 *     </restriction>
 *   </complexContent>
 * </complexType>
 * }</pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "geneLocationType", propOrder = {
    "name"
})
public class GeneLocationType {

    protected List<StatusType> name;
    @XmlAttribute(name = "type", required = true)
    protected String type;
    @XmlAttribute(name = "evidence")
    protected List<Integer> evidence;

    /**
     * Gets the value of the name property.
     * 
     * <p>This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the name property.</p>
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * </p>
     * <pre>
     * getName().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link StatusType }
     * </p>
     * 
     * 
     * @return
     *     The value of the name property.
     */
    public List<StatusType> getName() {
        if (name == null) {
            name = new ArrayList<>();
        }
        return this.name;
    }

    /**
     * Gets the value of the type property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getType() {
        return type;
    }

    /**
     * Sets the value of the type property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setType(String value) {
        this.type = value;
    }

    /**
     * Gets the value of the evidence property.
     * 
     * <p>This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the evidence property.</p>
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * </p>
     * <pre>
     * getEvidence().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link Integer }
     * </p>
     * 
     * 
     * @return
     *     The value of the evidence property.
     */
    public List<Integer> getEvidence() {
        if (evidence == null) {
            evidence = new ArrayList<>();
        }
        return this.evidence;
    }

}