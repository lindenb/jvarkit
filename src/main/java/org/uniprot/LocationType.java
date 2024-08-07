//
// This file was generated by the Eclipse Implementation of JAXB, v4.0.5 
// See https://eclipse-ee4j.github.io/jaxb-ri 
// Any modifications to this file will be lost upon recompilation of the source schema. 
//


package org.uniprot;

import jakarta.xml.bind.annotation.XmlAccessType;
import jakarta.xml.bind.annotation.XmlAccessorType;
import jakarta.xml.bind.annotation.XmlAttribute;
import jakarta.xml.bind.annotation.XmlType;


/**
 * Describes a sequence location as either a range with a begin and end or as a position. The 'sequence' attribute is only used when the location is not on the canonical sequence displayed in the current entry.
 * 
 * <p>Java class for locationType complex type</p>.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.</p>
 * 
 * <pre>{@code
 * <complexType name="locationType">
 *   <complexContent>
 *     <restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       <choice>
 *         <sequence>
 *           <element name="begin" type="{http://uniprot.org/uniprot}positionType"/>
 *           <element name="end" type="{http://uniprot.org/uniprot}positionType"/>
 *         </sequence>
 *         <element name="position" type="{http://uniprot.org/uniprot}positionType"/>
 *       </choice>
 *       <attribute name="sequence" type="{http://www.w3.org/2001/XMLSchema}string" />
 *     </restriction>
 *   </complexContent>
 * </complexType>
 * }</pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "locationType", propOrder = {
    "begin",
    "end",
    "position"
})
public class LocationType {

    protected PositionType begin;
    protected PositionType end;
    protected PositionType position;
    @XmlAttribute(name = "sequence")
    protected String sequence;

    /**
     * Gets the value of the begin property.
     * 
     * @return
     *     possible object is
     *     {@link PositionType }
     *     
     */
    public PositionType getBegin() {
        return begin;
    }

    /**
     * Sets the value of the begin property.
     * 
     * @param value
     *     allowed object is
     *     {@link PositionType }
     *     
     */
    public void setBegin(PositionType value) {
        this.begin = value;
    }

    /**
     * Gets the value of the end property.
     * 
     * @return
     *     possible object is
     *     {@link PositionType }
     *     
     */
    public PositionType getEnd() {
        return end;
    }

    /**
     * Sets the value of the end property.
     * 
     * @param value
     *     allowed object is
     *     {@link PositionType }
     *     
     */
    public void setEnd(PositionType value) {
        this.end = value;
    }

    /**
     * Gets the value of the position property.
     * 
     * @return
     *     possible object is
     *     {@link PositionType }
     *     
     */
    public PositionType getPosition() {
        return position;
    }

    /**
     * Sets the value of the position property.
     * 
     * @param value
     *     allowed object is
     *     {@link PositionType }
     *     
     */
    public void setPosition(PositionType value) {
        this.position = value;
    }

    /**
     * Gets the value of the sequence property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * Sets the value of the sequence property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setSequence(String value) {
        this.sequence = value;
    }

}
