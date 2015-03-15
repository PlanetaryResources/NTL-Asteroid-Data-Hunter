/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import javax.persistence.Basic;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Table;


/**
 * Represents a help topic.
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
@Entity
@Table(name = "help_topic")
public class HelpTopic extends IdentifiableEntity {

    /**
     * Represents the topic name.
     *
     * Required.
     */
    @Basic
    @Column(nullable = false)
    private String name;

    /**
     * <p>
     * The default constructor for class <code>HelpTopic</code>.
     * </p>
     */
    public HelpTopic() {
        // does nothing
    }

    /**
     * <p>
     * Gets the topic name.
     * </p>
     * @return the topic name.
     */
    public String getName() {
        return name;
    }

    /**
     * <p>
     * Sets the topic name.
     * </p>
     * @param name the topic name.
     */
    public void setName(String name) {
        this.name = name;
    }
}

