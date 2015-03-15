/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import javax.persistence.Basic;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Table;


/**
 * Represents a help item.
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
@Table(name = "help_item")
public class HelpItem extends IdentifiableEntity {

    /**
     * Represents the topic.
     *
     * Required.
     */
    @ManyToOne
    @JoinColumn(name = "help_topic_id", nullable = false)
    private HelpTopic topic;

    /**
     * Represents the title.
     *
     * Required.
     */
    @Basic
    @Column(nullable = false)
    private String title;

    /**
     * Represents the content.
     *
     * Required.
     */
    @Basic
    @Column(nullable = false)
    private String content;

    /**
     * <p>
     * The default constructor for class <code>HelpItem</code>.
     * </p>
     */
    public HelpItem() {
        // does nothing
    }

    /**
     * <p>
     * Gets the topic.
     * </p>
     * @return the topic.
     */
    public HelpTopic getTopic() {
        return topic;
    }

    /**
     * <p>
     * Sets the topic.
     * </p>
     * @param topic the topic.
     */
    public void setTopic(HelpTopic topic) {
        this.topic = topic;
    }

    /**
     * <p>
     * Gets the title.
     * </p>
     * @return the title.
     */
    public String getTitle() {
        return title;
    }

    /**
     * <p>
     * Sets the title.
     * </p>
     * @param title the title.
     */
    public void setTitle(String title) {
        this.title = title;
    }

    /**
     * <p>
     * Gets the content.
     * </p>
     * @return the content.
     */
    public String getContent() {
        return content;
    }

    /**
     * <p>
     * Sets the content.
     * </p>
     * @param content the content.
     */
    public void setContent(String content) {
        this.content = content;
    }
}

