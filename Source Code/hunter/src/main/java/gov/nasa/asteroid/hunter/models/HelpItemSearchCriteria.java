/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import javax.validation.constraints.Min;


/**
 * Represents help item search criteria.
 *
 * JSR303 validation is needed for this class.
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class HelpItemSearchCriteria extends BaseSearchParameters {
    /**
     * Represents the topic ID.
     *
     * Optional. Positive.
     */
    @Min(1)
    private Long topicId;

    /**
     * Represents the keyword.
     *
     * Optional.
     */
    private String keyword;

    /**
     * <p>
     * The default constructor for class <code>HelpItemSearchCriteria</code>.
     * </p>
     */
    public HelpItemSearchCriteria() {
        // does nothing
    }

    /**
     * <p>
     * Gets the topic ID.
     * </p>
     * @return the topic ID.
     */
    public Long getTopicId() {
        return topicId;
    }

    /**
     * <p>
     * Sets topic ID.
     * </p>
     * @param topicId the topic ID.
     */
    public void setTopicId(Long topicId) {
        this.topicId = topicId;
    }

    /**
     * <p>
     * Gets the keyword.
     * </p>
     * @return the keyword.
     */
    public String getKeyword() {
        return keyword;
    }

    /**
     * <p>
     * Sets the keyword.
     * </p>
     * @param keyword the keyword.
     */
    public void setKeyword(String keyword) {
        this.keyword = keyword;
    }
}

