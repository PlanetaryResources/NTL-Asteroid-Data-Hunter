/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import java.util.List;

/**
 * Represents search result.
 *
 * JSR303 validation is needed for this class.
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @param <T> the type of the entity in search result.
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class SearchResult<T> extends BaseSearchParameters {
    /**
     * Represents the total record count.
     *
     * Non-negative.
     */
    private int total;

    /**
     * Represents the total pages count.
     *
     * Non-negative.
     */
    private int totalPages;

    /**
     * Represents the values.
     *
     * Non-null.
     */
    private List<T> values;

    /**
     * <p>
     * The default constructor for class <code>SearchResult</code>.
     * </p>
     */
    public SearchResult() {
        // does nothing
    }

    /**
     * <p>
     * Gets the total record count.
     * </p>
     * @return the total record count.
     */
    public int getTotal() {
        return total;
    }

    /**
     * <p>
     * Sets the total record count.
     * </p>
     * @param total the total record count.
     */
    public void setTotal(int total) {
        this.total = total;
    }

    /**
     * <p>
     * Gets the total pages count.
     * </p>
     * @return the total pages count.
     */
    public int getTotalPages() {
        return totalPages;
    }

    /**
     * <p>
     * Sets the total pages count.
     * </p>
     * @param totalPages the total pages count.
     */
    public void setTotalPages(int totalPages) {
        this.totalPages = totalPages;
    }

    /**
     * <p>
     * Gets the values.
     * </p>
     * @return the values.
     */
    public List<T> getValues() {
        return values;
    }

    /**
     * <p>
     * Sets the values.
     * </p>
     * @param values the values.
     */
    public void setValues(List<T> values) {
        this.values = values;
    }
}

