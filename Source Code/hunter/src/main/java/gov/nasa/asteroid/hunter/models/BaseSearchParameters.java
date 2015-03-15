/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import javax.validation.constraints.Min;

/**
 * <p>
 * Represents base search parameters.
 * </p>
 *
 * <p>
 * JSR303 validation is needed for this class.
 * </p>
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public abstract class BaseSearchParameters {

    /**
     * Represents the page size.
     *
     * Optional. Positive.
     */
    @Min(1)
    private int pageSize;

    /**
     * Represents the page number.
     *
     * Positive integer or 0 (no paging).
     *
     * Required.
     */
    @Min(0)
    private int pageNumber;

    /**
     * Represents the sort by field.
     *
     * Optional. Valid field value of the entity to search.
     */
    private String sortBy;

    /**
     * Represents the sort type.
     *
     * Optional.
     */
    private SortType sortType;

    /**
     * <p>
     * The default constructor for class <code>BaseSearchParameters</code>.
     * </p>
     */
    protected BaseSearchParameters() {
        // does nothing
    }

    /**
     * <p>
     * Gets the page size.
     * </>
     * @return the page size.
     */
    public int getPageSize() {
        return pageSize;
    }

    /**
     * <p>
     * Sets the page size.
     * </p>
     * @param pageSize the page size.
     */
    public void setPageSize(int pageSize) {
        this.pageSize = pageSize;
    }

    /**
     * <p>
     * Gets the page number.
     * </>
     * @return the page number.
     */
    public int getPageNumber() {
        return pageNumber;
    }

    /**
     * <p>
     * Sets the page number.
     * </p>
     * @param pageNumber the page number.
     */
    public void setPageNumber(int pageNumber) {
        this.pageNumber = pageNumber;
    }

    /**
     * <p>
     * Gets the sort by field.
     * </>
     * @return the sort by field.
     */
    public String getSortBy() {
        return sortBy;
    }

    /**
     * <p>
     * Sets the sort by field.
     * </p>
     * @param sortBy the sort by field.
     */
    public void setSortBy(String sortBy) {
        this.sortBy = sortBy;
    }

    /**
     * <p>
     * Gets the sort by field.
     * </>
     * @return the sort by field.
     */
    public SortType getSortType() {
        return sortType;
    }

    /**
     * <p>
     * Sets the sort type.
     * </p>
     * @param sortType the sort type.
     */
    public void setSortType(SortType sortType) {
        this.sortType = sortType;
    }
}

