/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services;

import gov.nasa.asteroid.hunter.models.NEO;
import gov.nasa.asteroid.hunter.models.NEOSearchCriteria;
import gov.nasa.asteroid.hunter.models.SearchResult;

/**
 * <p>
 * This service provides a method to search known NEOs.
 * </p>
 * <p>
 * <b>Thread Safety:</b> Implementations must be effectively thread safe. Refer
 * to ADS 1.3.4 for general assumptions and notes on thread safety.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public interface NEOService {

    /**
     * <p>
     * Searches the known NEOs.
     * </p>
     * @param criteria the search criteria
     *
     * @return the search result
     *
     * @throws IllegalArgumentException if criteria is null or invalid
     * @throws ServiceException if any other error occurred during the operation
     *
     */
    public SearchResult<NEO> search(NEOSearchCriteria criteria) throws ServiceException;
}

