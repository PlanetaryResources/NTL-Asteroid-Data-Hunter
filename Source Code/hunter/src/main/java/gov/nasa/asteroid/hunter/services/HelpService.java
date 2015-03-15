/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services;

import gov.nasa.asteroid.hunter.models.HelpItem;
import gov.nasa.asteroid.hunter.models.HelpItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.HelpTopic;
import gov.nasa.asteroid.hunter.models.SearchResult;

import java.util.List;

/**
 * <p>
 * This service provides methods to access help contents.
 * </p>
 * <p>
 * <b>Thread Safety:</b> Implementations must be effectively thread safe. Refer
 * to ADS 1.3.4 for general assumptions and notes on thread safety.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public interface HelpService {

    /**
     * <p>
     * Get all help topics.
     * </p>
     * @return the topics
     *
     * @throws ServiceException if any error occurred during the operation
     */
    public List<HelpTopic> getHelpTopics() throws ServiceException;

    /**
     * <p>
     * This method is used to search help items.
     * </p>
     *
     * @param criteria the search criteria
     * @return the search result.
     *
     * @throws IllegalArgumentException if criteria is null or invalid
     * @throws ServiceException if any other error occurred during the operation
     */
    public SearchResult<HelpItem> searchHelpItems(HelpItemSearchCriteria criteria) throws ServiceException;

    /**
     * <p>
     * This method is used to retrieve a help item.
     * </p>
     *
     * @param id the ID of the help item
     *
     * @return the help item, null will be returned if there's no such entity.
     *
     * @throws IllegalArgumentException: if id is not positive
     * @throws ServiceException if any other error occurred during the operation
     *
     */
    public HelpItem getHelpItem(long id) throws ServiceException;
}

