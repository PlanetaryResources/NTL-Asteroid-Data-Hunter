/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.models.HelpItem;
import gov.nasa.asteroid.hunter.models.HelpItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.HelpTopic;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.services.HelpService;
import gov.nasa.asteroid.hunter.services.ServiceException;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

/**
 * <p>
 * This class is the implementation of HelpService.
 * </p>
 *
 * <p>
 * This service provides methods to access help contents.
 * </p>
 *
 * <p>
 * <b>Thread Safety:</b> This class is effectively thread safe (injected
 * configurations are not considered as thread safety factor).
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class HelpServiceImpl extends BasePersistenceService implements HelpService {

    /**
     * <p>
     * Represents the name of the class for logging.
     * </p>
     */
    private static final String CLASS_NAME = HelpServiceImpl.class.getName();

    /**
     * <p>
     * Represents the percent mark use for like query in JQL.
     * </p>
     */
    private static final String PERCENT_MARK = "%";

    /**
     * <p>
     * Represents the parameter name of the keyword.
     * </p>
     */
    private static final String KEYWORD_PARAMETER_NAME = "keyword";

    /**
     * <p>
     * Represents the parameter name of the topic id.
     * </p>
     */
    private static final String TOPIC_ID_PARAMETER_NAME = "topicId";

    /**
     * <p>
     * Represents the keyword where cause.
     * </p>
     */
    private static final String KEYWORD_WHERE_CAUSE = 
            " AND (UPPER(e.content) LIKE UPPER(:keyword) OR UPPER(e.title) LIKE UPPER(:keyword))";

    /**
     * <p>
     * Represents the topic id where cause.
     * </p>
     */
    private static final String TOPIC_ID_WHERE_CAUSE = " AND e.topic.id = :topicId";

    /**
     * <p>
     * Represents the prefix of the where cause.
     * </p>
     */
    private static final String WHERE_CAUSE_PREFIX = "1 = 1";

    /**
     * <p>
     * Represents the query to select the topics.
     * </p>
     */
    private static final String SELECT_TOPIC_QUERY = "SELECT t FROM HelpTopic t";

    /**
     * <p>
     * This is the default constructor of <code>HelpServiceImpl</code>.
     * </p>
     */
    public HelpServiceImpl() {
        // does nothing
    }

    /**
     * <p>
     * Get all help topics.
     * </p>
     * @return the topics
     *
     * @throws ServiceException if any error occurred during the operation
     */
    @Override
    public List<HelpTopic> getHelpTopics() throws ServiceException {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".getHelpTopics()";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);

        try {
            List<HelpTopic> result = Helper.getResultList(logger, signature,
                    getEntityManager().createQuery(SELECT_TOPIC_QUERY, HelpTopic.class));

            // log the exit
            LoggingHelper.logExit(logger, signature, new Object[] {result});

            return result;
        } catch (ServiceException e) {
            // log the error and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }

    }

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
     *
     */
    @Override
    public SearchResult<HelpItem> searchHelpItems(HelpItemSearchCriteria criteria)
        throws ServiceException {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".searchHelperItems(HelpItemSearchCriteria)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature,
                new String[] {"criteria"}, new Object[] {criteria});

        // validate the parameters
        Helper.checkBaseSearchParameters(logger, signature, "criteria", criteria,
                Arrays.asList("id", "title", "content"));

        StringBuffer sb = new StringBuffer(WHERE_CAUSE_PREFIX);
        Map<String, Object> queryParameters = new HashMap<String, Object>();

        // Append criteria
        if (criteria.getTopicId() != null) {
            sb.append(TOPIC_ID_WHERE_CAUSE);
            queryParameters.put(TOPIC_ID_PARAMETER_NAME, criteria.getTopicId());
        }
        if (criteria.getKeyword() != null && criteria.getKeyword().trim().length() > 0) {
            sb.append(KEYWORD_WHERE_CAUSE);
            queryParameters.put(KEYWORD_PARAMETER_NAME, PERCENT_MARK + criteria.getKeyword() + PERCENT_MARK);
        }

        try {
            SearchResult<HelpItem> result = search(criteria, sb.toString(), queryParameters, HelpItem.class);

            // log the exit
            LoggingHelper.logExit(logger, signature, new Object[] {result});

            return result;
        } catch (ServiceException e) {
            // log the exception and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
}

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
    @Override
    public HelpItem getHelpItem(long id) throws ServiceException {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".getHelpItem(long)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "id" }, new Object[] { id });

        // check the parameters
        Helper.checkPositive(logger, signature, "id", id);

        try {
            HelpItem result = getEntityManager().find(HelpItem.class, id);

            // log the exit
            LoggingHelper.logExit(logger, signature, new Object[] { result });

            return result;
        } catch (IllegalArgumentException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException("Failed to get the result.", e));
        }
    }
}

