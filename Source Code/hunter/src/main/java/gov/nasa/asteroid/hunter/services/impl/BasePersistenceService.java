/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import gov.nasa.asteroid.hunter.ConfigurationException;
import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.models.BaseSearchParameters;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.models.SortType;
import gov.nasa.asteroid.hunter.services.ServiceException;

import java.util.List;
import java.util.Map;

import javax.annotation.PostConstruct;
import javax.persistence.EntityExistsException;
import javax.persistence.EntityManager;
import javax.persistence.PersistenceContext;
import javax.persistence.TransactionRequiredException;
import javax.persistence.TypedQuery;

import org.apache.log4j.Logger;

/**
 * <p>
 * Base class of service implementations that need to access database
 * persistence. It holds an EntityManager.
 * </p>
 *
 * <p>
 * <strong>Thread Safety:</strong> This class is effectively thread safe
 * (injected configurations are not considered as thread safety factor).
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public abstract class BasePersistenceService extends BaseService {

    /**
     * <p>
     * Represents the name of this class for logging.
     * </p>
     */
    private static final String CLASS_NAME = BasePersistenceService.class.getName();

    /**
     * <p>
     * Represents the query to search items.
     * </p>
     */
    private static final String QUERY_SEARCH = "SELECT e FROM %1$s e WHERE %2$s ORDER BY %3$s %4$s";

    /**
     * <p>
     * Represents the query to count the items.
     * </p>
     */
    private static final String QUERY_COUNT = "SELECT COUNT(e) FROM %1$s e WHERE %2$s";

    /**
     * <p>
     * Represents the EntityManager used to access database persistence.
     * </p>
     *
     * Required. Not null.
     */
    @PersistenceContext(unitName = "HunterPersistenceUnit")
    private EntityManager entityManager;



    /**
     * <p>
     * This is the default class for <code>BasePersistenceService</code>.
     * </p>
     */
    protected BasePersistenceService() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     * </p>
     * Note, in this class, entityManager is required.
     *
     * @throws ConfigurationException
     *             if any required field is not initialized properly.
     */
    @Override
    @PostConstruct
    protected void checkConfiguration() {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".checkConfiguration()";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);

        // call the supper class
        super.checkConfiguration();

        // validate the entity manager
        Helper.checkConfigurationNull(logger, signature, "entityManager", entityManager);

        // log the exit
        LoggingHelper.logExit(logger, signature, null);
    }

    /**
     * <p>
     * Gets the EntityManager used to access database persistence.
     * </p>
     * @return the EntityManager used to access database persistence.
     */
    protected EntityManager getEntityManager() {
        return entityManager;
    }

    /**
     * <p>
     * Sets the entity manager.
     * </p>
     * @param entityManager the entity manager.
     */
    public void setEntityManager(EntityManager entityManager) {
        this.entityManager = entityManager;
    }

    /**
     * <p>
     * This method is used to search entities based on the search criteria.
     * </p>
     * @param <T> the type of the entity
     * @param searchParameters
     *            the base search parameters
     * @param whereClause
     *            the JQL where clause
     * @param queryParameters
     *            the query parameters map
     * @param entityClass
     *            the entity class
     * @return the search result
     * @throws ServiceException
     *             if any error occurred during the operation.
     * @throws IllegalArgumentException
     *            if the searchParameters is null.
     */
    protected <T> SearchResult<T> search(BaseSearchParameters searchParameters, String whereClause,
            Map<String, Object> queryParameters, Class<T> entityClass) throws ServiceException {

        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".search(BaseSearchParameters, String, Map, Class)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature,
                new String[] {"searchParameters", "whereClause", "queryParameters", "entityClass"},
                new Object[] {searchParameters, whereClause, queryParameters, entityClass});

        Helper.checkNull(logger, signature, "searchParameters", searchParameters);

        // Set default sorting option
        if (searchParameters.getSortBy() == null) {
            searchParameters.setSortBy("id");
        }
        if (searchParameters.getSortType() == null) {
            searchParameters.setSortType(SortType.ASC);
        }

        EntityManager entityManager = getEntityManager();

        String searchQueryString = String.format(QUERY_SEARCH, entityClass.getSimpleName(), whereClause,
                searchParameters.getSortBy(), searchParameters.getSortType().name());

        // Create query
        TypedQuery<T> query = entityManager.createQuery(searchQueryString, entityClass);

        // set the parameters
        setQueryParameters(logger, signature, query, queryParameters);


        // set paging options
        if (searchParameters.getPageNumber() > 0) {
            query.setMaxResults(searchParameters.getPageSize());
            query.setFirstResult((searchParameters.getPageNumber() - 1) * searchParameters.getPageSize());
        }

        List<T> list = Helper.getResultList(logger, signature, query);

        int totalCount = list.size();
        int totalPageCount = totalCount > 0 ? 1 : 0;

        if (searchParameters.getPageNumber() > 0) {

            String countQueryString = String.format(QUERY_COUNT, entityClass.getSimpleName(), whereClause);

            // Create query
            TypedQuery<Long> countQuery = entityManager.createQuery(countQueryString, Long.class);

            // Set query parameters
            setQueryParameters(logger, signature, countQuery, queryParameters);

            long totalCountLong = Helper.getSingleResult(logger, signature, countQuery);
            totalCount = (int) totalCountLong;

            totalPageCount = (totalCount + searchParameters.getPageSize() - 1) / searchParameters.getPageSize();
        }

        SearchResult<T> result = new SearchResult<T>();

        result.setSortBy(searchParameters.getSortBy());
        result.setSortType(searchParameters.getSortType());
        result.setValues(list);
        result.setTotal(totalCount);
        result.setTotalPages(totalPageCount);
        result.setPageNumber(searchParameters.getPageNumber());
        result.setPageSize(searchParameters.getPageSize());

        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] {result});

        return result;
    }

    protected <T> SearchResult<T> searchWithQuery(BaseSearchParameters searchParameters, String queryString,
            String countQueryString, Map<String, Object> queryParameters, Class<T> entityClass) throws ServiceException {

        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".searchWithQuery(BaseSearchParameters, String, String, Map, Class)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature,
                new String[] {"searchParameters", "queryString", "countQueryString", "queryParameters", "entityClass"},
                new Object[] {searchParameters, queryString, countQueryString, queryParameters, entityClass});

        Helper.checkNull(logger, signature, "searchParameters", searchParameters);


        EntityManager entityManager = getEntityManager();


        // Create query
        TypedQuery<T> query = entityManager.createQuery(queryString, entityClass);

        // set the parameters
        setQueryParameters(logger, signature, query, queryParameters);


        // set paging options
        if (searchParameters.getPageNumber() > 0) {
            query.setMaxResults(searchParameters.getPageSize());
            query.setFirstResult((searchParameters.getPageNumber() - 1) * searchParameters.getPageSize());
        }

        List<T> list = Helper.getResultList(logger, signature, query);

        int totalCount = list.size();
        int totalPageCount = totalCount > 0 ? 1 : 0;

        if (searchParameters.getPageNumber() > 0) {
            
            // Create query
            TypedQuery<Long> countQuery = entityManager.createQuery(countQueryString, Long.class);

            // Set query parameters
            setQueryParameters(logger, signature, countQuery, queryParameters);

            long totalCountLong = Helper.getSingleResult(logger, signature, countQuery);
            totalCount = (int) totalCountLong;

            totalPageCount = (totalCount + searchParameters.getPageSize() - 1) / searchParameters.getPageSize();
        }

        SearchResult<T> result = new SearchResult<T>();

        result.setSortBy(searchParameters.getSortBy());
        result.setSortType(searchParameters.getSortType());
        result.setValues(list);
        result.setTotal(totalCount);
        result.setTotalPages(totalPageCount);
        result.setPageNumber(searchParameters.getPageNumber());
        result.setPageSize(searchParameters.getPageSize());

        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] {result});

        return result;
    }
    
    /**
     * <p>
     * Sets the parameters to the query.
     * </p>
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param query the query.
     * @param queryParameters the parameters.
     * @throws ServiceException if there are any errors when setting the parameters.
     */
    private static void setQueryParameters(Logger logger, String signature, TypedQuery<?> query,
            Map<String, Object> queryParameters) throws ServiceException {
        try {
            for (Map.Entry<String, Object> entry : queryParameters.entrySet()) {
                query.setParameter(entry.getKey(), entry.getValue());
            }
        } catch (IllegalArgumentException e) {
            // parameter not found, wrap it to ServiceException
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Failed to set the query parameters", e));
        }
    }

    /**
     * <p>
     * Merge the entity. This is a helper wrapper method of JPA merge.
     * </p>
     *
     * @param <T>
     *            the type of the entity.
     * @param logger
     *            the logger for logging.
     * @param signature
     *            the signature of the caller method for logging.
     * @param entity
     *            the entity to merge.
     * @throws ServiceException
     *             if there are any error.
     */
    protected <T> void merge(Logger logger, String signature, T entity) throws ServiceException {
        try {
            entityManager.merge(entity);
        } catch (IllegalArgumentException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException("Invalid entity", e));
        } catch (TransactionRequiredException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Transaction is required.", e));
        }
    }

    /**
     * <p>
     * Merge the entity. This is a helper wrapper method of JPA persist.
     * </p>
     *
     * @param <T>
     *            the type of the entity.
     * @param logger
     *            the logger for logging.
     * @param signature
     *            the signature of the caller method for logging.
     * @param entity
     *            the entity to persist.
     * @throws ServiceException
     *             if there are any error.
     */
    protected <T> void persist(Logger logger, String signature, T entity) throws ServiceException {
        try {
            entityManager.persist(entity);
        } catch (EntityExistsException e) {
            throw LoggingHelper
                    .logException(logger, signature, new ServiceException("Enity already exists", e));
        } catch (IllegalArgumentException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Enity is not a Hibernate entity", e));
        } catch (TransactionRequiredException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Transaction is required.", e));
        }
    }
}
