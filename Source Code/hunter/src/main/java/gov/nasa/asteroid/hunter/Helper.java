/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter;

import gov.nasa.asteroid.hunter.models.BaseSearchParameters;
import gov.nasa.asteroid.hunter.services.ServiceException;

import java.util.List;

import javax.persistence.LockTimeoutException;
import javax.persistence.NoResultException;
import javax.persistence.NonUniqueResultException;
import javax.persistence.PersistenceException;
import javax.persistence.PessimisticLockException;
import javax.persistence.QueryTimeoutException;
import javax.persistence.TransactionRequiredException;
import javax.persistence.TypedQuery;

import org.apache.log4j.Logger;


/**
 * <p>
 * This is a helper class for the package: gov.nasa.asteroid.hunter.services.impl.
 * </p>
 * <p>
 * This class provides the methods for configuration validation and parameters validation.
 * </p>
 *
 * <p>
 * <strong>Thread-Safety</strong>: This class is immutable so it is thread-safe.
 * </p>
 *
 * @author TCSASSEMBLER
 * @version 1.0
 *
 */
public final class Helper {

    /**
     * <p>
     * This is the default constructor for <code>ServiceHelper.</code>
     * </p>
     * <p>
     * This is a private constructor to prevent creating instances.
     * </p>
     */
    private Helper() {
        // does nothing
    }

    /**
     * <p>
     * Checks if a property configuration is null or not.
     * </p>
     *
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param name the name of the property.
     * @param value the value of the property.
     *
     * @throws ConfigurationException if the property is null.
     *
     */
    public static void checkConfigurationNull(Logger logger, String signature, String name, Object value) {
        if (value == null) {
            throw LoggingHelper.logException(logger, signature, new ConfigurationException("Property " + name
                    + " is not properly injected: it cannot be null."));
        }
    }

    /**
     * <p>
     * Checks if a string property configuration is null and empty or not.
     * </p>
     *
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param name the name of the property.
     * @param value the value of the property.
     *
     * @throws ConfigurationException if the property is null or empty.
     *
     */
    public static void checkConfigurationNullOrEmpty(Logger logger, String signature, String name, String value) {
        checkConfigurationNull(logger, signature, name, value);
        if (value.trim().length() == 0) {
            throw LoggingHelper.logException(logger, signature, new ConfigurationException("Property " + name
                    + " is not properly injected: it cannot be empty."));
        }
    }


    /**
     * <p>
     * Checks the parameter is null or not.
     * </p>
     *
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param name the name of the parameter.
     * @param value the value of the parameter.
     *
     * @throws IllegalArgumentException if the parameter is null.
     */
    public static void checkNull(Logger logger, String signature, String name, Object value) {
        if (value == null) {
            throw LoggingHelper.logException(logger, signature, new IllegalArgumentException("Parameter " + name
                    + " cannot be null."));
        }
    }

    /**
     * <p>
     * Checks the parameter is null and empty or not.
     * </p>
     *
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param name the name of the parameter.
     * @param value the value of the parameter.
     *
     * @throws IllegalArgumentException if the parameter is null or empty.
     */
    public static void checkNullOrEmpty(Logger logger, String signature, String name, String value) {
        checkNull(logger, signature, name, value);
        if (value.trim().length() == 0) {
            throw LoggingHelper.logException(logger, signature, new IllegalArgumentException("Parameter " + name
                    + " cannot be empty."));
        }
    }

    /**
     * <p>
     * Checks the parameter is positive or not.
     * </p>
     *
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param name the name of the parameter.
     * @param value the value of the parameter.
     *
     * @throws IllegalArgumentException if the parameter is not positive.
     */
    public static void checkPositive(Logger logger, String signature, String name, long value) {
        if (value <= 0) {
            throw LoggingHelper.logException(logger, signature, new IllegalArgumentException("Property " + name
                    + " must be positive. current value " + value));
        }
    }

    /**
     * <p>
     * A helper method to get result list from a query.
     * </p>
     * @param <T> the type of the query entity
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param query the query.
     * @return the result list
     * @throws ServiceException if there are any error occurs.
     */
    public static <T> List<T> getResultList(Logger logger, String signature, TypedQuery<T> query)
        throws ServiceException {
        try {
            return query.getResultList();
        } catch (IllegalStateException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (QueryTimeoutException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (TransactionRequiredException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (PessimisticLockException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (LockTimeoutException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (PersistenceException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        }
    }

    /**
     * <p>
     * A helper method to get single result from a query.
     * </p>
     * @param <T> the type of the entity.
     * @param logger the logger for logging.
     * @param signature the signature of the caller method for logging.
     * @param query the query.
     * @return the single result
     * @throws ServiceException if there are any error occurs.
     */
    public static <T> T getSingleResult(Logger logger, String signature, TypedQuery<T> query)
        throws ServiceException {
        try {
            return query.getSingleResult();
        } catch (NoResultException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (NonUniqueResultException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (IllegalStateException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (QueryTimeoutException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (TransactionRequiredException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (PessimisticLockException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (LockTimeoutException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        } catch (PersistenceException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException(
                    "Error occurs when query the list.", e));
        }
    }

    /**
     * <p>
     * Validates if the base search criteria is valid or not.
     * </p>
     * @param logger the logger for logging.
     * @param signature the signature for logging.
     * @param name the parameter name.
     * @param criteria the criteria to check.
     * @param validSortByFields the valid field names for sort by.
     * @throws IllegalArgumentException if the criteria is invalid.
     * (You can see the <code>BaseSearchParameters</code> document to see the valid criteria properties).
     *
     */
    public static void checkBaseSearchParameters(Logger logger, String signature,
            String name, BaseSearchParameters criteria, List<String> validSortByFields) {
        checkNull(logger, signature, name, criteria);
        if (criteria.getPageNumber() < 0) {
            throw LoggingHelper.logException(logger, signature,
                    new IllegalArgumentException("The pageNumber must be non-negative."));
        }
        if (criteria.getPageNumber() > 0 && criteria.getPageSize() <= 0) {
            throw LoggingHelper.logException(logger, signature,
                    new IllegalArgumentException("The pageSize must be positive if paging enabled."));
        }
        if (criteria.getSortBy() != null && !validSortByFields.contains(criteria.getSortBy())) {
            throw LoggingHelper.logException(logger, signature,
                    new IllegalArgumentException("The sort by field is invalid."));
        }
    }
}
