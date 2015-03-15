/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services;

import gov.nasa.asteroid.hunter.models.DetectionItem;
import gov.nasa.asteroid.hunter.models.DetectionItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.tester.ImageAlreadyExistException;

import java.io.File;
import java.util.List;

/**
 * <p>
 * This service provides methods to detect asteroids and access asteroid
 * detection results.
 * </p>
 * <p>
 * <b>Thread Safety:</b> Implementations must be effectively thread safe.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public interface AsteroidDetectionService {

    /**
     * <p>
     * This method is used to search asteroid detection items.
     * </p>
     *
     * @param criteria
     *            the search criteria.
     *
     * @return the search result.
     *
     * @throws IllegalArgumentException
     *             if criteria is null or invalid (see
     *             <code>DetectionItemSearchCriteria</code> document for
     *             details).
     * @throws ServiceException
     *             if any other error occurred during the operation
     *
     */
    public SearchResult<DetectionItem> searchDetectionItems(DetectionItemSearchCriteria criteria)
            throws ServiceException;

    /**
     * <p>
     * This method is used to retrieve a detection item.
     * </p>
     *
     * @param id
     *            the ID
     * @return the detection item, null will be returned if there is no such
     *         entity.
     * @throws ServiceException
     *             if there are any errors.
     * @throws IllegalArgumentException
     *             if the id is not positive.
     */
    public DetectionItem getDetectionItem(long id) throws ServiceException;

    /**
     * <p>
     * This method is used to submit a detection item (mark the detection item
     * as submitted).
     * </p>
     *
     * @param id
     *            the ID
     *
     * @throws IllegalArgumentException
     *             if id is not positive
     * @throws EntityNotFoundException
     *             if there is no such detection item
     * @throws ServiceException
     *             if any other error occurred during the operation
     *
     */
    public void submitDetectionItem(long id) throws ServiceException;

    /**
     * <p>
     * This method is used to mark a detection item as known by MPC.
     * </p>
     *
     * @param id
     *            the ID
     *
     * @throws IllegalArgumentException
     *             if id is not positive
     * @throws EntityNotFoundException
     *             if there is no such detection item
     * @throws ServiceException
     *             if any other error occurred during the operation
     *
     */
    public void markDetectionItemKnownByMPC(long id) throws ServiceException;

    /**
     * <p>
     * This method is used to detect asteroids in FITS images.
     * </p>
     *
     * <p>
     * This method should use READ_UNCOMMITTED transaction isolation level so
     * that progress information can be read while this method is executing.
     * </p>
     *
     * @param fitsImages
     *            the FITS image files, must be an array of exactly 4 items
     * @param detectionSessionId
     *            the ID of the detection session.
     * @param forced 
     * @param observatoryCode 
     * @return the detection result
     *
     * @throws IllegalArgumentException
     *             if fitsImages is null, or does not contain 4 elements, or any
     *             item is null or detectionSessionId is not positive.
     * @throws ServiceException
     *             if any other error occurred during the operation
     * @throws ImageAlreadyExistException 
     *
     */
    public List<DetectionItem> detectAsteroids(List<File> fitsImages, long detectionSessionId, boolean forced, String observatoryCode)
            throws ServiceException, ImageAlreadyExistException;

    void deleteAllDetectionItems() throws ServiceException;

}
