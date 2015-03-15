/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services;

/**
 * <p>
 * This exception will be thrown by services to indicate that relevant entity
 * can not be found when updating or deleting an entity.
 * </p>
 *
 * <p>
 * <strong>Thread Safety: </strong> Thread Safety is not an issue for Exceptions
 * since exceptions is used only in one thread.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class EntityNotFoundException extends ServiceException {

    /**
     * <p>
     * This is the constructor of class <code>EntityNotFoundException</code> with message.
     * </p>
     * @param message the error message
     */
    public EntityNotFoundException(String message) {
        super(message);
    }

    /**
     * <p>
     * This is the constructor of class <code>EntityNotFoundException</code> with message
     * and error cause.
     * </p>
     * @param message the error message
     * @param cause the error cause.
     */
    public EntityNotFoundException(String message, Throwable cause) {
        super(message, cause);
    }
}
