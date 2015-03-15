/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */

package gov.nasa.asteroid.hunter.services.impl;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import gov.nasa.asteroid.hunter.models.DetectionSession;
import gov.nasa.asteroid.hunter.services.DetectionSessionService;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

/**
 * <p>
 * This class is the unit tests for class <code>DetectionSessionServiceImpl</code>.
 * </p>
 * @author TCSASSEMBLER
 * @version 1.0
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(locations = { "classpath:applicationContext.xml" })

public class DetectionSessionServiceImplTest extends PersistenceBaseTest {

    /**
     * <p>
     * Represents the detection session service.
     * </p>
     */
    @Autowired
    private DetectionSessionService detectionSessionService;

    /**
     * <p>
     * Tests to get detection session.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testGetDetectionSession() throws Exception {
        DetectionSession session = detectionSessionService.newDetectionSession();
        session = detectionSessionService.getDetectionSession(session.getId());
        assertNotNull("Session should exist.", session);
    }

    /**
     * <p>
     * Tests the new detection session.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testNewDetectionSession() throws Exception {
        DetectionSession session = detectionSessionService.newDetectionSession();

        session.setProgress(0.5f);
        merge(session);

        session = detectionSessionService.getDetectionSession(session.getId());
        assertTrue("should be equal", Double.compare(0.5, session.getProgress()) == 0);

        assertNotNull("Session should created.", session);
    }

    /**
     * <p>
     * Tests for the detection update.
     * </p>
     * @throws Exception
     */
    @Test
    public void testUpdateDetectionSession() throws Exception {
        DetectionSession session = detectionSessionService.newDetectionSession();
        detectionSessionService.updateDetectionSessionProgress(session.getId(), 0.5f);

        session = detectionSessionService.getDetectionSession(session.getId());
        assertTrue("should be equal", Double.compare(0.5, session.getProgress()) == 0);

        assertNotNull("Session should created.", session);
    }

}
