/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import gov.nasa.asteroid.hunter.models.DetectionItem;
import gov.nasa.asteroid.hunter.models.DetectionItemFrame;
import gov.nasa.asteroid.hunter.models.DetectionItemFramePk;
import gov.nasa.asteroid.hunter.models.DetectionItemSearchCriteria;
import gov.nasa.asteroid.hunter.models.DetectionSession;
import gov.nasa.asteroid.hunter.models.SearchResult;
import gov.nasa.asteroid.hunter.services.AsteroidDetectionService;
import gov.nasa.asteroid.hunter.services.DetectionSessionService;
import gov.nasa.asteroid.hunter.services.ServiceException;

import java.io.File;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.List;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

import com.fasterxml.jackson.databind.ObjectMapper;

/**
 * <p>
 * This class is the unit tests for class <code>AsteroidDetectionServiceImpl</code>.
 * </p>
 * @author TCSASSEMBLER
 * @version 1.0
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(locations = { "classpath:applicationContext.xml" })
public class AsteroidDetectionServiceImplTest extends PersistenceBaseTest {
    /**
     * <p>
     * Represents the folder storing the test fits images.
     * </p>
     */
    private static final String FITS_IMAGE_PATH = "/fits_images/";

    /**
     * <p>
     * Represents the asteroid detection service for unit tests.
     * </p>
     */
    @Autowired
    private AsteroidDetectionService asteroidDetectionService;

    /**
     * <p>
     * Represents the detection session service.
     * </p>
     */
    @Autowired
    private DetectionSessionService detectionSessionService;

    /**
     * <p>
     * Represents the search detection items.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchDetectionItems() throws Exception {
        DetectionItem item = new DetectionItem();
        item.setSubmitted(true);
        item.setTimestamp(new Date());
        item.setNeo(true);
        List<DetectionItemFrame> frames = new ArrayList<DetectionItemFrame>();
        for (int i = 0; i < 4; i++) {
            DetectionItemFrame frame = new DetectionItemFrame();
            DetectionItemFramePk detectionItemFramePk = new DetectionItemFramePk();
            detectionItemFramePk.setFrame(i + 1);
            frame.setDetectionItemFramePk(detectionItemFramePk);
            frame.setDetectionItem(item);
            frame.setDeclination(111);
            frame.setObservationLatitude(123);
            frame.setObservationLongitude(456);
            frame.setRightAscension(1);
            frame.setRoughMagnitude(2);
            frame.setVisualizationImage("abc");
            frames.add(frame);
        }
        item.setFrames(frames);
        persist(item);

        DetectionItemSearchCriteria criteria = new DetectionItemSearchCriteria();
        criteria.setSubmitted(true);
        SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(criteria);
        assertTrue("Should have result", result.getValues().size() == 1);
        item = result.getValues().get(0);
        item.getTimestamp();
        assertTrue("should be neo", item.isNeo());
        assertEquals("Should have 4 frames.", 4, item.getFrames().size());
        for (int i = 0; i < 4; i++) {
            DetectionItemFrame frame = item.getFrames().get(i);
            frame.getDetectionItemFramePk().equals(frame.getDetectionItemFramePk());
            frame.getDetectionItemFramePk().hashCode();
            assertEquals("should equal.", item, frame.getDetectionItem());
            assertEquals("should equal.", i + 1, frame.getDetectionItemFramePk().getFrame());
            assertEquals("should equal.", "abc", frame.getVisualizationImage());
            assertTrue("should be equal", Double.compare(111, frame.getDeclination()) == 0);
            assertTrue("should be equal", Double.compare(123, frame.getObservationLatitude()) == 0);
            assertTrue("should be equal", Double.compare(456, frame.getObservationLongitude()) == 0);
            assertTrue("should be equal", Double.compare(1, frame.getRightAscension()) == 0);
            assertTrue("should be equal", Double.compare(2, frame.getRoughMagnitude()) == 0);
        }
    }

    /**
     * <p>
     * Tests the search detection items and search by date start.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchDetectionItemsSearchByDateStart() throws Exception {
        DetectionItem item = new DetectionItem();
        item.setTimestamp(new Date());
        persist(item);
        item = new DetectionItem();
        item.setTimestamp(new Date(new Date().getTime() - 86400 * 20 * 1000L));
        persist(item);

        DetectionItemSearchCriteria criteria = new DetectionItemSearchCriteria();
        criteria.setDateStart(new Date(new Date().getTime() - 86400 * 10 * 1000L));
        SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(criteria);
        assertEquals("Should have valid result.", 1, result.getValues().size());
    }

    /**
     * <p>
     * Tests the search detection items and search by date end.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchDetectionItemsSearchByDateEnd() throws Exception {
        DetectionItem item = new DetectionItem();
        item.setTimestamp(new Date());
        persist(item);
        item = new DetectionItem();
        item.setTimestamp(new Date(new Date().getTime() - 86400 * 20 * 1000L));
        persist(item);

        DetectionItemSearchCriteria criteria = new DetectionItemSearchCriteria();
        criteria.setDateEnd(new Date(new Date().getTime() - 86400 * 10 * 1000L));
        SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(criteria);
        assertEquals("Should have valid result.", 1, result.getValues().size());
    }

    /**
     * <p>
     * Tests the search detection items and search by hour start.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchDetectionItemsSearchByHourStart() throws Exception {
        DetectionItem item = new DetectionItem();
        Calendar calendar = Calendar.getInstance();
        calendar.setTime(new Date());
        calendar.set(Calendar.HOUR_OF_DAY, 1);
        item.setTimestamp(calendar.getTime());
        persist(item);
        item = new DetectionItem();
        calendar.setTime(new Date());
        calendar.set(Calendar.HOUR_OF_DAY, 10);
        item.setTimestamp(calendar.getTime());
        persist(item);

        DetectionItemSearchCriteria criteria = new DetectionItemSearchCriteria();
        criteria.setHourStart(5);
        SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(criteria);
        assertEquals("Should have valid result.", 1, result.getValues().size());
    }

    /**
     * <p>
     * Tests the search detection items and search by hour end.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchDetectionItemsSearchByHourEnd() throws Exception {
        DetectionItem item = new DetectionItem();
        Calendar calendar = Calendar.getInstance();
        calendar.setTime(new Date());
        calendar.set(Calendar.HOUR_OF_DAY, 1);
        item.setTimestamp(calendar.getTime());
        persist(item);
        item = new DetectionItem();
        calendar.setTime(new Date());
        calendar.set(Calendar.HOUR_OF_DAY, 10);
        item.setTimestamp(calendar.getTime());
        persist(item);

        DetectionItemSearchCriteria criteria = new DetectionItemSearchCriteria();
        criteria.setHourEnd(5);
        SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(criteria);
        assertEquals("Should have valid result.", 1, result.getValues().size());
    }

    /**
     * <p>
     * Tests the search detection items and search by minute start.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchDetectionItemsSearchByMinuteStart() throws Exception {
        DetectionItem item = new DetectionItem();
        Calendar calendar = Calendar.getInstance();
        calendar.setTime(new Date());
        calendar.set(Calendar.MINUTE, 1);
        item.setTimestamp(calendar.getTime());
        persist(item);
        item = new DetectionItem();
        calendar.setTime(new Date());
        calendar.set(Calendar.MINUTE, 10);
        item.setTimestamp(calendar.getTime());
        persist(item);

        DetectionItemSearchCriteria criteria = new DetectionItemSearchCriteria();
        criteria.setMinuteStart(5);
        SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(criteria);
        assertTrue("Should have valid result.", result.getValues().size() > 0);
    }

    /**
     * <p>
     * Tests the search detection items and search by minute end.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSearchDetectionItemsSearchByMinuteEnd() throws Exception {
        DetectionItem item = new DetectionItem();
        Calendar calendar = Calendar.getInstance();
        calendar.setTime(new Date());
        calendar.set(Calendar.MINUTE, 1);
        item.setTimestamp(calendar.getTime());
        persist(item);
        item = new DetectionItem();
        calendar.setTime(new Date());
        calendar.set(Calendar.MINUTE, 10);
        item.setTimestamp(calendar.getTime());
        persist(item);

        DetectionItemSearchCriteria criteria = new DetectionItemSearchCriteria();
        criteria.setMinuteEnd(5);
        SearchResult<DetectionItem> result = asteroidDetectionService.searchDetectionItems(criteria);
        assertTrue("Should have valid result.", result.getValues().size() > 0);
    }

    /**
     * <p>
     * Tests the get detection item.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testGetDetectionItem() throws Exception {
        DetectionItem item = new DetectionItem();
        item.setSubmitted(true);
        item.setTimestamp(new Date());
        persist(item);

        long itemId = item.getId();
        item = asteroidDetectionService.getDetectionItem(itemId);
        assertNotNull("Item should exists", item);
    }

    /**
     * <p>
     * Tests submit detection item.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testSubmitDetectionItem() throws Exception {
        DetectionItem item = new DetectionItem();
        item.setTimestamp(new Date());
        persist(item);

        asteroidDetectionService.submitDetectionItem(item.getId());

        item = asteroidDetectionService.getDetectionItem(item.getId());

        assertTrue("item should mark submitted.", item.isSubmitted());
    }

    /**
     * <p>
     * Tests the detect asteroids.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testDetectAsteroids() throws Exception {
        List<File> fitsImages = scanFitsImages(FITS_IMAGE_PATH);
        for (File file : fitsImages) {
            System.out.println("The fits file:" + file.getAbsolutePath());
        }
        // create a new session
        final DetectionSession session = detectionSessionService.newDetectionSession();

        // start another thread to check the progress
        new Thread(new Runnable() {

            @Override
            public void run() {
                System.out.println("Progress checker threadto check the progress...");
                float prevProgress = 0;
                while (true) {
                    // read the progress
                    DetectionSession newSession = null;
                    try {
                        newSession = detectionSessionService.getDetectionSession(session.getId());
                    } catch (ServiceException e) {
                        // ignore
                    }
                    if (newSession == null) {
                        // ignore
                        continue;
                    }
                    if (Double.compare(newSession.getProgress(), prevProgress) != 0 ) {
                        System.out.println("Detected Progress updated to " + newSession.getProgress());
                        prevProgress = newSession.getProgress();

                    }
                    if (Double.compare(newSession.getProgress(), 1) >= 0) {
                        break;
                    }

                    try {
                        Thread.sleep(500);
                    } catch (InterruptedException e) {
                        // ignore
                    }
                }
                System.out.println("Progress checker thread is done.");
            }

        }).start();



        List<DetectionItem> result = asteroidDetectionService.detectAsteroids(
                fitsImages, session.getId(), false, "703");
        assertTrue("should contains result", result.size() > 0);
        // check if already save into db
        for (DetectionItem item : result) {
            DetectionItem dbResult = asteroidDetectionService.getDetectionItem(item.getId());
            assertNotNull("should exists in db.", dbResult);
        }

        // print the result
        ObjectMapper om = new ObjectMapper();
        String json = om.writeValueAsString(result);
        System.out.println("The detection result: " + json);
    }



    /**
     * <p>
     * Tests the mark detection item known by MPC.
     * </p>
     * @throws Exception to JUnit.
     */
    @Test
    public void testMarkDetectionItemKnownByMPC() throws Exception {
        DetectionItem item = new DetectionItem();
        item.setTimestamp(new Date());
        persist(item);

        asteroidDetectionService.markDetectionItemKnownByMPC(item.getId());

        item = asteroidDetectionService.getDetectionItem(item.getId());

        assertTrue("item should mark known by MPC.", item.isKnownByMPC());
    }

    /**
     * <p>
     * Helper method to scan fits images in a folder.
     * </p>
     * @param path the folder.
     * @return the list of fits images.
     */
    private List<File> scanFitsImages(String path) {
        path = getClass().getResource(path).getFile();
        File file = new File(path);
        List<File> result = new ArrayList<File>();
        for (String fileName : file.list()) {
            result.add(new File(path + File.separator + fileName));
        }
        return result;
    }
}
