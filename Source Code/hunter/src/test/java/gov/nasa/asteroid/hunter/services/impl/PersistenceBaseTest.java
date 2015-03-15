/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import java.io.BufferedReader;
import java.io.FileReader;

import javax.persistence.EntityManager;
import javax.persistence.EntityManagerFactory;
import javax.persistence.Query;

import org.junit.After;
import org.junit.Before;
import org.springframework.beans.factory.annotation.Autowired;

/**
 * <p>
 * This class is the helper base class for unit tests.
 * </p>
 * @author TCSASSEMBLER
 * @version 1.0
 */
public class PersistenceBaseTest {

    /**
     * <p>
     * Represents the sqls to clear the db.
     * </p>
     */
    private static final String CLEAR_DB_FILE = "/clear.sql";

    /**
     * <p>
     * Represents the entity manager factory. to create the entity manger.
     * </p>
     */
    @Autowired
    private EntityManagerFactory entityManagerFactory;

    /**
     * <p>
     * Represents the entity manager.
     * </p>
     */
    private EntityManager entityManager;

    /**
     * <p>
     * Sets up the test environment.
     * </p>
     * @throws Exception to JUnit.
     */
    @Before
    public void setUp() throws Exception {
        if (entityManager == null) {
            entityManager = entityManagerFactory.createEntityManager();
        }

        clearDB();
    }

    /**
     * <p>
     * Tears down. the test environment.
     * </p>
     * @throws Exception to JUnit.
     */
    @After
    public void tearDown() throws Exception {
        // clear the db
        clearDB();
    }

    /**
     * <p>
     * Clears the db.
     * </p>
     * @throws Exception to JUnit.
     */
    protected void clearDB() throws Exception {
        BufferedReader br = null;

        entityManager.getTransaction().begin();
        try {
            br = new BufferedReader(new FileReader(getClass().getResource(CLEAR_DB_FILE).getFile()));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().length() == 0) {
                    continue;
                }
                if (line.trim().charAt(0) == '-') {
                    // comment
                    continue;
                }
                // execute the sql

                Query query = entityManager.createNativeQuery(line.trim());

                // execute the update
                query.executeUpdate();
            }

            entityManager.getTransaction().commit();

        } finally {
            if (br != null) {
                br.close();
            }
        }
    }

    /**
     * <p>
     * Persist the entity.
     * </p>
     * @param entity entity to persist.
     */
    protected <T> void persist(T entity) {
        entityManager.getTransaction().begin();
        entityManager.persist(entity);
        entityManager.getTransaction().commit();
    }

    /**
     * <p>
     * Merges the entity.
     * </p>
     * @param entity  the entity to merge.
     */
    protected <T> void merge(T entity) {
        entityManager.getTransaction().begin();
        entityManager.merge(entity);
        entityManager.getTransaction().commit();
    }
}
