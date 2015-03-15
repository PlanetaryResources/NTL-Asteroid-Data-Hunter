/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter;

import java.awt.Desktop;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.concurrent.Executors;

import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.webapp.WebAppContext;

/**
 * <p>
 * This is the entry-point of the application. It will launch embedded Jetty
 * server and launch browser.
 * </p>
 * 
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class AsteroidHunterApplication {

    /**
     * <p>
     * Represents the prefix of the URL.
     * </p>
     */
    private static final String URL_PREFIX = "http://localhost:";
    
    /**
     * <p>
     * Represents the default listen port.
     * </p>
     */
    private static final int DEFAULT_PORT = 8080;
    
    /**
     * <p>
     * Represents the context path.
     * </p>
     */
    private static final String CONTEXT_PATH = "/";

    /**
     * <p>
     * This is the entry-point of the application. It will launch embedded Jetty
     * server and launch browser.
     * </p>
     * 
     * @param args the main args.
     */
    public static void main(String[] args) {
        WebAppContext context = new WebAppContext();
        
        // default port is 8080
        int port = DEFAULT_PORT;
        
        if (args.length < 1) {
            System.err.println("You should set the war file path.");
            return;
        }
        context.setWar(args[0]);
        if (args.length < 2) {
            System.err.println("You should set the current data path.");
            return;
        }
        String dataDir = args[1] + File.separator;
        if (args.length < 3) {
            System.err.println("You should set detector file path.");
            return;
        }
        
        String jdbcURL = "jdbc:h2:file:" + dataDir + "db" + File.separator + "hunter;DB_CLOSE_DELAY=10";
        System.setProperty("jdbc.url", jdbcURL);
        System.setProperty("detector.asteroids_detection_exe", args[2]);
        System.setProperty("detector.base_directory",dataDir + "hunting-data");
        
        
        // check if port is specified as an argument
        if (args.length > 3) {
            port = Integer.parseInt(args[3]);
        }
        
        System.out.println("dataDIR=" + dataDir);
        System.out.println("jdbc.url=" + jdbcURL);
        
        // Startup embedded Jetty server
        Server server = new Server(port);
        
        context.setContextPath(CONTEXT_PATH);
        context.setParentLoaderPriority(true);
        
        server.setHandler(context);
        try {
            server.start();
        } catch (Exception e) {
            e.printStackTrace();
        }
        final int finalPort = port;
        // Launch browser
        Executors.newSingleThreadExecutor().submit(new Runnable() {
            public void run() {
                try {
                    Desktop.getDesktop().browse(new URI(URL_PREFIX + finalPort));
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (URISyntaxException e) {
                    e.printStackTrace();
                }
            }
        });

        // Wait for server thread to join
        try {
            server.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }
}

