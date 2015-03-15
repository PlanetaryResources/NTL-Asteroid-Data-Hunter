package gov.nasa.asteroid.tester;


public interface ImageService {
    public long addImage(AsteroidDetectorTester.FITSImage image);
    public boolean imageExists(AsteroidDetectorTester.FITSImage image);
    public <T> T getLatestImage(Class<T> clazz);
}
