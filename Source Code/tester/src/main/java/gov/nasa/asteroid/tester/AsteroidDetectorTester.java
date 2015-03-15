package gov.nasa.asteroid.tester;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.security.SecureRandom;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Random;

import javax.imageio.ImageIO;

public class AsteroidDetectorTester
{
    private boolean debug = true;

    private String trainDataFolder = null;
    private boolean visualize = false;
    private int visnr = 0;
    private String execCommand = null;
    private String trainFile = null;
    private String testFile = null;
    private String folder = "";
    private int maxNumberOfAnswer = 10;
    private boolean forced = false;
    private Date imageCreateTime = null;
    private ImageSet currentImageSet = null;
    private String observatoryCode = null;
    private String fileExtension = null;
    private static final int NUM_OF_TESTS = 1;
    
    /**
     * Represents the detection progress listener.
    */
    private DetectionProgressListener detectionProgressListener;
    
    private ImageService imageService;
    
    static double MATCH_DISTANCE = 0.001;
    static String[] WCSfields = {
      "CRPIX1",
      "CRPIX2",
      "CRVAL1",
      "CRVAL2",
      "CD1_1",
      "CD1_2",
      "CD2_1",
      "CD2_2"
    };
    public static long seed = 1;

    public void printMessage(String s) {
        if (debug) {
            System.out.println(s);
        }
    }

    static final String _ = File.separator;

    public class FITSImage {
        int[] pixData;
        String[] header;
        public int W, H;
        double[] WCS = new double[8];
        
        /**
         * Represents the observation latitude.
         */
        private double latitude;

        /**
         * Represents the observation longitude.
         */
        private double longitude;

        private String timeOBS;
        private String dateOBS;
        private String imageCreateTime;
        private String RA;
        private String DEC;
        private String observatoryCode;
        
        public double getLatitude() {
            return latitude;
        }

        public void setLatitude(double latitude) {
            this.latitude = latitude;
        }

        public double getLongitude() {
            return longitude;
        }

        public void setLongitude(double longitude) {
            this.longitude = longitude;
        }

        public String getTimeOBS() {
            return timeOBS;
        }

        public void setTimeOBS(String timeOBS) {
            this.timeOBS = timeOBS;
        }

        public String getDateOBS() {
            return dateOBS;
        }

        public void setDateOBS(String dateOBS) {
            this.dateOBS = dateOBS;
        }

        public String getImageCreateTime() {
            return imageCreateTime;
        }

        public void setImageCreateTime(String imageCreateTime) {
            this.imageCreateTime = imageCreateTime;
        }

        public String getRA() {
            return RA;
        }

        public void setRA(String rA) {
            RA = rA;
        }

        public String getDEC() {
            return DEC;
        }

        public void setDEC(String dEC) {
            DEC = dEC;
        }

        public String getObservatoryCode() {
            return observatoryCode;
        }

        public void setObservatoryCode(String observatoryCode) {
            this.observatoryCode = observatoryCode;
        }
        
       
    }

    public static class Entity {
        
        String ID;
        public double[] RA = new double[4];
        public double[] DEC = new double[4];
        public double[] mag = new double[4];
        public int[] x = new int[4];
        public int[] y = new int[4];
        public boolean isNeo;
        public double significance;
        boolean matched;
        
        /**
         * Represents the visualization image paths. Should contain four items.
         */
        private String[] visualizationImages = new String[4];

        /**
         * Represents the observation latitudes in four frames. Should contain four items.
         */
        public double[] latitudes = new double[4];

        /**
         * Represents the observation longitudes in four frames. Should contain four items.
         */
        public double[] longitudes = new double[4];

        public Entity() {}
        
        public String[] getVisualizationImages() {
            return visualizationImages;
        }

        public void setVisualizationImages(String[] visualizationImages) {
            this.visualizationImages = visualizationImages;
        }

        public double[] getLatitudes() {
            return latitudes;
        }

        public void setLatitudes(double[] latitudes) {
            this.latitudes = latitudes;
        }

        public double[] getLongitudes() {
            return longitudes;
        }

        public void setLongitudes(double[] longitudes) {
            this.longitudes = longitudes;
        }
        
    };

    public class ImageSet {
        public FITSImage[] img = new FITSImage[4];
        ArrayList<Entity> detections = new ArrayList<Entity>();
        int numNeos;
        String name;
    };

    class Wcs {

        double crpix1;
        double crpix2;
        double crval1;
        double crval2;
        double cd11;
        double cd12;
        double cd21;
        double cd22;
        double invcd11;
        double invcd12;
        double invcd21;
        double invcd22;
        double D2R = Math.PI/180.0;
        double R2D = 180.0/Math.PI;

        public void init(double[] wcs) {
            crpix1 = wcs[0];
            crpix2 = wcs[1];
            crval1 = wcs[2];
            crval2 = wcs[3];
            cd11 = wcs[4];
            cd12 = wcs[5];
            cd21 = wcs[6];
            cd22 = wcs[7];
            update_invcd();
        }

        void update_invcd() {
            double inv_det = cd11*cd22 - cd12*cd21;
            inv_det = 1.0 / inv_det;
            invcd11 = cd22 * inv_det;
            invcd12 = -cd12 * inv_det;
            invcd21 = -cd21 * inv_det;
            invcd22 = cd11 * inv_det;
        }


        public double[] sphs2x(double lng, double lat) {
            double coslat, coslat3, coslat4, coslng, dlng, dphi, sinlat, sinlat3, sinlat4, sinlng, x, y, z;
            double[] eul = new double[5];
            double[] phiTheta = new double[2];
            eul[0] = crval1;
            eul[1] = 90.0 - crval2;
            eul[2] = 180.0;
            eul[3] = Math.cos( D2R * eul[1] );
            eul[4] = Math.sin( D2R * eul[1] );
            /* Do lng dependency. */
            dlng = lng - eul[0];
            phiTheta[0] = dlng;
            /* Do lat dependency. */
            sinlat = Math.sin(lat*D2R);
            coslat = Math.cos(lat*D2R);
            coslat3 = coslat*eul[3];
            coslat4 = coslat*eul[4];
            sinlat3 = sinlat*eul[3];
            sinlat4 = sinlat*eul[4];
            dlng = phiTheta[0];
            sinlng = Math.sin(dlng*D2R);
            coslng = Math.cos(dlng*D2R);
            /* Compute the native longitude. */
            x = sinlat4 - coslat3*coslng;
            y = -coslat*sinlng;
            dphi = R2D*Math.atan2(y, x);
            phiTheta[0] = (eul[2] + dphi) % 360.0;
            /* Normalize the native longitude. */
            if (phiTheta[0] > 180.0) {
                phiTheta[0] -= 360.0;
            } else if (phiTheta[0] < -180.0) {
                phiTheta[0] += 360.0;
            }
            /* Compute the native latitude. */
            z = sinlat3 + coslat4*coslng;
            if (Math.abs(z) > 0.99) {
                if (z < 0.0)
                    phiTheta[1] = -Math.abs(R2D * Math.acos(Math.sqrt(x*x+y*y)));
                else
                    phiTheta[1]= Math.abs(R2D * Math.acos(Math.sqrt(x*x+y*y)));
            } else {
                phiTheta[1] = R2D * Math.asin(z);
            }
            return phiTheta;
        }

        public double[] sphx2s(double phi, double theta) {
            double cosphi, costhe, costhe3, costhe4, dlng, dphi, sinphi, sinthe, sinthe3, sinthe4, x, y, z;
            double[] lngLat = new double[2];
            double[] eul = new double[5];
            eul[0] = crval1;
            eul[1] = 90.0 - crval2;
            eul[2] = 180.0;
            eul[3] = Math.cos( D2R * eul[1] );
            eul[4] = Math.sin( D2R * eul[1] );
            /* Do phi dependency. */
            dphi = phi - eul[2];
            /* Do theta dependency. */
            sinthe = Math.sin( theta * D2R );
            costhe = Math.cos( theta * D2R );
            costhe3 = costhe*eul[3];
            costhe4 = costhe*eul[4];
            sinthe3 = sinthe*eul[3];
            sinthe4 = sinthe*eul[4];
            sinphi = Math.sin( dphi * D2R );
            cosphi = Math.cos( dphi * D2R );
            /* Compute the celestial longitude. */
            x = sinthe4 - costhe3*cosphi;
            y = -costhe*sinphi;
            dlng = R2D * Math.atan2( y, x);
            lngLat[0] = eul[0] + dlng;
            /* Normalize the celestial longitude. */
            if (eul[0] >= 0.0) {
                if (lngLat[0] < 0.0)
                    lngLat[0] += 360.0;
            } else {
                if (lngLat[0] > 0.0)
                    lngLat[0] -= 360.0;
            }
            if (lngLat[0] > 360.0) {
                lngLat[0] -= 360.0;
            } else if (lngLat[0] < -360.0) {
                lngLat[0] += 360.0;
            }
            /* Compute the celestial latitude. */
            z = sinthe3 + costhe4*cosphi;
            if (Math.abs(z) > 0.99) {
            /* Use an alternative formula for greater accuracy. */
                if (z<0.0)
                    lngLat[1] = -Math.abs( R2D * Math.acos(Math.sqrt(x*x+y*y)) );
                else
                    lngLat[1] = Math.abs( R2D * Math.acos(Math.sqrt(x*x+y*y)) );
            } else {
                lngLat[1] = R2D * Math.asin(z);
            }
            return lngLat;
        }


        public double[] tans2x(double phi, double theta) {
            double[] xy = new double[2];
            double cotan_theta = 1.0 / Math.tan( D2R * theta );
            xy[0] = R2D * Math.sin( D2R * phi ) * cotan_theta ;
            xy[1] = -R2D * Math.cos( D2R * phi ) * cotan_theta;
            return xy;
        }

        public double[] tanx2s(double x, double y) {
            double[] phiTheta = new double[2];
            phiTheta[0] = R2D * Math.atan2(x, -y);
            phiTheta[1] = R2D * Math.atan2(R2D, Math.sqrt(x*x + y*y) );
            return phiTheta;
        }

        public double[] convertRADEC2XY(FITSImage fits, double RA, double DEC) {
            init(fits.WCS);
            double[] XY = new double[2];
            double[] phiTheta = sphs2x(RA, DEC);
            double[] XXYY = tans2x(phiTheta[0], phiTheta[1]);
            XY[0] = invcd11*XXYY[0] + invcd12*XXYY[1] + crpix1;
            XY[1] = invcd21*XXYY[0] + invcd22*XXYY[1] + crpix2;
            return XY;
        }

        public double[] convertXY2RADEC(FITSImage fits, double x, double y) {
            init(fits.WCS);
            double dx = x-crpix1;
            double dy = y-crpix2;
            double xx = cd11*dx + cd12*dy;
            double yy = cd21*dx + cd22*dy;
            double[] phiTheta = tanx2s(xx, yy);
            double[] RD = sphx2s(phiTheta[0], phiTheta[1]);
            return RD;
        }
    };

    ArrayList<String> knownNEOSList = new ArrayList<String>();
    String[] knownNEOS;
    Wcs WcsConvert = new Wcs();

    int bits_to_go;
    int buffer;
    int ridx;

    public int input_nbits(byte[] data, int n)
    {
        int c;
        if (bits_to_go < n)
        {
            buffer <<= 8;
            c = (int)(data[ridx++]&0xFF);
            buffer |= c;
            bits_to_go += 8;
        }
        bits_to_go -= n;
        return ( (buffer>>bits_to_go) & ((1<<n)-1) );
    }
    public int input_bit(byte[] data)
    {
        if (bits_to_go == 0)
        {
            buffer = (int)(data[ridx++]&0xFF);
            bits_to_go = 8;
        }
        bits_to_go -= 1;
        return((buffer>>bits_to_go) & 1);
    }

    public void read_bdirect(byte[] data, int[] pixData, int pidx, int n, int nqx, int nqy, int[] scratch, int bit)
    {
        for (int i = 0; i < ((nqx+1)/2) * ((nqy+1)/2); i++)
        {
            scratch[i] = input_nbits(data, 4);
        }
        qtree_bitins(scratch, nqx, nqy, pixData, pidx, n, bit);
    }

    public int input_huffman(byte[] data)
    {
        int c;

        c = input_nbits(data, 3);
        if (c < 4)
        {
           return (1<<c);
        }
        c = input_bit(data) | (c<<1);
        if (c < 13)
        {
            switch (c)
            {
                case  8 : return(3);
                case  9 : return(5);
                case 10 : return(10);
                case 11 : return(12);
                case 12 : return(15);
            }
        }
        c = input_bit(data) | (c<<1);
        if (c < 31)
        {
            switch (c)
            {
                case 26 : return(6);
                case 27 : return(7);
                case 28 : return(9);
                case 29 : return(11);
                case 30 : return(13);
            }
        }
        c = input_bit(data) | (c<<1);
        if (c == 62)
        {
                return(0);
        }
        return(14);
    }

    public void qtree_copy(int[] a, int nx,int ny, int[] b, int n)
    {
        int i, j, k, nx2, ny2;
        int s00, s10;

        nx2 = (nx+1)/2;
        ny2 = (ny+1)/2;
        k = ny2*(nx2-1)+ny2-1;          /* k   is index of a[i,j]           */
        for (i = nx2-1; i >= 0; i--) {
                s00 = 2*(n*i+ny2-1);        /* s00 is index of b[2*i,2*j]       */
                for (j = ny2-1; j >= 0; j--) {
                        b[s00] = a[k];
                        k -= 1;
                        s00 -= 2;
                }
        }
        for (i = 0; i<nx-1; i += 2) {
                s00 = n*i;                  /* s00 is index of b[i,j]   */
                s10 = s00+n;                /* s10 is index of b[i+1,j] */
                for (j = 0; j<ny-1; j += 2) {
                        b[s10+1] =  b[s00]     & 1;
                        b[s10  ] = (b[s00]>>1) & 1;
                        b[s00+1] = (b[s00]>>2) & 1;
                        b[s00  ] = (b[s00]>>3) & 1;
                        s00 += 2;
                        s10 += 2;
                }
                if (j < ny) {
                        b[s10  ] = (b[s00]>>1) & 1;
                        b[s00  ] = (b[s00]>>3) & 1;
                }
        }
        if (i < nx) {
                s00 = n*i;
                for (j = 0; j<ny-1; j += 2) {
                        b[s00+1] = (b[s00]>>2) & 1;
                        b[s00  ] = (b[s00]>>3) & 1;
                        s00 += 2;
                }
                if (j < ny) {
                        b[s00  ] = (b[s00]>>3) & 1;
                }
        }
    }

    public void qtree_expand(byte[] data, int[] a, int nx, int ny, int[] b)
    {
        qtree_copy(a, nx, ny, b, ny);
        for (int i = nx*ny-1; i >= 0; i--)
        {
            if (b[i] != 0) b[i] = input_huffman(data);
        }
    }
    public void qtree_bitins(int[] a, int nx, int ny, int[] b, int pidx, int n, int bit)
    {
            int i, j, k;
            int s00, s10;

            /*
             * expand each 2x2 block
             */
            k = 0;                          /* k   is index of a[i/2,j/2]   */
            for (i = 0; i<nx-1; i += 2) {
                    s00 = n*i;                  /* s00 is index of b[i,j]       */
                    s10 = s00+n;                /* s10 is index of b[i+1,j]     */
                    for (j = 0; j<ny-1; j += 2) {
                            b[s10+1+pidx] |= ( a[k]     & 1) << bit;
                            b[s10  +pidx] |= ((a[k]>>1) & 1) << bit;
                            b[s00+1+pidx] |= ((a[k]>>2) & 1) << bit;
                            b[s00  +pidx] |= ((a[k]>>3) & 1) << bit;
                            s00 += 2;
                            s10 += 2;
                            k += 1;
                    }
                    if (j < ny) {
                            /*
                             * row size is odd, do last element in row
                             * s00+1, s10+1 are off edge
                             */
                            b[s10  +pidx] |= ((a[k]>>1) & 1) << bit;
                            b[s00  +pidx] |= ((a[k]>>3) & 1) << bit;
                            k += 1;
                    }
            }
            if (i < nx) {
                    /*
                     * column size is odd, do last row
                     * s10, s10+1 are off edge
                     */
                    s00 = n*i;
                    for (j = 0; j<ny-1; j += 2) {
                            b[s00+1+pidx] |= ((a[k]>>2) & 1) << bit;
                            b[s00  +pidx] |= ((a[k]>>3) & 1) << bit;
                            s00 += 2;
                            k += 1;
                    }
                    if (j < ny) {
                            /*
                             * both row and column size are odd, do corner element
                             * s00+1, s10, s10+1 are off edge
                             */
                            b[s00  +pidx] |= ((a[k]>>3) & 1) << bit;
                            k += 1;
                    }
            }
    }

  public void qtree_decode(byte[] data, int[] pixData, int pidx, int n, int nqx, int nqy, int nbitplanes)
  {
        int nqmax = (nqx>nqy) ? nqx : nqy;
        int log2n = (int)(Math.log((float) nqmax) / Math.log(2.0) + 0.5);
        if (nqmax > (1<<log2n))
        {
          log2n += 1;
        }
        int nqx2=(nqx+1)/2;
        int nqy2=(nqy+1)/2;
        int[] scratch = new int[nqx2*nqy2];
        for (int bit = nbitplanes-1; bit >= 0; bit--)
        {
                int b = input_nbits(data, 4);
                if (b == 0)
                {
                        read_bdirect(data, pixData, pidx, n, nqx, nqy, scratch, bit);
                } else if (b != 0xf)
                {
                        printMessage("qtree_decode: bad format code" + b);
                        return;
                } else
                {
                        scratch[0] = input_huffman(data);
                        int nx = 1;
                        int ny = 1;
                        int nfx = nqx;
                        int nfy = nqy;
                        int c = 1<<log2n;
                        for (int k = 1; k<log2n; k++)
                        {
                                c = c>>1;
                                nx = nx<<1;
                                ny = ny<<1;
                                if (nfx <= c) { nx -= 1; } else { nfx -= c; }
                                if (nfy <= c) { ny -= 1; } else { nfy -= c; }
                                qtree_expand(data, scratch, nx, ny, scratch);
                        }
                        qtree_bitins(scratch, nqx, nqy, pixData, pidx, n, bit);
                }
        }
  }

public void hinv(int[] a, int nx, int ny, int smooth, int scale)
{
        int nmax, log2n, i, j, k;
        int nxtop,nytop,nxf,nyf,c;
        int oddx,oddy;
        int shift, bit0, bit1, bit2, mask0, mask1, mask2,
                prnd0, prnd1, prnd2, nrnd0, nrnd1, nrnd2, lowbit0, lowbit1;
        int h0, hx, hy, hc;
        int s10, s00;
        nmax = (nx>ny) ? nx : ny;
        log2n = (int)(Math.log((float) nmax) / Math.log(2.0) + 0.5);
        if ( nmax > (1<<log2n) )
        {
                log2n += 1;
        }
        int[] tmp = new int[(nmax+1)/2];
        shift  = 1;
        bit0   = 1 << (log2n - 1);
        bit1   = bit0 << 1;
        bit2   = bit0 << 2;
        mask0  = -bit0;
        mask1  = mask0 << 1;
        mask2  = mask0 << 2;
        prnd0  = bit0 >> 1;
        prnd1  = bit1 >> 1;
        prnd2  = bit2 >> 1;
        nrnd0  = prnd0 - 1;
        nrnd1  = prnd1 - 1;
        nrnd2  = prnd2 - 1;
        a[0] = (a[0] + ((a[0] >= 0) ? prnd2 : nrnd2)) & mask2;
        nxtop = 1;
        nytop = 1;
        nxf = nx;
        nyf = ny;
        c = 1<<log2n;
        for (k = log2n-1; k>=0; k--) {
                c = c>>1;
                nxtop = nxtop<<1;
                nytop = nytop<<1;
                if (nxf <= c) { nxtop -= 1; } else { nxf -= c; }
                if (nyf <= c) { nytop -= 1; } else { nyf -= c; }
                if (k == 0) {
                        nrnd0 = 0;
                        shift = 2;
                }
                for (i = 0; i<nxtop; i++) {
                        unshuffle(a, ny*i,nytop,1,tmp);
                }
                for (j = 0; j<nytop; j++) {
                        unshuffle(a, j,nxtop,ny,tmp);
                }
//                if (smooth) hsmooth(a,nxtop,nytop,ny,scale);
                oddx = nxtop % 2;
                oddy = nytop % 2;
                for (i = 0; i<nxtop-oddx; i += 2) {
                        s00 = ny*i;             /* s00 is index of a[i,j]   */
                        s10 = s00+ny;           /* s10 is index of a[i+1,j] */
                        for (j = 0; j<nytop-oddy; j += 2) {
                                h0 = a[s00  ];
                                hx = a[s10  ];
                                hy = a[s00+1];
                                hc = a[s10+1];
                                hx = (hx + ((hx >= 0) ? prnd1 : nrnd1)) & mask1;
                                hy = (hy + ((hy >= 0) ? prnd1 : nrnd1)) & mask1;
                                hc = (hc + ((hc >= 0) ? prnd0 : nrnd0)) & mask0;
                                lowbit0 = hc & bit0;
                                hx = (hx >= 0) ? (hx - lowbit0) : (hx + lowbit0);
                                hy = (hy >= 0) ? (hy - lowbit0) : (hy + lowbit0);
                                lowbit1 = (hc ^ hx ^ hy) & bit1;
                                h0 = (h0 >= 0)
                                        ? (h0 + lowbit0 - lowbit1)
                                        : (h0 + ((lowbit0 == 0) ? lowbit1 : (lowbit0-lowbit1)));
                                a[s10+1] = (h0 + hx + hy + hc) >> shift;
                                a[s10  ] = (h0 + hx - hy - hc) >> shift;
                                a[s00+1] = (h0 - hx + hy - hc) >> shift;
                                a[s00  ] = (h0 - hx - hy + hc) >> shift;
                                s00 += 2;
                                s10 += 2;
                        }
                        if (oddy!=0) {
                                h0 = a[s00  ];
                                hx = a[s10  ];
                                hx = ((hx >= 0) ? (hx+prnd1) : (hx+nrnd1)) & mask1;
                                lowbit1 = hx & bit1;
                                h0 = (h0 >= 0) ? (h0 - lowbit1) : (h0 + lowbit1);
                                a[s10  ] = (h0 + hx) >> shift;
                                a[s00  ] = (h0 - hx) >> shift;
                        }
                }
                if (oddx!=0) {
                        s00 = ny*i;
                        for (j = 0; j<nytop-oddy; j += 2) {
                                h0 = a[s00  ];
                                hy = a[s00+1];
                                hy = ((hy >= 0) ? (hy+prnd1) : (hy+nrnd1)) & mask1;
                                lowbit1 = hy & bit1;
                                h0 = (h0 >= 0) ? (h0 - lowbit1) : (h0 + lowbit1);
                                a[s00+1] = (h0 + hy) >> shift;
                                a[s00  ] = (h0 - hy) >> shift;
                                s00 += 2;
                        }
                        if (oddy!=0) {
                                /*
                                 * do corner element if both row and column lengths are odd
                                 * s00+1, s10, s10+1 are off edge
                                 */
                                h0 = a[s00  ];
                                a[s00  ] = h0 >> shift;
                        }
                }
                /*
                 * divide all the masks and rounding values by 2
                 */
                bit2 = bit1;
                bit1 = bit0;
                bit0 = bit0 >> 1;
                mask1 = mask0;
                mask0 = mask0 >> 1;
                prnd1 = prnd0;
                prnd0 = prnd0 >> 1;
                nrnd1 = nrnd0;
                nrnd0 = prnd0 - 1;
        }
}

public void unshuffle(int[] a, int aidx, int n, int n2, int[] tmp)
{
        int i;
        int nhalf;
        int p1, p2, pt;

        /*
         * copy 2nd half of array to tmp
         */
        nhalf = (n+1)>>1;
        pt = 0;
        p1 = n2*nhalf+aidx;//&a[n2*nhalf];              /* pointer to a[i]          */
        for (i=nhalf; i<n; i++) {
                tmp[pt] = a[p1];
                //*pt = *p1;
                p1 += n2;
                pt += 1;
        }
        /*
         * distribute 1st half of array to even elements
         */
        p2 = n2*(nhalf-1)+aidx;//&a[ n2*(nhalf-1) ];        /* pointer to a[i]          */
        p1 = ((n2*(nhalf-1))<<1)+aidx;//&a[(n2*(nhalf-1))<<1];      /* pointer to a[2*i]        */
        for (i=nhalf-1; i >= 0; i--) {
                a[p1] = a[p2];
                //*p1 = *p2;
                p2 -= n2;
                p1 -= (n2+n2);
        }
        /*
         * now distribute 2nd half of array (in tmp) to odd elements
         */
        pt = 0;//tmp;
        p1 = n2+aidx;//&a[n2];                  /* pointer to a[i]          */
        for (i=1; i<n; i += 2) {
               // *p1 = *pt;
               a[p1] = tmp[pt];
                p1 += (n2+n2);
                pt += 1;
        }
}

  // decompress .arch.H
   public FITSImage decompress(byte[] data) throws Exception
   {
        ridx = 0;
        char[] line = new char[80];

        FITSImage img = new FITSImage();

        // reading header
        int headSz = 0;
        for (;;)
        {
            headSz++;
            for (int i=0;i<80;i++) line[i] = (char)data[ridx++];
            String s = new String(line);
            if (s.substring(0, 4).equals("END ")) break;
        }
        img.header = new String[headSz];
        ridx = 0;
        headSz = 0;
        for (;;)
        {
            for (int i=0;i<80;i++) line[i] = (char)data[ridx++];
            String s = new String(line);
            img.header[headSz++] = s;
            if (s.substring(0, 4).equals("END ")) break;
        }

        // check magic code
        int magic1 = (int)(data[ridx++]&0xFF);
        int magic2 = (int)(data[ridx++]&0xFF);
        if (magic1!=0xDD && magic2!=0x99)
        {
            printMessage("Bad magic code!");
            return null;
        }

        int nx,ny,scale;
        nx = ((data[ridx++] & 0xff) << 24) | ((data[ridx++] & 0xff) << 16) | ((data[ridx++] & 0xff) << 8) | (data[ridx++] & 0xff);
        ny = ((data[ridx++] & 0xff) << 24) | ((data[ridx++] & 0xff) << 16) | ((data[ridx++] & 0xff) << 8) | (data[ridx++] & 0xff);
        scale = ((data[ridx++] & 0xff) << 24) | ((data[ridx++] & 0xff) << 16) | ((data[ridx++] & 0xff) << 8) | (data[ridx++] & 0xff);
      //  System.out.println(nx + " " + ny + " " + scale);
        int nel = nx*ny;

        img.H = nx;
        img.W = ny;
        img.pixData = new int[nel];
        int sumOfAll = ((data[ridx++] & 0xff) << 24) | ((data[ridx++] & 0xff) << 16) | ((data[ridx++] & 0xff) << 8) | (data[ridx++] & 0xff);
        int[] nbitplanes = new int[3];
        for (int i=0;i<3;i++)
        {
            nbitplanes[i] = (int)(data[ridx++]&0xff);
        }

        int nx2 = (nx+1)/2;
        int ny2 = (ny+1)/2;
        for (int i=0;i<nel;i++)
        {
            img.pixData[i] = 0;
        }
        bits_to_go = 0;
        buffer = 0;
        qtree_decode(data, img.pixData, 0,          ny, nx2,  ny2,  nbitplanes[0]);
        qtree_decode(data, img.pixData, ny2,        ny, nx2,  ny/2, nbitplanes[1]);
        qtree_decode(data, img.pixData, ny*nx2,     ny, nx/2, ny2,  nbitplanes[1]);
        qtree_decode(data, img.pixData, ny*nx2+ny2, ny, nx/2, ny/2, nbitplanes[2]);

        if (input_nbits(data, 4) != 0)
        {
            printMessage("dodecode: bad bit plane values");
        }
        bits_to_go = 0;
        buffer = 0;
        for (int i=0; i<nel; i++)
        {
            if (img.pixData[i] != 0)
            {
               if (input_bit(data) != 0)
               {
                   img.pixData[i] = -img.pixData[i];
               }
            }
        }
        img.pixData[0] = sumOfAll;

        if (scale>1)
        {
            for (int i=0;i<nel;i++)
            {
                img.pixData[i] *= scale;
            }
        }

        hinv(img.pixData, nx, ny, 0, scale);
        for (int i=0;i<nel;i++)
        {
            img.pixData[i] += 32768; // P = BZERO + BSCALE*P;  // assume BZERO = 32768 and BSCALE = 1
        }

        return img;
   }

   public FITSImage extract(byte[] data) throws Exception
   {
        ridx = 0;
        char[] line = new char[80];

        FITSImage img = new FITSImage();

        // reading header
        int headSz = 0;
        for (;;)
        {
            headSz++;
            for (int i=0;i<80;i++) line[i] = (char)data[ridx++];
            String s = new String(line);
            if (s.substring(0, 4).equals("END ")) break;
        }
        img.header = new String[headSz];
        ridx = 0;
        headSz = 0;
        for (;;)
        {
            for (int i=0;i<80;i++) line[i] = (char)data[ridx++];
            String s = new String(line);
            img.header[headSz++] = s;
            if (s.substring(0, 4).equals("END ")) break;
        }

        int nx=0,ny=0,scale=1,zero=0;
        // extract WCS values from header
        for (String row : img.header) {
            if (row.contains("NAXIS1")) {
                String[] tok = row.split("=");
                String[] tok2 = tok[1].split("/");
                ny = Integer.parseInt(tok2[0].trim());
            }
            if (row.contains("NAXIS2")) {
                String[] tok = row.split("=");
                String[] tok2 = tok[1].split("/");
                nx = Integer.parseInt(tok2[0].trim());
            }
            if (row.contains("BSCALE")) {
                String[] tok = row.split("=");
                if (tok[0].contains("BSCALE"))
                {
                    String[] tok2 = tok[1].split("/");
                    scale = (int)Double.parseDouble(tok2[0].trim());
                }
            }
            if (row.contains("BZERO")) {
                String[] tok = row.split("=");
                if (tok[0].contains("BZERO"))
                {
                    String[] tok2 = tok[1].split("/");
                    zero = (int)Double.parseDouble(tok2[0].trim());
                }
            }
        }
        int nel = nx*ny;
        img.H = nx;
        img.W = ny;
        img.pixData = new int[nel];
        // eat spaces, always 640 of them ? (80*8)
        for (int i=0;i<10000;i++)
        {
            int v = data[ridx++]&0xff;
            if (v!=32)
            {
                ridx--;
                break;
            }
        }
        for (int i=0;i<nel;i++)
        {
            // Only 16 bit supported for now
            int v = ((data[ridx++]&0xff) << 8) | (data[ridx++]&0xff);
            v -= zero;
            img.pixData[i] = v*scale;
        }

        return img;
   }

    public FITSImage readArchFits(String archFitsFile) throws Exception
    {
        printMessage("Reading FITS from [" + archFitsFile +"]");

        RandomAccessFile f = new RandomAccessFile(archFitsFile, "r");
        int length = (int)f.length();
        byte[] data = new byte[length];
        f.readFully(data);

        //printMessage("Length = " + length);
        FITSImage img = null;
        if (archFitsFile.substring(archFitsFile.length()-4).toLowerCase().equals("fits") ||
            archFitsFile.substring(archFitsFile.length()-3).toLowerCase().equals("fit"))
            img = extract(data);
        else
            img = decompress(data);

        // extract WCS values from header
        for (String row : img.header) {
            // System.out.println("The header: " + row);
            for (int i=0;i<WCSfields.length;i++) {
                if (row.contains(WCSfields[i])) {
                    String[] tok = row.split("=");
                    String[] tok2 = tok[1].split("/");
                    img.WCS[i] = Double.parseDouble(tok2[0]);
                }
            }
            if (row.contains("TIME-OBS")) {
                String[] tok = row.split("=");
                String[] tok2 = tok[1].split("/");
                img.setTimeOBS(trimString(tok2[0]));
            }
            if (row.contains("DATE-OBS")) {
                String[] tok = row.split("=");
                String[] tok2 = tok[1].split("/");
                img.setDateOBS(trimString(tok2[0]));
            }
            if (row.contains("RA")) {
                String[] tok = row.split("=");
                if (tok.length >= 2 && tok[0].trim().equals("RA")) {
                    String[] tok2 = tok[1].split("/");
                    img.setRA(trimString(tok2[0]));
                }
            }
            if (row.contains("DEC")) {
                String[] tok = row.split("=");
                if (tok.length >= 2 && tok[0].trim().equals("DEC")) {
                    String[] tok2 = tok[1].split("/");
                    img.setDEC(trimString(tok2[0]));
                }
            }
            if (row.contains("DATE")) {
                String[] tok = row.split("=");
                if (tok.length >= 2 && tok[0].trim().equals("DATE")) {
                    String[] tok2 = tok[1].split("/ ");
                    img.setImageCreateTime(trimString(tok2[0]));
                }
            }
        }
        f.close();
        
        return img;
    }
    
    public String trimString(String str) {
        return str.trim().replace('\'', ' ').trim();
    }

    public double scoreAnswer(ArrayList<Entity> userAns, ArrayList<Entity> modelAns, int numOfNEOs) throws Exception
    {
        boolean[] matched = new boolean[modelAns.size()];
        for (int i=0;i<matched.length;i++) matched[i] = false;
        double score = 0;
        double detected = 0;
        int neo_count = 0;
        double neo_detected = 0;
        for (int i=0;i<userAns.size();i++) {
            Entity userA = userAns.get(i);
            userA.matched = false;
            if (userA.isNeo) neo_count++;
            for (int j=0;j<modelAns.size();j++) {
                if (!matched[j] && userA.ID.equals(modelAns.get(j).ID)) {
                    Entity modelA = modelAns.get(j);
                    double sum = 0;
                    for (int f=0;f<4;f++) {
                        sum += (userA.RA[f] - modelA.RA[f]) * (userA.RA[f] - modelA.RA[f]);
                        sum += (userA.DEC[f] - modelA.DEC[f]) * (userA.DEC[f] - modelA.DEC[f]);
                    }
                    if (sum < MATCH_DISTANCE) {
                        userA.matched = true;
                        matched[j] = true;
                        detected += 1.0;
                        score += (1000000.0 / modelAns.size()) * (detected / (i + 1));;
                        if (modelA.isNeo && userA.isNeo) {
                            neo_detected += 1.0;
                            score += (100000.0 / numOfNEOs) * (neo_detected / neo_count);
                        }
                        break;
                    }

                }
            }
        }
        return score;
    }
    private float currentProgress = 0;

    private long imageId;
    
    public void visualizeSet(ImageSet imgset, int frame, String fileName, ArrayList<Entity> userAnsViz) throws Exception
    {
            int W = imgset.img[frame].W;
            int H = imgset.img[frame].H;
            BufferedImage bi = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
            // Graphics2D g = (Graphics2D)bi.getGraphics();
            for (int y=0;y<H;y++)
                for (int x=0;x<W;x++)
                    bi.setRGB(x, y, 0xffffff);
            
            int MARGIN = 5;
            int B16 = 1<<16;
            int WH = W*H;
            int[] histo = new int[B16];
            for(int i=0; i<B16; i++)
                    histo[i] = 0;
            for(int i=0; i<WH; i++) {
                    int x = i%W;
                    if(x<MARGIN || x>=W-MARGIN)
                            continue;
                    int p = imgset.img[frame].pixData[i];
                    p = Math.max(0, Math.min(B16-1, p) );
                    ++histo[p];
            }
            for(int i=1; i<B16; ++i)
                    histo[i] += histo[i-1];
            int center = 0;
            int cmax = 0;
            int spread = 1000;
            for(int i=spread+1; i<B16; i++) {
                    int c = histo[i]-histo[i-spread-1];
                    if(c>cmax) {
                            cmax = c;
                            center = i-spread/2;
                    }
            }
            int dmin = Math.max(0, center-spread/2);
            int dmax = Math.min(B16, center+spread/2);
            for (int y=0;y<H;y++)
            for (int x=0;x<W;x++)
            {
                int ival = imgset.img[frame].pixData[x+y*W];
                ival = Math.max(0, Math.min( 255, (255*(ival-dmin))/(dmax-dmin) ) );
                if (ival<0) ival = 0;
                if (ival>255) ival = 255;
                int cr = ival;
                int cg = ival;
                int cb = ival;
                int rgb = cr + (cg<<8) + (cb<<16);
                bi.setRGB(x,y,rgb);
            }
            
            String visualizationImage = fileName + "_VIZ.png";
            ImageIO.write(bi, "jpeg", new File(visualizationImage));
            bi.flush();
            bi = null;
            for (Entity ent : userAnsViz) {
                double[] XY = WcsConvert.convertRADEC2XY(imgset.img[0], ent.RA[frame], ent.DEC[frame]);
                ent.x[frame] = (int) XY[0];
                ent.y[frame] = (int) XY[1];
                ent.visualizationImages[frame] = visualizationImage;
            }
    }

    public boolean isDetectionInImage(Entity ent, int W, int H) {
     return (ent.x[0] >= 0 && ent.x[0] < W && ent.y[0] >= 0 && ent.y[0] < H &&
             ent.x[1] >= 0 && ent.x[1] < W && ent.y[1] >= 0 && ent.y[1] < H &&
             ent.x[2] >= 0 && ent.x[2] < W && ent.y[2] >= 0 && ent.y[2] < H &&
             ent.x[3] >= 0 && ent.x[3] < W && ent.y[3] >= 0 && ent.y[3] < H);
    }
    
    public void readMpcd(String mpcdFile, ImageSet imgset, boolean isNeoFile) throws Exception
    {
        printMessage("Reading from [" + mpcdFile + "]");
        try {
            BufferedReader br = new BufferedReader(new FileReader(mpcdFile));
            boolean bRead = true;
            while (bRead) {
                String s = br.readLine();
                if (s == null) {
                    break;
                }
                s = s.replaceAll("\\s+", " "); // remove duplicate spaces
                Entity dt = new Entity();
                for (int i=0;i<4;i++)
                {
                    if (i!=0)
                    {
                        s = br.readLine();
                        if (s==null)
                        {
                           bRead = false;
                           break;
                        }
                        s = s.replaceAll("\\s+", " "); // remove duplicate spaces
                    }
                    String[] tok = s.split(" ");
                    if (tok.length<5)
                    {
                       bRead = false;
                       break;
                    }
                    double Hours = Double.parseDouble(tok[4]);
                    double Minutes = Double.parseDouble(tok[5]);
                    double Secs = Double.parseDouble(tok[6]);
                    dt.RA[i] = Hours * 15.0 + Minutes / 4.0 + Secs / 240.0;
                    double Degrees = Double.parseDouble(tok[7]);
                    Minutes = Double.parseDouble(tok[8]);
                    Secs = Double.parseDouble(tok[9]);
                    dt.DEC[i] = Degrees + Math.copySign(Minutes / 60.0 + Secs / 3600.0, Degrees);
                    dt.mag[i] = Double.parseDouble(tok[10]);
                    System.out.println("MAG in Mpcd: " + dt.mag[i]);
                    double[] XY = WcsConvert.convertRADEC2XY(imgset.img[i], dt.RA[i], dt.DEC[i]);
                    dt.x[i] = (int)XY[0];
                    dt.y[i] = (int)XY[1];
//                    double[] INVXY = WcsConvert.convertXY2RADEC(imgset.img[i], dt.x[i], dt.y[i]);
  //                  printMessage(dt.RA[i] + " " + dt.DEC[i] + " " + INVXY[0] + " " + INVXY[1]);
                    dt.isNeo = isNeoFile;
                }
                if (bRead)
                {
                    if (isDetectionInImage(dt, imgset.img[0].W, imgset.img[0].H)) {
                       imgset.detections.add(dt);
                       if (isNeoFile) imgset.numNeos++;
                    }
                }
            }
            br.close();
        } catch (Exception e) {
          //  System.out.println("FAILURE: " + e.getMessage());
        }
     }

    public void readEphm(String ephmFile, ImageSet imgset) throws Exception
    {
      //  String[] ephmData = null;
        printMessage("Reading from [" + ephmFile + "]");
        try {
            BufferedReader br = new BufferedReader(new FileReader(ephmFile));
          //  ArrayList<String> fileData = new ArrayList<String>();
            boolean bRead = true;
            while (bRead) {
                String s = br.readLine();
                if (s == null) {
                    break;
                }
                s = s.replaceAll("\\s+", " "); // remove duplicate spaces
              //  fileData.add(s);

                Entity dt = new Entity();
                String name = "";
                for (int i=0;i<4;i++)
                {
                    if (i!=0)
                    {
                        s = br.readLine();
                        if (s==null)
                        {
                           bRead = false;
                           break;
                        }
                        s = s.replaceAll("\\s+", " "); // remove duplicate spaces
                       // fileData.add(s);
                    }
                    String[] tok = s.split(" ");
                    if (tok.length<5)
                    {
                       bRead = false;
                       break;
                    }
                    dt.RA[i] = Double.parseDouble(tok[1]);// / 15.0;
                    dt.DEC[i] = Double.parseDouble(tok[2]);
                    dt.mag[i] = Double.parseDouble(tok[3]);
                    System.out.println("MAG in Ephm: " + dt.mag[i]);
                    double[] XY = WcsConvert.convertRADEC2XY(imgset.img[i], dt.RA[i], dt.DEC[i]);
                    dt.x[i] = (int)XY[0];
                    dt.y[i] = (int)XY[1];
                    if (i==3)
                    {
                        String[] tokn = s.split("\"");
                        if (tokn.length>3)
                        {
                            name = tokn[3];
                        }
                    }
                   // printMessage(dt.RA[i] + " " + dt.DEC[i] + " " + dt.x[i] + " " + dt.y[i]);
                }
                if (bRead)
                {
                    if (isDetectionInImage(dt, imgset.img[0].W, imgset.img[0].H)) {
                       // search in known NEO list
                       String[] tok = name.split(" ");
                       if (tok.length>1) {
                            name = "";
                            if (tok[0].length()==4 && isNumeric(tok[0]))
                                name = tok[0] + " " + tok[1];
                            if (tok[1].length()==4 && isNumeric(tok[1]) && tok.length>2)
                                name = tok[1] + " " + tok[2];
                            dt.isNeo = (Arrays.binarySearch(knownNEOS, name) >= 0);
                            if (dt.isNeo) {
                                imgset.numNeos++;
                            }
                       }

                       imgset.detections.add(dt);
                    }
                }
            }
            br.close();
            //ephmData = new String[fileData.size()];
            //fileData.toArray(ephmData);
        } catch (Exception e) {
           // System.out.println("FAILURE: " + e.getMessage());
        }
     }

    public ImageSet loadDataSet(String folder, String name, String extension) throws Exception {
        if (extension == null) {
            extension = ".arch.H";
        }
        String basename = folder + name;
        ImageSet imgset = new ImageSet();
        imgset.numNeos = 0;
        imgset.name = name;
        // read arch.H files
        float progress = 0;
        for (int i=0;i<4;i++)
        {
            String fileName = basename + "_000" + (i+1) + extension;
            imgset.img[i] = readArchFits(fileName);
            imgset.img[i].setObservatoryCode(observatoryCode);
            progress += 0.2 / 4;
            this.detectionProgressListener.updateProgress(progress);
        }
        readEphm(basename + "_0001.ephm", imgset);
        readMpcd(basename + "_0001.mpcd", imgset, false);
        readMpcd(basename + "_0001.neos", imgset, true);

        return imgset;

    }

    public boolean isNumeric(String input) {
      try {
        Integer.parseInt(input);
        return true;
      }
      catch (NumberFormatException e) {
        // s is not numeric
        return false;
      }
    }

    public void addKnownObjects(String fileName) throws Exception {

        try {
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            boolean bRead = true;
            while (bRead) {
                String s = br.readLine();
                if (s == null) {
                    break;
                }
                s = s.replaceAll("\\s+", " "); // remove duplicate spaces
                String[] tok = s.split(" ");
                if (tok.length>3)
                {
                    if (tok[1].length()==4 && isNumeric(tok[1]))
                        knownNEOSList.add(tok[1] + " " + tok[2]);
                    if (tok[2].length()==4 && isNumeric(tok[2]))
                        knownNEOSList.add(tok[2] + " " + tok[3]);
                }
            }
            br.close();
        } catch (Exception e) {
           System.out.println("FAILURE: " + e.getMessage());
        }

    }

    public void loadKnownNEOS() throws Exception {
        String knownNEOsPath = trainDataFolder;
        addKnownObjects(knownNEOsPath +  "Amors.txt");
        addKnownObjects(knownNEOsPath + "Apollos.txt");
        addKnownObjects(knownNEOsPath + "Atens.txt");
        knownNEOS = new String[knownNEOSList.size()];
        knownNEOSList.toArray(knownNEOS);
        Arrays.sort(knownNEOS);
    }

    private List<String> splitCommandLine(String str) {
        List<String> result = new ArrayList<String>();
        boolean inQuote = false;
        StringBuffer word = new StringBuffer();
        for (int i = 0; i < str.length(); i++) {
            char c = str.charAt(i);
            boolean endWord = false;
            if (c == '\"') {
                inQuote = !inQuote;
                if (!inQuote) {
                    endWord = true;
                }
            } else if (c == ' ' && !inQuote) {
                endWord = true;
            } else {
                word.append(c);
            }
            if (endWord && word.length() > 0) {
                result.add(word.toString());
                word = new StringBuffer();
            }
        }
        if (word.length() > 0) {
            result.add(word.toString());
        }
        for (String s : result) {
            System.out.println("The Command: " + s);
        }
        return result;
    }
    
    public List<Entity> doExec() throws Exception {
        // launch solution
    	System.out.println("OS TYPE: " + System.getProperty("os.name"));
        printMessage("Executing your solution: " + execCommand + ".");
        
        Process solution = Runtime.getRuntime().exec(splitCommandLine(execCommand).toArray(new String[0]));

        BufferedReader reader = new BufferedReader(new InputStreamReader(solution.getInputStream()));
        PrintWriter writer = new PrintWriter(solution.getOutputStream());
        // PrintWriter writer = new PrintWriter(new FileOutputStream(new File("c:\\tc\\AST-front\\fits_output.txt")));
        new ErrorStreamRedirector(solution.getErrorStream()).start();

        Random rnd = null;
        try {
            rnd = SecureRandom.getInstance("SHA1PRNG");
        } catch (Exception e) {
            System.err.println("ERROR: unable to generate test case.");
            System.exit(1);
        }

        rnd.setSeed(seed);

        loadKnownNEOS();

        // read training file
        
        if (trainFile!=null)
        {
            printMessage("Load training data.");
            BufferedReader br = new BufferedReader(new FileReader(trainFile));
            while (true) {
                String s = br.readLine();
                if (s == null) {
                    break;
                }
                // load set data
                printMessage("load set:" + trainDataFolder + s);
                ImageSet imgset = loadDataSet(trainDataFolder, s, fileExtension);
                printMessage("Number of detections = "+imgset.detections.size());
                printMessage("Number of neos = "+imgset.numNeos);

                // call trainingData(..)
                  writer.println(imgset.img[0].W);
                  writer.println(imgset.img[0].H);
                  for (int frame=0;frame<4;frame++) {
                    // image data
                    for (int pix : imgset.img[frame].pixData) {
                        writer.println(pix);
                    }
                    // header data
                    writer.println(imgset.img[frame].header.length);
                    for (String row : imgset.img[frame].header) {
                        writer.println(row);
                    }
                    // WCS data
                    for (double d : imgset.img[frame].WCS) {
                        writer.println(d);
                    }
                    writer.flush();
                  }
                  // detections
                  String[] strDetect = new String[imgset.detections.size()*4];
                  int dIdx = 0;
                  int lIdx = 0;
                  for (Entity ent : imgset.detections) {
                    for (int frame=0;frame<4;frame++) {
                        strDetect[lIdx] = dIdx + " " + (frame+1) + " " + ent.RA[frame] + " " + ent.DEC[frame] +
                                          " " + ent.x[frame] + " " + ent.y[frame] + " " + ent.mag[frame] + " " + (ent.isNeo ? "1" : "0");
                        lIdx++;
                    }
                    dIdx++;
                  }
                  writer.println(strDetect.length);
                  for (String row : strDetect) {
                      writer.println(row);
                  }
                  writer.flush();

               
                String trainResp = reader.readLine();
                while (trainResp.startsWith("PROGRESS:")) {
                    trainResp = reader.readLine();
                }

                int iResp = Integer.parseInt(trainResp);
                if (iResp == 1)
                    break;
                
            }
            
            // write the end of the training data (separate with the test data)
            writer.println("================");
            writer.flush();
            
            br.close();
        } else
        {
            printMessage("Skipping training phase");
        }
        
        if (testFile == null) {
            // train only
            System.out.println("This execution is for training only.");
            writer.close();
            return null;
        }
        
        // read testing file
        printMessage("Load testing data.");
        String[] testing_sets = new String[NUM_OF_TESTS];
        ArrayList<String> all_sets = new ArrayList<String>();
        {
            BufferedReader br = new BufferedReader(new FileReader(testFile));
            while (true) {
                String s = br.readLine();
                //printMessage(s);
                if (s == null) {
                    break;
                }
                all_sets.add(s);
            }
            br.close();
        }

        // select test images
        boolean[] selected = new boolean[all_sets.size()];
        for (int i=0;i<selected.length;i++) selected[i] = false;
        for (int i=0;i<testing_sets.length;i++) {
            int idx;
            do { 
                idx = rnd.nextInt(all_sets.size());
            } while (selected[idx]);
            testing_sets[i] = all_sets.get(idx);
            selected[idx] = true;
            //printMessage(testing_sets[i]);
        }

        ImageSet[] visSets = null;
        if (visualize)
        {
            visSets = new ImageSet[testing_sets.length];
        }

        ArrayList<Entity> modelAns = new ArrayList<Entity>();
        // int totalNumNeos = 0;
        ImageSet imgset = null;
        this.detectionProgressListener.updateProgress(0);
        Thread progressThread = null;
        for (int i=0;i<testing_sets.length;i++) {
            // load set data
            printMessage("load set:" + folder + testing_sets[i]);
            imgset = loadDataSet(folder, testing_sets[i], fileExtension);
            
            currentImageSet = imgset;
            printMessage("Number of detections = "+imgset.detections.size());
            printMessage("Number of neos = "+imgset.numNeos);
            if (this.imageService != null && !forced) {
                // check if the image set is already exists
                if (this.imageService.imageExists(imgset.img[0])) {
                    // throw the already exists exception
                    writer.close();
                    throw new ImageAlreadyExistException("Image Already Exist");
                }
            }
            // read the image create time
            DateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");
            
            try {
                this.imageCreateTime = sdf.parse(imgset.img[0].getImageCreateTime());
            } catch (java.text.ParseException e) {
                sdf = new SimpleDateFormat("yy/MM/dd");
                this.imageCreateTime = sdf.parse(imgset.img[0].getImageCreateTime());
            }
            
            // totalNumNeos += imgset1.numNeos;
            if (visualize)
            {
                if (visnr==i || visnr==-1)
                    visSets[i] = imgset;
                else
                    visSets[i] = null;
            }

             // detections
             for (Entity ent : imgset.detections) {
                 ent.ID = testing_sets[i];
                 modelAns.add(ent);
                 //printMessage(ent.ID + " " + ent.RA[0] + " " + ent.DEC[0] + " " + ent.x[0] + " " + ent.y[0] + " " + ent.mag[0]);
             }

            // start a new thread to read the progress
            final BufferedReader theReader = reader; 
            progressThread = new Thread(new Runnable() {
                @Override
                public void run() {
                    System.out.println("waiting for response...");
                    try {
                        String cmd = theReader.readLine();
                        System.out.println("readed: " + cmd);
                        while (cmd.startsWith("PROGRESS:")) {
                            detectionProgressListener.updateProgress(Float.parseFloat(cmd.substring("PROGRESS:".length())));
                            cmd = theReader.readLine();
                            System.out.println("readed: " + cmd);
                        }
                    } catch (NumberFormatException e) {
                        e.printStackTrace();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            });
            progressThread.start();
            
            // call testingData(..)
            writer.println(testing_sets[i]);
            writer.println(imgset.img[0].W);
            writer.println(imgset.img[0].H);
            for (int frame=0;frame<4;frame++) {
                // image data
                for (int pix : imgset.img[frame].pixData) {
                    writer.println(pix);
                }
                // header data
                writer.println(imgset.img[frame].header.length);
                for (String row : imgset.img[frame].header) {
                    writer.println(row);
                }
                // WCS data
                for (double d : imgset.img[frame].WCS) {
                    writer.println(d);
                }
                writer.flush();
            }
            
        }
        
        progressThread.join();
        writer.close();
        
        // get response from solution
        String cmd = reader.readLine();
        System.out.println("readed: " + cmd);
        while (cmd.startsWith("PROGRESS:")) {
            this.detectionProgressListener.updateProgress(Float.parseFloat(cmd.substring("PROGRESS:".length())));
            cmd = reader.readLine();
            System.out.println("readed: " + cmd);
        }
        
        
        int n = Integer.parseInt(cmd);
        System.out.println("The returned answer num: " + n + " max answer num: " + maxNumberOfAnswer);
        ArrayList<Entity> userAns = new ArrayList<Entity>();
        for (int i=0;i<n && i < maxNumberOfAnswer;i++) {
            String s = reader.readLine();
            String[] tok = s.split(" ");
            if (tok.length!=12) {
                throw new Exception("ERROR: Incomplete detection returned. Expected 12 space delimited fields, only "+tok.length+" received:\n" + s);
            }
            Entity ent = new Entity();
            ent.ID = tok[0];
            for (int frame=0;frame<4;frame++)
            {
                ent.RA[frame] = Double.parseDouble(tok[1+frame*2]);
                ent.DEC[frame] = Double.parseDouble(tok[2+frame*2]);
                ent.latitudes[frame] = imgset.img[frame].latitude;
                ent.longitudes[frame] = imgset.img[frame].longitude;
                ent.mag[frame] = Double.parseDouble(tok[11]);
            }
            ent.isNeo = tok[9].equals("1");
            ent.significance = Double.parseDouble(tok[10]);
            userAns.add(ent);
            if (userAns.size()>maxNumberOfAnswer) {
                throw new Exception("ERROR: Too many detections returned. Maximum allowed is "+maxNumberOfAnswer);
            }
        }

        
        
        ArrayList<Entity> userAnsViz = new ArrayList<Entity>();
        this.currentProgress = 0.8f;
        if (visualize)
        {
            System.out.println("Visualizing....");
            for (int i=0;i<visSets.length;i++)
            if (visSets[i] != null)
            {
                
                for (Entity ent : userAns) {
                    if (ent.ID.equals(testing_sets[i])) {
                        userAnsViz.add(ent);
                    }
                }

                for (int frame=0;frame<4;frame++) {
                    this.detectionProgressListener.updateProgress(this.currentProgress);
                    visualizeSet(visSets[i], frame, folder + testing_sets[i] + "_" + (frame+1), userAnsViz);
                    this.currentProgress += 0.025;
                }
            }
        }

        this.detectionProgressListener.updateProgress(0.9f);
        
        // persist the image info
        if (this.imageService != null) {
            setImageId(this.imageService.addImage(visSets[0].img[0]));
        }
        
        return userAnsViz;
    }

    public static void main(String[] args) throws Exception {
        AsteroidDetectorTester tester = new AsteroidDetectorTester();
        String outputFolder = null;
        String execCommand = null;
        for (int i = 0; i < args.length; i++) {
             if (args[i].equals("-train")) {
                 tester.setTrainFile(args[++i]);
             } else if (args[i].equals("-exec")) {
                 execCommand = args[++i];
             } else if (args[i].equals("-folder")) {
                 String folder = args[++i];
                 tester.setTrainFolder(folder);
             } else if (args[i].equals("-output")) {
                 outputFolder = args[++i];
             }  
         } 
        
        if (outputFolder == null || execCommand == null) {
            throw new Exception("You must set the output folder and the exec command");
        }

         tester.setExecCommand(execCommand + " --folder " + outputFolder);
         
         // run it
         tester.doExec();
     }
    

    public void setMaxNumberOfAnswer(int maxNumberOfAnswer) {
        this.maxNumberOfAnswer = maxNumberOfAnswer;
    }
    
    public void setDebug(boolean debug) {
        this.debug = debug;
    }
    public void setVisualize(boolean visualize) {
        this.visualize = visualize;
    }
    public void setVisnr(int visnr) {
        this.visnr = visnr;
    }
    public void setExecCommand(String execCommand) {
        this.execCommand = execCommand;
    }
    public void setTrainFile(String trainFile) {
        this.trainFile = trainFile;
    }
    public void setTestFile(String testFile) {
        this.testFile = testFile;
    }
    public void setFolder(String folder) {
        this.folder = folder;
    }
    public void setTrainFolder(String trainDataFolder) {
        this.trainDataFolder = trainDataFolder;
    }
    public void setDetectionProgressListener(DetectionProgressListener detectionProgressListener) {
        this.detectionProgressListener = detectionProgressListener;
    }


    public ImageService getImageService() {
        return imageService;
    }
    public void setImageService(ImageService imageService) {
        this.imageService = imageService;
    }

    public void setForced(boolean forced) {
        this.forced = forced;
    }



    public Date getImageCreateTime() {
        return imageCreateTime;
    }
    public void setImageCreateTime(Date imageCreateTime) {
        this.imageCreateTime = imageCreateTime;
    }



    public ImageSet getCurrentImageSet() {
        return currentImageSet;
    }
    public void setCurrentImageSet(ImageSet currentImageSet) {
        this.currentImageSet = currentImageSet;
    }



    public String getObservatoryCode() {
        return observatoryCode;
    }
    public void setObservatoryCode(String observatoryCode) {
        this.observatoryCode = observatoryCode;
    }



    public long getImageId() {
        return imageId;
    }
    public void setImageId(long imageId) {
        this.imageId = imageId;
    }



    public String getFileExtension() {
        return fileExtension;
    }
    public void setFileExtension(String fileExtension) {
        this.fileExtension = fileExtension;
    }



    class ErrorStreamRedirector extends Thread {
        public BufferedReader reader;

        public ErrorStreamRedirector(InputStream is) {
            reader = new BufferedReader(new InputStreamReader(is));
        }

        public void run() {
            while (true) {
                String s;
                try {
                    s = reader.readLine();
                } catch (Exception e) {
                    // e.printStackTrace();
                    return;
                }
                if (s == null) {
                    break;
                }
                System.out.println(s);
            }
        }
    }
}


