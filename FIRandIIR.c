#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#ifndef PI
# define PI	3.14159265358979323846264338327950288
#endif

#define q	12		            /* for 2^12 points */
#define sampleSize	(1<<q)		/* N-point FFT, iFFT */

/* Choose If you want to use FFTs -> Would recommened leaving it on */
#define FFT
#define miscellaneous
#define MODSIG      // Define if you want AM mod signal
#define DEMODSIG    // Define if you want to Demodulate the ouput signal

/* Choose Between FIR and IIR Implementation */
/* Run main.py */
// #define FIR
#define IIR_PASCALS

/* Take this away for IIR if you want to see the output*/
/* Run the script filterCoeffs.py in contnuous mode and stop running main.py*/
// #define PLOTFILTERCOEFFS

#ifdef FIR
    /* Choose only one FIR Window Method Here */
    // #define HAMMING
    // #define HANNING
    #define BLACKMAN

    /* --- Choose a method of implementation --- */
    #define TRANSVERSAL //Use UFFT For this implementation
    // #define FASTCONV
#endif

#ifdef IIR_PASCALS
    /* --- Choose a method --- */
    /* BANDPASS means non-cascaded version will be used n=2n*/
    /* LP_HP means cascaded version with Low pass and High pass will be used n=n */
    #define BANDPASS
    // #define LP_HP
    // #define BANDPASS_CASCADED

    /* --- Choose analog type --- */
    #define BUTTERWORTH
    // #define ELLIPTICAL

    /* --- Choose the order --- */
    /* Remember that a second order is actually a 4th order n=2n */
    #ifdef BANDPASS 
        // #define FIRST_ORDER_DIRECT_BANDPASS
        #define SECOND_ORDER_DIRECT_BANDPASS // Which is the Only setting that can be used on elliptical
        // #define THIRD_ORDER_DIRECT_BANDPASS
        // #define FOURTH_ORDER_DIRECT_BANDPASS
        // #define TENTH_ORDER_DIRECT_BANDPASS // Tenth order is extremely unstable
    #endif

    /* --- Choose the order --- */
    /* Remember that a second order is actually a 4th order n=2n */
    #ifdef BANDPASS_CASCADED
        // #define FOURTH_ORDER_2NDORDER_BANDPASS_CASCADED
        // #define EIGHTH_ORDER_2NDORDER_BANDPASS_CASCADED
        // #define TENTH_ORDER_2NDORDER_BANDPASS_CASCADED
        // #define TWENTIETH_ORDER_4THORDER_BANDPASS_CASCADED //Technically 40th (n=8)
        // #define FOURTHIETH_ORDER_10THORDER_BANDPASS_CASCADED //Tenth Order is extremely unstable
    #endif

    /* --- Choose the order --- */
    #ifdef LP_HP
        // #define SECOND_ORDER_2ND_ORDER_LP_HP_CASCADED 
        // #define FOURTH_ORDER_2ND_ORDER_LP_HP_CASCADED
        // #define FOURTH_ORDER_SAME_2ND_ORDER_LP_HP_CASCADED //special case
        // #define SIXTH_ORDER_LP_HP_CASCADED_SPECIAL //Change n=6 for this one
        // #define EIGHTH_ORDER_2ND_ORDER_LP_HP_CASCADED
    #endif

    /* define n Order of the fiter*/
    /* If using LP_HP n = n */
    /* If using bandpass n = 2n */
    #define n 4                     // When you update this update the following
    #define ROW (int)(n/2 + 1)      // n/2 + 1
    #define COL n+1                 // n + 1
#endif

#ifdef FFT
    /* --- Choose a method --- */
    #define UFFT
    // #define FFTW
    /* --- Enables Persist --- */
    // #define PERSIST
#endif


/* ===================== General Specifications ===================== */
// This was used for the example of 
double cuttoff_bw = 8000;                               // Cutoff Bandwidth
double channel_bw = 9000;                               // Channel Bandwidth
const double fs = 3.2e6;                                // Sampling Frequency
int channel_selection = 105;                            // Channel of interest
double starting_channel_freq = 526500;                  // Start of spectrum

/* Referred to adcOutput because it was assumed that the information would
be received from an ADC of some sort and then the information would be processed */
double adcOutput[sampleSize];                           
double adcout_FFT_mag_double[sampleSize];
double y_FFT_mag_double[sampleSize];
double y_out[sampleSize];

#ifdef UFFT
    /* --- FFT Variables and Functions --- */
    typedef double real;
    typedef struct{real Re; real Im;} cmplx;

    void ifft(cmplx [], int, cmplx []);
    void fft(cmplx [], int, cmplx []);                  
    void dataPrep(double *,cmplx*);                       // Preps the data for fft functions
    void fftNormalize(cmplx *, cmplx*, int );             // Display FFT
    cmplx scratch[sampleSize];                            // Array need for FFT
    cmplx filterFFT[sampleSize];                          // Temporary array to be passed into FFT
#endif

#ifdef FFTW
    /* --- FFTW --- */
    //macros for the real and imag parts
    #define REAL 0
    #define IMAG 1
    void fftw_FFT(fftw_complex *, fftw_complex *);
    void ifftw_FFT(fftw_complex *, fftw_complex *);
    fftw_complex filterFFTfftw[sampleSize];
#endif
/* ===================== FIR Implementations ===================== */
#ifdef FIR
    void firImplemenation(double* , double* , double*, double*);

    /* --- Filter Coefficents Variables and Functions --- */
    void coeffCalc();
    void zeroPad(double *, double*, int);
    void transversalFilterFunction(double*, double* , double* );

    double h_ZeroPad[sampleSize];                         // Zero padded array filters
    // Passband Ripple
    double A_p_min = 1;
    // Minimum stopband attenuation, A_s_min [dB]:
    double A_s_min = 30;
    // Transition Width
    const double delta_f = 1600;
    
    // -- Hanning
    #ifdef HANNING
        // N = round(3.1/delta_fn)
        int N = 6200;
        double h_final[6200];
        // Half array = (N/2)
        double h[3100];    
    #endif
    // -- Hamming
    #ifdef HAMMING
        // N = round(3.3/delta_fn)
        int N = 6600;
        double h_final[6600];
        // Half array = (N/2)
        double h[3300]; 
    #endif
    // -- Blackman
    #ifdef BLACKMAN
        // N = round(5.5/delta_fn) 
        int N = 11000;
        double h_final[11000];
        // Half array = (N/2)
        double h[5500];
    #endif

    double bpfFilter(int, double, double); 
    void filterCoeffCalculator(double, double);
#endif

/* ===================== IIR Implementation ===================== */
#ifdef IIR_PASCALS
    #ifdef BUTTERWORTH
        #ifdef LP_HP
            // ======================= LP and HP 2nd Order cascades ========================== //    
            #ifdef SECOND_ORDER_2ND_ORDER_LP_HP_CASCADED
                double a_i_1[] = {1,0,0};
                double b_i_1[] = {1,1.414214,1};
            #endif

            #ifdef FOURTH_ORDER_SAME_2ND_ORDER_LP_HP_CASCADED
                // 4th order
                double a_i_1[] = {1,0,0};
                double b_i_1[] = {1,1.414214,1};

                double a_i_2[] = {1,0,0};
                double b_i_2[] = {1,1.414214,1};
            #endif


            #ifdef FOURTH_ORDER_2ND_ORDER_LP_HP_CASCADED
                // 4th order
                double a_i_1[] = {1,0,0};
                double b_i_1[] = {1,0.765367,1};

                double a_i_2[] = {1,0,0};
                double b_i_2[] = {1,1.847759,1};
            #endif

            #ifdef EIGHTH_ORDER_2ND_ORDER_LP_HP_CASCADED
                // 8th order
                double a_i_1[] = {1,0,0};
                double b_i_1[] = {1,0.390181,1};

                double a_i_2[] = {1,0,0};
                double b_i_2[] = {1,1.11114,1};

                double a_i_3[] = {1,0,0};
                double b_i_3[] = {1,1.662939,1};

                double a_i_4[] = {1,0,0};
                double b_i_4[] = {1,1.961571,1};
            #endif
        #endif
        // ======================= BANDPASS ========================== //
        #ifdef BANDPASS
            #ifdef FIRST_ORDER_DIRECT_BANDPASS
                // First order butterworth
                // Make n=2 for this for bandpass
                double a_i_bp[] = {1,0};
                double b_i_bp[] = {1,1};
            #endif
            
            #ifdef SECOND_ORDER_DIRECT_BANDPASS
                // Second order butterworth
                // Make n=4 for this for bandpass
                double a_i_bp[] = {1,0,0};
                double b_i_bp[] = {1,1.414214,1};
            #endif

            #ifdef THIRD_ORDER_DIRECT_BANDPASS
                // Third order butterworth
                // Make n=6 for this for bandpass
                double a_i_bp[] = {1,0,0,0};
                double b_i_bp[] = {1,2,2,1};
            #endif

            #ifdef FOURTH_ORDER_DIRECT_BANDPASS
                // Fourth Order ButterWorth Filter 
                // Make n=8 for this for bandpass
                double a_i_bp[] = {1,0,0,0,0};
                double b_i_bp[] = {1,2.613126,3.414214,2.613126,1}; 
            #endif

            #ifdef TENTH_ORDER_DIRECT_BANDPASS
                // 10th Order ButterWorth Filter 
                // Make n=20 for this for bandpass
                // Make n=10 for this for lowpass and highpass
                double a_i_bp[] = {1,0,0,0,0,0,0,0,0,0,0};
                double b_i_bp[] = {1,6.392453,20.431729,42.802061,64.882396,74.233429,64.882396,42.802061,20.431729,6.392453,1}; 
            #endif
        #endif
        // ======================= BANDPASS CASCADED ========================== //
        #ifdef BANDPASS_CASCADED
            #if defined (FOURTH_ORDER_2NDORDER_BANDPASS_CASCADED)|| defined (EIGHTH_ORDER_2NDORDER_BANDPASS_CASCADED)
                // First order butterworth
                // Make n=2 for this for bandpass
                double a_i_bp_c[] = {1,0};
                double b_i_bp_c[] = {1,1};
            #endif

            #ifdef TENTH_ORDER_2NDORDER_BANDPASS_CASCADED
                // First order butterworth
                // Make n=2 for this for bandpass
                double a_i_bp_c[] = {1,0};
                double b_i_bp_c[] = {1,1};
            #endif

            #ifdef TWENTIETH_ORDER_4THORDER_BANDPASS_CASCADED
                // Fourth Order ButterWorth Filter 
                // Make n=8 for this for bandpass
                // Too much attenuation
                double a_i_bp_c[] = {1,0,0,0,0};
                double b_i_bp_c[] = {1,2.613126,3.414214,2.613126,1}; 
            #endif
            
            #ifdef FOURTHIETH_ORDER_10THORDER_BANDPASS_CASCADED
                // 10th Order ButterWorth Filter 
                // Make n=20 for this for bandpass
                // super unstable
                double a_i_bp_c[] = {1,0,0,0,0,0,0,0,0,0,0};
                double b_i_bp_c[] = {1,6.392453,20.431729,42.802061,64.882396,74.233429,64.882396,42.802061,20.431729,6.392453,1}; 
            #endif
        #endif
    #endif

    #ifdef ELLIPTICAL
        // ======================= BANDPASS ========================== //
        // 60dB (Table A.16) - Has a specific form
        #ifdef BANDPASS
            //TODO: ADD
            #ifdef SECOND_ORDER_DIRECT_BANDPASS
                // First order butterworth
                // Make n=2 for this for bandpass
                double a_i_bp[] = {3.16227766e-2, 0, 0.9979372256};
                double b_i_bp[] = {1,1.07759259,1.11970398};
            #endif
        #endif
    #endif

    void iir_coef_bandpass(int ,double* ,double*, double*, double*); 
    void iir_coef_highpass(int ,double* ,double*, double*, double*); 
    void iirImplementation(double*, double* ,double*, double*);
    
#endif
// --- used for IIR_PASCALS and used for DEMODSIG
double gain; 
double factorial(int);
double dot(double *, double *, int);
void iir_coef_lowpass(int ,double* ,double*, double*, double*); 

// These help with generating signals that would be recieved in real world
void signalGenerator(double*, int);
// Used with the persist function to demonstrate the filter response shape
void chirpGenerator(double*, double);

#ifdef DEMODSIG
    /* -- Lowpass setup for signal demodulator -- */
    int demodOrder = 4;
    // Lowpass signal varibles to retrieve the signal
    double f_c_sig_demod = 4000;
    // Filter coefficient length (demodOrder+1)
    int zNumLenDemod = 5;
    int zDenLenDemod = 5;
    // Fourth order lowpass butterworth coefficients
    double a_i_demod[] = {1,0,0,0,0};
    double b_i_demod[] = {1,2.6131,3.4142,2.6131,1}; 

    /* -- Function declartaions to be used -- */
    // Function used for signal demodulation
    void signalDemodulator(double*, double*, double *); 
    int intervalLen = 1; //Viewing window on how far back something should be detected
    // Global variables used for signal demodulation
    double demodSig_double[sampleSize];
    double demodSigFFT_double[sampleSize];
    // Max function used in demod of signal
    double max(double num1, double num2);

#endif

#ifdef miscellaneous
    /* --- Textfile Variables and Functions --- */
    char nameOfFile[27];
    void toTextFile(double*, int, char*);
#endif
void Display();

// ================================================================== //
// ================================================================== //
// ================================================================== //
//                  CAN LEAVE EVERYTHING BELOW HERE                   //
// ================================================================== //
// ================================================================== //
// ================================================================== //
int main() {

#ifdef PERSIST
    double timeAverageArray[sampleSize];
    /* ------ PERSIST ------ */
    for (int persist_num=0; persist_num < sampleSize; persist_num+=1) {
#endif
    /*_______________________*/
    // =============================================================//
    //                      Coefficient Calculator                  //
    // =============================================================//
    // Calculates the Coefficients that will be used in the FIR     //
    // _____________________________________________________________//

    #ifdef IIR_PASCALS
        //Bandpass Only
        double zfilterNum[COL];
        double zfilterDen[COL];

        //BP LP
        // Low Pass
        double zfilterNum1[COL];
        double zfilterDen1[COL];
        double zfilterNum2[COL];
        double zfilterDen2[COL];
        double zfilterNum3[COL];
        double zfilterDen3[COL];
        double zfilterNum4[COL];
        double zfilterDen4[COL];

        // High Pass
        double zfilterNum5[COL];
        double zfilterDen5[COL];
        double zfilterNum6[COL];
        double zfilterDen6[COL];
        double zfilterNum7[COL];
        double zfilterDen7[COL];
        double zfilterNum8[COL];
        double zfilterDen8[COL];
    #endif

    #ifdef FIR
        coeffCalc(channel_selection);    
    #endif

    #ifdef IIR_PASCALS
        /* ===================== BANDPASS ===================== */
        #ifdef BANDPASS
            /* Single Channel Test */
            iir_coef_bandpass(channel_selection, zfilterNum, zfilterDen, a_i_bp, b_i_bp);
            #ifdef PLOTFILTERCOEFFS
                toTextFile(zfilterNum,2,"ZfilterNum1");
                toTextFile(zfilterDen,2,"ZfilterDen1"); 
            #endif
        #endif
        /* ===================== BANDPASS CASCADE ===================== */
        #ifdef BANDPASS_CASCADED
            #ifdef FOURTH_ORDER_2NDORDER_BANDPASS_CASCADED
                /* Single Channel Test */
                iir_coef_bandpass(channel_selection, zfilterNum1, zfilterDen1, a_i_bp_c, b_i_bp_c);
                iir_coef_bandpass(channel_selection, zfilterNum2, zfilterDen2, a_i_bp_c, b_i_bp_c);
                #ifdef PLOTFILTERCOEFFS
                    toTextFile(zfilterNum1,2,"ZfilterNum1");
                    toTextFile(zfilterDen1,2,"ZfilterDen1"); 

                    toTextFile(zfilterNum2,2,"ZfilterNum2");
                    toTextFile(zfilterDen2,2,"ZfilterDen2"); 
                #endif
            #endif

            /* ------------------------- 8th Order or 40th Order ------------------------------ */
            #if defined (EIGHTH_ORDER_2NDORDER_BANDPASS_CASCADED) || defined (FOURTHIETH_ORDER_10THORDER_BANDPASS_CASCADED)
                iir_coef_bandpass(channel_selection,zfilterNum1,zfilterDen1,a_i_bp_c,b_i_bp_c);
                iir_coef_bandpass(channel_selection,zfilterNum2,zfilterDen2,a_i_bp_c,b_i_bp_c);

                iir_coef_bandpass(channel_selection,zfilterNum3,zfilterDen3,a_i_bp_c,b_i_bp_c);
                iir_coef_bandpass(channel_selection,zfilterNum4,zfilterDen4,a_i_bp_c,b_i_bp_c);

                #ifdef PLOTFILTERCOEFFS
                    toTextFile(zfilterNum1,2,"ZfilterNum1");
                    toTextFile(zfilterDen1,2,"ZfilterDen1"); 

                    toTextFile(zfilterNum2,2,"ZfilterNum2");
                    toTextFile(zfilterDen2,2,"ZfilterDen2"); 

                    toTextFile(zfilterNum3,2,"ZfilterNum3");
                    toTextFile(zfilterDen3,2,"ZfilterDen3"); 

                    toTextFile(zfilterNum4,2,"ZfilterNum4");
                    toTextFile(zfilterDen4,2,"ZfilterDen4"); 
                #endif
            #endif
            /* ------------------------- 10th Order or 20th Order ------------------------------ */
            #if defined (TENTH_ORDER_2NDORDER_BANDPASS_CASCADED)||defined (TWENTIETH_ORDER_4THORDER_BANDPASS_CASCADED)
                iir_coef_bandpass(channel_selection,zfilterNum1,zfilterDen1,a_i_bp_c,b_i_bp_c);
                iir_coef_bandpass(channel_selection,zfilterNum2,zfilterDen2,a_i_bp_c,b_i_bp_c);
                iir_coef_bandpass(channel_selection,zfilterNum3,zfilterDen3,a_i_bp_c,b_i_bp_c);
                iir_coef_bandpass(channel_selection,zfilterNum4,zfilterDen4,a_i_bp_c,b_i_bp_c);
                iir_coef_bandpass(channel_selection,zfilterNum5,zfilterDen5,a_i_bp_c,b_i_bp_c);

                #ifdef PLOTFILTERCOEFFS
                    toTextFile(zfilterNum1,2,"ZfilterNum1");
                    toTextFile(zfilterDen1,2,"ZfilterDen1"); 

                    toTextFile(zfilterNum2,2,"ZfilterNum2");
                    toTextFile(zfilterDen2,2,"ZfilterDen2"); 

                    toTextFile(zfilterNum3,2,"ZfilterNum3");
                    toTextFile(zfilterDen3,2,"ZfilterDen3"); 

                    toTextFile(zfilterNum4,2,"ZfilterNum4");
                    toTextFile(zfilterDen4,2,"ZfilterDen4"); 

                    toTextFile(zfilterNum5,2,"ZfilterNum5");
                    toTextFile(zfilterDen5,2,"ZfilterDen5"); 
                #endif
            #endif
        #endif

        /* ===================== LP HP CASCADE ===================== */
        #ifdef LP_HP
            /* ------------------------- 2th Order ------------------------------ */
            #ifdef SECOND_ORDER_2ND_ORDER_LP_HP_CASCADED
                iir_coef_lowpass(channel_selection,zfilterNum1,zfilterDen1,a_i_1,b_i_1);
                iir_coef_highpass(channel_selection,zfilterNum2,zfilterDen2,a_i_1,b_i_1);

                #ifdef PLOTFILTERCOEFFS
                    toTextFile(zfilterNum1,2,"ZfilterNum1");
                    toTextFile(zfilterDen1,2,"ZfilterDen1"); 

                    toTextFile(zfilterNum2,2,"ZfilterNum2");
                    toTextFile(zfilterDen2,2,"ZfilterDen2"); 
                #endif
            #endif
            /* ------------------------- 4th Order ------------------------------ */
            #if defined FOURTH_ORDER_2ND_ORDER_LP_HP_CASCADED || FOURTH_ORDER_SAME_2ND_ORDER_LP_HP_CASCADED
                iir_coef_lowpass(channel_selection,zfilterNum1,zfilterDen1,a_i_1,b_i_1);
                iir_coef_lowpass(channel_selection,zfilterNum2,zfilterDen2,a_i_2,b_i_2);

                iir_coef_highpass(channel_selection,zfilterNum3,zfilterDen3,a_i_1,b_i_1);
                iir_coef_highpass(channel_selection,zfilterNum4,zfilterDen4,a_i_2,b_i_2);

                #ifdef PLOTFILTERCOEFFS
                    toTextFile(zfilterNum1,2,"ZfilterNum1");
                    toTextFile(zfilterDen1,2,"ZfilterDen1"); 

                    toTextFile(zfilterNum2,2,"ZfilterNum2");
                    toTextFile(zfilterDen2,2,"ZfilterDen2"); 

                    toTextFile(zfilterNum3,2,"ZfilterNum3");
                    toTextFile(zfilterDen3,2,"ZfilterDen3"); 

                    toTextFile(zfilterNum4,2,"ZfilterNum4");
                    toTextFile(zfilterDen4,2,"ZfilterDen4"); 
                #endif
            #endif
            /* ------------------------- 6th Order ------------------------------ */
            /* Works with sixth order coefficients, but is a cascaded LP and HP   */
            
            #ifdef SIXTH_ORDER_LP_HP_CASCADED_SPECIAL
                iir_coef_lowpass(channel_selection,zfilterNum1,zfilterDen1,a_i_bp,b_i_bp);
                iir_coef_highpass(channel_selection,zfilterNum2,zfilterDen2,a_i_bp,b_i_bp);
                #ifdef PLOTFILTERCOEFFS
                    toTextFile(zfilterNum1,2,"ZfilterNum1");
                    toTextFile(zfilterDen1,2,"ZfilterDen1"); 

                    toTextFile(zfilterNum2,2,"ZfilterNum2");
                    toTextFile(zfilterDen2,2,"ZfilterDen2"); 
                #endif
            #endif
            /* ------------------------- 8th Order Cascaded Bandpass ------------------------------ */
            #ifdef EIGHTH_ORDER_2ND_ORDER_LP_HP_CASCADED
                // Lowpass Coefficients
                iir_coef_lowpass(channel_selection,zfilterNum1,zfilterDen1,a_i_1,b_i_1);
                iir_coef_lowpass(channel_selection,zfilterNum2,zfilterDen2,a_i_2,b_i_2);
                iir_coef_lowpass(channel_selection,zfilterNum3,zfilterDen3,a_i_3,b_i_3);
                iir_coef_lowpass(channel_selection,zfilterNum4,zfilterDen4,a_i_4,b_i_4);
                // Highpass Coefficients
                iir_coef_highpass(channel_selection,zfilterNum5,zfilterDen5,a_i_1,b_i_1);
                iir_coef_highpass(channel_selection,zfilterNum6,zfilterDen6,a_i_2,b_i_2);
                iir_coef_highpass(channel_selection,zfilterNum7,zfilterDen7,a_i_3,b_i_3);
                iir_coef_highpass(channel_selection,zfilterNum8,zfilterDen8,a_i_4,b_i_4);

                #ifdef PLOTFILTERCOEFFS
                    toTextFile(zfilterNum1,2,"ZfilterNum1");
                    toTextFile(zfilterDen1,2,"ZfilterDen1"); 

                    toTextFile(zfilterNum2,2,"ZfilterNum2");
                    toTextFile(zfilterDen2,2,"ZfilterDen2"); 

                    toTextFile(zfilterNum3,2,"ZfilterNum3");
                    toTextFile(zfilterDen3,2,"ZfilterDen3"); 

                    toTextFile(zfilterNum4,2,"ZfilterNum4");
                    toTextFile(zfilterDen4,2,"ZfilterDen4"); 

                    toTextFile(zfilterNum5,2,"ZfilterNum5");
                    toTextFile(zfilterDen5,2,"ZfilterDen5"); 

                    toTextFile(zfilterNum6,2,"ZfilterNum6");
                    toTextFile(zfilterDen6,2,"ZfilterDen6"); 

                    toTextFile(zfilterNum7,2,"ZfilterNum7");
                    toTextFile(zfilterDen7,2,"ZfilterDen7"); 

                    toTextFile(zfilterNum8,2,"ZfilterNum8");
                    toTextFile(zfilterDen8,2,"ZfilterDen8"); 
                #endif
            #endif
        #endif
    #endif
 
    /* -------------------------------------------------------------*/
    // =============================================================//
    //           ADC Data Needs to be read in here                  //
    // =============================================================//
    /* -------------------------------------------------------------*/
    // Read in ADC Data over here
    // please call the array that comes from ADC = adcOutput

    // =============================================================//
    //                      Signal Generator
    // =============================================================//
    #ifndef PERSIST
        #ifdef MODSIG
            signalGenerator(adcOutput, 2);        
        #endif
        #ifndef MODSIG
            signalGenerator(adcOutput, 1);        
        #endif
    #endif

    #ifdef PERSIST
        float temp = starting_channel_freq+((double)(persist_num)*(fs/2)*(1/(double)sampleSize));
        chirpGenerator(adcOutput, temp);
    #endif 
    // =============================================================//

    clock_t t;
    t = clock();
    // Calculate the time taken by FIR or IIR

    // =============================================================//
    //                      Filter Implementation
    // =============================================================//
    #ifdef FIR
        firImplemenation(adcOutput, adcout_FFT_mag_double, y_FFT_mag_double, y_out);
    #endif

    #ifdef IIR_PASCALS
        /* ===================== BANDPASS ===================== */
        #ifdef BANDPASS
            // Filter the signal //
            iirImplementation(adcOutput, zfilterNum, zfilterDen, y_out);
        #endif
        /* ===================== BANDPASS CASCADE ===================== */
        #ifdef BANDPASS_CASCADED
            /* ------------------------- 4th Order ------------------------------ */
            #ifdef FOURTH_ORDER_2NDORDER_BANDPASS_CASCADED
                double temp1[sampleSize];

                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,y_out);
            #endif
            /* ------------------------- 8th Order ------------------------------ */
            #ifdef EIGHTH_ORDER_2NDORDER_BANDPASS_CASCADED
                double temp1[sampleSize];
                double temp2[sampleSize];
                double temp3[sampleSize];

                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,temp2);
                iirImplementation(temp2,zfilterNum3,zfilterDen3,temp3);
                iirImplementation(temp3,zfilterNum4,zfilterDen4,y_out);
            #endif

            /* ------------------------- 10th Order or 20th Order ------------------------------ */
            #if defined (TENTH_ORDER_2NDORDER_BANDPASS_CASCADED)|| defined (TWENTIETH_ORDER_4THORDER_BANDPASS_CASCADED)
                double temp1[sampleSize];
                double temp2[sampleSize];
                double temp3[sampleSize];
                double temp4[sampleSize];

                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,temp2);
                iirImplementation(temp2,zfilterNum3,zfilterDen3,temp3);
                iirImplementation(temp3,zfilterNum4,zfilterDen4,temp4);
                iirImplementation(temp4,zfilterNum5,zfilterDen5,y_out);
            #endif

            /* ------------------------- 40th Order ------------------------------ */
            /* Super Unstable */
            #ifdef FOURTHIETH_ORDER_10THORDER_BANDPASS_CASCADED   
                double temp1[sampleSize];
                double temp2[sampleSize];
                double temp3[sampleSize];

                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,temp2);
                iirImplementation(temp2,zfilterNum3,zfilterDen3,temp3);
                iirImplementation(temp3,zfilterNum4,zfilterDen4,y_out);   
            #endif
        #endif
        /* ===================== LP HP CASCADE ===================== */
        #ifdef LP_HP
            /* ------------------------- 2th Order ------------------------------ */
            #ifdef SECOND_ORDER_2ND_ORDER_LP_HP_CASCADED
                double temp1[sampleSize];

                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,y_out);
            #endif
            
            /* ------------------------- 4th Order ------------------------------ */
            #ifdef FOURTH_ORDER_2ND_ORDER_LP_HP_CASCADED
                double temp1[sampleSize];
                double temp2[sampleSize];
                double temp3[sampleSize];

                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,temp2);
                iirImplementation(temp2,zfilterNum3,zfilterDen3,temp3);
                iirImplementation(temp3,zfilterNum4,zfilterDen4,y_out);         
            #endif

            /* ------------------------- 4th Order Special Case ------------------------------ */
            // Do not use in general - Would not recommend
            #ifdef FOURTH_ORDER_SAME_2ND_ORDER_LP_HP_CASCADED
                double temp1[sampleSize];
                double temp2[sampleSize];
                double temp3[sampleSize];

                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,temp2);
                iirImplementation(temp2,zfilterNum3,zfilterDen3,temp3);
                iirImplementation(temp3,zfilterNum4,zfilterDen4,y_out);         
            #endif

            /* ------------------------- 6th Order ------------------------------ */
            #ifdef SIXTH_ORDER_LP_HP_CASCADED_SPECIAL
                double temp1[sampleSize];

                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,y_out);
            #endif
            /* ------------------------- 8th Order ------------------------------ */
            #ifdef EIGHTH_ORDER_2ND_ORDER_LP_HP_CASCADED
                double temp1[sampleSize];
                double temp2[sampleSize];
                double temp3[sampleSize];
                double temp4[sampleSize];
                double temp5[sampleSize];
                double temp6[sampleSize];
                double temp7[sampleSize];
                iirImplementation(adcOutput,zfilterNum1,zfilterDen1,temp1);
                iirImplementation(temp1,zfilterNum2,zfilterDen2,temp2);
                iirImplementation(temp2,zfilterNum3,zfilterDen3,temp3);
                iirImplementation(temp3,zfilterNum4,zfilterDen4,temp4);
                iirImplementation(temp4,zfilterNum5,zfilterDen5,temp5);
                iirImplementation(temp5,zfilterNum6,zfilterDen6,temp6);
                iirImplementation(temp6,zfilterNum7,zfilterDen7,temp7);
                iirImplementation(temp7,zfilterNum8,zfilterDen8,y_out);   
            #endif
        #endif
    #endif
    
    // Time
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    #ifdef PERSIST
        timeAverageArray[persist_num]=time_taken;
    #endif
    printf("Implementation took %lf seconds to execute \n", time_taken);

    // =============================================================//
    //                              Display                         // 
    // =============================================================//
    /* At the input */
    //_______________//
    // FFT Magnitude Array of input signal - adcout_FFT_mag_double(double) 
    // Time Array - This will be directly from the ADC

    /* At the output */
    //_______________//
    // FFT Magnitude Array - y_FFT_mag_double(double) 
    // Time Array - y_out(double(We disregard the imaginary in the time domain because it is so small))

    Display();

#ifdef PERSIST
    /* ------ PERSIST ------ */
    int num = persist_num;
    char snum[5];
    snprintf( snum, 5, "%d", persist_num );
    toTextFile(y_FFT_mag_double,0,snum);
};

/* --- Calculate the average time taken of a specific Implementation --- */
double sum = 0;
for (int i=0;i<sampleSize;i++){
    sum += timeAverageArray[i];
}
printf("Average time taken to execute implementation %lf",sum/(double)sampleSize);

#endif

}

// =============================================================//
//                              Display                         // 
// =============================================================//
void Display() {
    #ifdef FFTW
    /* --------------------- ADC data ------------------------- //
        /* FFT on ADC values */
        fftw_complex adcOutputTime[sampleSize];
        fftw_complex adcOutFFT[sampleSize];                          // Temporary array to be passed into FFT
        
        /* DataPrep */
        for (int i=0; i<sampleSize; i++) {
            adcOutputTime[i][REAL] = adcOutput[i];
            adcOutputTime[i][IMAG] = 0;
        }
        
        /* ADC FFT */
        fftw_FFT(adcOutputTime, adcOutFFT);

        /* ADC Display Data */
        for (int i=0;i<sampleSize;i++) {
            adcout_FFT_mag_double[i] = 20*log10( ( sqrt ( pow( adcOutFFT[i][REAL], 2) + pow(adcOutFFT[i][IMAG], 2) ) ) );
        }
        
        fftw_complex filterOutputTime[sampleSize];
        fftw_complex filterOutFFT[sampleSize];                        // Temporary array to be passed into FFT

        /* DataPrep */
        for (int i=0; i<sampleSize; i++) {
            filterOutputTime[i][REAL] = y_out[i];
            filterOutputTime[i][IMAG] = 0;
        }

        /* Filter Output FFT */
        fftw_FFT(filterOutputTime, filterOutFFT);

        /* Gets the FFT magnitude of the output of filter data */     
        for (int i=0;i<sampleSize;i++) {
            y_FFT_mag_double[i] = 20*log10( ( sqrt ( pow( filterOutFFT[i][REAL], 2) + pow(filterOutFFT[i][IMAG], 2) ) ) );
        }

    #endif
    #ifdef UFFT
        /* FFT on ADC values */
        cmplx adcOutFFT[sampleSize];                          // Temporary array to be passed into FFT
        dataPrep(adcOutput, adcOutFFT);
        fft(adcOutFFT, sampleSize, scratch );
        
        /* Gets the FFT magnitude of the adc data */
        cmplx adcMagOut[sampleSize];
        fftNormalize(adcOutFFT,adcMagOut,sampleSize);

        /* Converting the FFT Array to a double */
        for (int i=0;i<sampleSize;i++) adcout_FFT_mag_double[i] = adcMagOut[i].Re;

        /* FFT on Output */
        cmplx fftOut[sampleSize];         
        dataPrep(y_out, fftOut);
        fft(fftOut, sampleSize, scratch);

        /* Gets the FFT magnitude of the adc data */
        cmplx y_FFT_mag[sampleSize];
        fftNormalize(fftOut, y_FFT_mag, sampleSize); // Displays the whole FFT(can only use half the array for display)
        
        /* Converting the FFT Array to a double */
        for (int i=0; i<sampleSize; i++) y_FFT_mag_double[i] = y_FFT_mag[i].Re;
    #endif

    #ifdef DEMODSIG
        signalDemodulator(y_out, demodSig_double, demodSigFFT_double);
    #endif


    #ifdef miscellaneous
        // TODO : DELETE
        /* Save to text file functions */
        #ifndef PERSIST
            toTextFile(y_FFT_mag_double,0,"FFT Filter Output");
            toTextFile(y_out,1,"Filter Ouput");
            toTextFile(adcOutput,1,"ADC Output");
            toTextFile(adcout_FFT_mag_double,0," FFT ADC Output");
            #ifdef DEMODSIG
                toTextFile(demodSig_double,1,"Demodulated Signal");
                toTextFile(demodSigFFT_double,0," FFT demodulated signal");
            #endif
        #endif
    #endif
}


#ifdef miscellaneous
// ==================================================================== //
//                      Send to a textfile in C                         //
// ==================================================================== //
void toTextFile(double* arrToTF, int nameChr, char nameOfFile[]) {
    FILE * fp;

    /* Dynamic Names */
    char buffer[32];
    int sizeOfArrtoTF = 0;
    /* open the file for writing*/
    // Frequency domain file
    if (nameChr == 0) {
        snprintf(buffer,sizeof(char)*32,"f_%s.txt",nameOfFile);
        fp = fopen (buffer,"w");    
        sizeOfArrtoTF = sampleSize;
    }
    // Time domain file
    if (nameChr == 1) {
        snprintf(buffer,sizeof(char)*32,"t_%s.txt",nameOfFile);
        fp = fopen (buffer,"w");    
        sizeOfArrtoTF = sampleSize;

    }
    #ifdef IIR_PASCALS
        if (nameChr == 2) {
            snprintf(buffer,sizeof(char)*32,"z_%s.txt",nameOfFile);
            fp = fopen (buffer,"w");    
            sizeOfArrtoTF = COL;
        }
    #endif

    /* write 10 lines of text into the file stream*/
    for(int i = 0; i < sizeOfArrtoTF;i++){
        fprintf(fp, "%lf\n",arrToTF[i]);
    }

    /* close the file*/  
    fclose (fp);
}
#endif
/*                                                                          */
// These are just functions so that can be placed anywhere after the main loop
/*                                                                          */
// =============================================================//
//                      Signal Generator                        //
// =============================================================//
// This generates a test signal that we can use to test         //
// -------------------------------------------------------------//
void signalGenerator(double* result, int dynamicSigGenerator) {
    double period = 1/fs;

    /* Normal Signal Generator */
    if (dynamicSigGenerator == 1) {
        double freq_sine;
        printf("\nFrequency of Signal:\n");
        scanf("%lf",&freq_sine);
        double sine[sampleSize];
        for(int i = 0; i < sampleSize; i++){
            sine[i] = 1*sin(2.0*PI*freq_sine*(period)*i);
            result[i] = sine[i]; 
        }   
    } 
    /* Modulated Signal Generator */
    else if (dynamicSigGenerator == 2) {
        double freq_sig;
        /* Frequency of Message Signal */
        printf("Frequency of Message Signal\n");
        scanf("%lf",&freq_sig);
        double freq_carrier;
        /* Frequency of Carrier Signal  */
        printf("Frequency of Carrier Signal\n");
        scanf("%lf",&freq_carrier);
        double messageAmplitude;
        /* Message Amplitude */
        printf("Message Amplitude\n");
        scanf("%lf",&messageAmplitude);
        double carrierAmplitude;
        /* Carrier Amplitude  */
        printf("Carrier Amplitude\n");
        scanf("%lf",&carrierAmplitude);
        double dcOffset;
        /* DC Offset */
        printf("DC Offset\n");
        scanf("%lf",&dcOffset);
        double modIndex;
        /* Modulation Index */
        printf("Modulation Index\n");
        scanf("%lf",&modIndex);
    
        double messageSignal[sampleSize];
        double carrierSignal[sampleSize];
        for(int i = 0; i < sampleSize; i++){
            messageSignal[i] = messageAmplitude * cos(2.0*PI*freq_sig*(period)*i);
            carrierSignal[i] = carrierAmplitude * cos(2.0*PI*freq_carrier*(period)*i);
            result[i] = (dcOffset + messageSignal[i])*carrierSignal[i]; 
        }   
    }
    else {
        double freq_sine;
        printf("\nFrequency of Message Signal:\n");
        scanf("%lf",&freq_sine);
        double sine[sampleSize];
        for(int i = 0; i < sampleSize; i++){
            sine[i] = 1*sin(2.0*PI*freq_sine*(period)*i);
            result[i] = sine[i]; 
        }
    }
}

void chirpGenerator(double* result, double freq_sine) {
    double period = 1/fs;
    double sine[sampleSize];
    for(int i = 0; i < sampleSize; i++){
        sine[i] = 1*sin(2.0*PI*freq_sine*(period)*i);
        result[i] = sine[i]; 
    }   
    
}

#ifdef DEMODSIG
    void signalDemodulator(double* iirOutputSignalTime, double* demodSignal, double* fftDemodSignal) {
        // Getting the absolute value of the signal
        double absSignal[sampleSize],tempSig[sampleSize];

        for (int i=0;i<sampleSize;i++){
            absSignal[i] = (iirOutputSignalTime[i]*iirOutputSignalTime[i]);
        }   
        // Peak detection - Like a Zero-Order-Hold
        for (int i=sampleSize; i--;){	
            //first look for the first maxima, if found, set all the next values to max until next max is found
            if ((absSignal[i-1]>=absSignal[i]) && (absSignal[i-1]>=absSignal[i-2])){
                demodSignal[i-1] = absSignal[i-1];
            } 
            else {
                demodSignal[i-1] = demodSignal[i];
            }
        }
        for (int i=0;i<sampleSize;i++){
            demodSignal[i] = sqrt(demodSignal[i]);
        }

        toTextFile(demodSignal,1,"Pre-LP Demod");

        // Need to lowpass the signal to get smooth output 
        double zfilterNumDemod[zNumLenDemod];
        double zfilterDenDemod[zDenLenDemod];

        /* Single Channel Test */
        iir_coef_lowpass(-1, 
                          zfilterNumDemod, 
                          zfilterDenDemod,
                          a_i_demod, 
                          b_i_demod);
        // ================================================= //
        // IIR Low Pass Filter Implementation                //
        // ================================================= //
        double y_demod[sampleSize];

        for (int i=0;i<sampleSize;i++){
            y_demod[i]=0.000;
        }
        // Temporary storage variables for intermediate values
        double firValDemod;
        double iirValDemod;

        // n represents the current output y(n)
        for (int i=0; i<sampleSize; i++){
            
            // Compute the FIR term:
            firValDemod = 0.0;
            for (int k=0; k<zNumLenDemod; k++){
                // Make sure to avoid index errors 
                if ((i-k) >= 0) {
                    firValDemod += zfilterNumDemod[k]*demodSignal[i-k];
                } 
            }
        
            // Compute the IIR term:
            iirValDemod = 0.0;
            // printf("%lf \n",iirVal);
            for (int k=1; k<zDenLenDemod; k++){
                // Make sure to avoid index errors 
                if ((i-k) >= 0) {
                    iirValDemod += zfilterDenDemod[k]*y_demod[i-k];
                } 
            }     

            // Put the two terms together  
            y_demod[i] = (firValDemod - iirValDemod);
        }

        // Apply the gain to the output:
        for (int i=0; i<sampleSize; i++) {
            y_demod[i] *= gain;
            demodSignal[i]=y_demod[i];
        }
        
        // ================================================= //
        // ================================================= //

        #ifdef PLOTFILTERCOEFFS
            toTextFile(zfilterNum,2,"ZfilterNum1");
            toTextFile(zfilterDen,2,"ZfilterDen1"); 
        #endif

        #ifdef UFFT
            /* FFT on ADC values */
            cmplx demodSignalFFT[sampleSize];                          // Temporary array to be passed into FFT
            dataPrep(demodSignal, demodSignalFFT);
            fft(demodSignalFFT, sampleSize, scratch );
            
            /* Gets the FFT magnitude of the Ouput data */
            cmplx demodSigMagOut[sampleSize];
            fftNormalize(demodSignalFFT,demodSigMagOut,sampleSize);

            /* Converting the FFT Array to a double */
            for (int i=0;i<sampleSize;i++) fftDemodSignal[i] = demodSigMagOut[i].Re;

        #endif
        #ifdef FFTW
            fftw_complex demodSignalTime[sampleSize];
            fftw_complex demodSignalFFT[sampleSize];                          // Temporary array to be passed into FFT
            
            /* DataPrep */
            for (int i=0; i<sampleSize; i++) {
                demodSignalTime[i][REAL] = demodSignal[i];
                demodSignalTime[i][IMAG] = 0;
            }
            
            /* Demod FFT */
            fftw_FFT(demodSignalTime, demodSignalFFT);

            /* Display Data */
            for (int i=0;i<sampleSize;i++) {
                fftDemodSignal[i] = 20*log10( ( sqrt ( pow( demodSignalFFT[i][REAL], 2) + pow(demodSignalFFT[i][IMAG], 2) ) ) );
            }
        #endif
    }
#endif

#ifdef FIR
    // ==================================================================== //
    //                      FIR Implementation Function
    // ==================================================================== //
    /* firImplemenation
    double* adcOutput -> Time Domain output from the ADC 
    double* adcout_FFT_mag_double -> Freq Domain output from the ADC 
    double* y_FFT_mag_double -> Freq Domain output final
    double* y_out -> TFreq Domain output final */
    void firImplemenation(double* adcOutput, double* adcout_FFT_mag_double, double* y_FFT_mag_double, double* y_out) {
        /* Array needed for FFT */
        /* --------------------------  Transversal Filter ------------------------------ */
        #ifdef TRANSVERSAL
            /* --------------------- ADC data ------------------------- */
            /* FFT on ADC values */
            cmplx adcOutFFT[sampleSize];                          // Temporary array to be passed into FFT
            dataPrep(adcOutput, adcOutFFT);
            fft(adcOutFFT, sampleSize, scratch);
            
            /* Gets the FFT magnitude of the adc data */
            cmplx adcMagOut[sampleSize];
            fftNormalize(adcOutFFT,adcMagOut,sampleSize);
            
            /* Converting the FFT Array to a double */
            for (int i=0;i<sampleSize;i++) adcout_FFT_mag_double[i] = adcMagOut[i].Re;
            
            /* ----------------- Time domain output ------------------ */
            transversalFilterFunction(y_out, h_ZeroPad, adcOutput);   

            /* ------------------- FFT output ------------------------ */
            /* FFT on Output values */
            cmplx y_FFT[sampleSize];                          // Temporary array to be passed into FFT
            dataPrep(y_out, y_FFT);
            fft(y_FFT, sampleSize, scratch );
            
            /* Gets the FFT magnitude of the Output data */
            cmplx y_FFT_mag[sampleSize];
            fftNormalize(y_FFT,y_FFT_mag,sampleSize);
            
            /* Converting the FFT Array to a double */
            for (int i=0;i<sampleSize;i++) y_FFT_mag_double[i] = y_FFT_mag[i].Re;
        #endif
        /* ----------------------------------------------------------------------------- */
        #ifdef FASTCONV
            // =============================================================//
            //                              FFT 
            // =============================================================//
            #ifdef UFFT
                    
                /* --------------------- ADC data ------------------------- //
                /* FFT on ADC values */
                cmplx adcOutFFT[sampleSize];                          // Temporary array to be passed into FFT
                dataPrep(adcOutput, adcOutFFT);
                fft(adcOutFFT, sampleSize, scratch );
                
                /* Gets the FFT magnitude of the adc data */
                cmplx adcMagOut[sampleSize];
                fftNormalize(adcOutFFT,adcMagOut,sampleSize);

                /* Converting the FFT Array to a double */
                for (int i=0;i<sampleSize;i++) adcout_FFT_mag_double[i] = adcMagOut[i].Re;
                
                // =============================================================//
                //          Combining the two arrays topether 
                // =============================================================//
                // Multiply the two complex FFTs:
                cmplx output[sampleSize];
                for(int i = 0; i <sampleSize; i++){
                    output[i].Re =  filterFFT[i].Re*adcOutFFT[i].Re - filterFFT[i].Im*adcOutFFT[i].Im;
                    output[i].Im =  filterFFT[i].Re*adcOutFFT[i].Im + filterFFT[i].Im*adcOutFFT[i].Re;
                }
                
                /* Gets the FFT magnitude of the adc data */
                cmplx y_FFT_mag[sampleSize];
                fftNormalize(output, y_FFT_mag, sampleSize); // Displays the whole FFT(can only use half the array for display)
                
                /* Converting the FFT Array to a double */
                for (int i=0; i<sampleSize; i++) y_FFT_mag_double[i] = y_FFT_mag[i].Re;
            
                // =============================================================//
                //                              iFFT 
                // =============================================================//
                // Divide by the sample size to get the correct output
                ifft(output, sampleSize, scratch );
                for (int i=0;i<sampleSize;i++) y_out[i] = output[i].Re /sampleSize;
            #endif

            // =============================================================//
            //                              FFTW 
            // =============================================================//
            #ifdef FFTW
                /* --------------------- ADC data ------------------------- //
                /* FFT on ADC values */
                fftw_complex adcOutputTime[sampleSize];
                fftw_complex adcOutFFT[sampleSize];                          // Temporary array to be passed into FFT
                
                /* DataPrep */
                for (int i=0; i<sampleSize; i++) {
                    adcOutputTime[i][REAL] = adcOutput[i];
                    adcOutputTime[i][IMAG] = 0;
                }
                
                /* ADC FFT */
                fftw_FFT(adcOutputTime, adcOutFFT);

                /* ADC Display Data */
                for (int i=0;i<sampleSize;i++) {
                    adcout_FFT_mag_double[i] = 20*log10( ( sqrt ( pow( adcOutFFT[i][REAL], 2) + pow(adcOutFFT[i][IMAG], 2) ) ) );
                }

                /* Multiplying H and X together */
                fftw_complex output[sampleSize];
                for(int i = 0; i <sampleSize; i++){
                    output[i][REAL] =  filterFFTfftw[i][REAL]*adcOutFFT[i][REAL] - filterFFTfftw[i][IMAG]*adcOutFFT[i][IMAG];
                    output[i][IMAG] =  filterFFTfftw[i][REAL]*adcOutFFT[i][IMAG] + filterFFTfftw[i][IMAG]*adcOutFFT[i][REAL];
                }
                
                /* Gets the FFT magnitude of the adc data */     
                for (int i=0;i<sampleSize;i++) {
                    y_FFT_mag_double[i] = 20*log10( ( sqrt ( pow( output[i][REAL], 2) + pow(output[i][IMAG], 2) ) ) );
                }

                // =============================================================//
                //                              iFFTW 
                // =============================================================//
                /* Time domain ouput */
                fftw_complex y_out_fftw[sampleSize];
                ifftw_FFT(output,y_out_fftw);
                for (int i=0; i<sampleSize; i++) {
                    y_out[i] = y_out_fftw[i][REAL] ;
                }        
            #endif 
        #endif
    }

    // =========================================================== //
    //             Filter coefficients calculator functions        //
    // =========================================================== //
    void coeffCalc(int ch_number) {
        /// Getting the center frequency
        double f_center = ((channel_bw/2)+ starting_channel_freq) + channel_bw*(ch_number-1);
        // Upper frequency
        double fu = f_center + cuttoff_bw/2;
        // Lower frequency
        double fl = f_center - cuttoff_bw/2;
        // Passband ripple and stopband ripple
        double del_p = pow(10,(A_p_min/20)) - 1; 
        double del_s = pow(10,(-A_s_min/20));  
        // Normalised Transition Width(delta_f/fs)
        double delta_fn = delta_f/fs;
        // --- N for Hanning --- //
        #ifdef HANNING
            N = round(3.1/delta_fn);
        #endif
        // --- N for Blackman --- //
        #ifdef HAMMING
            N = round(3.3/delta_fn);
        #endif
        // --- N for Hamming --- //
        #ifdef BLACKMAN
            N = round(5.5/delta_fn);
        #endif
        printf("N = %d",N);
        // Normalized cutoff frequencies        
        double normalised_lower_cutoff = (fl - delta_f/2)/fs;
        double normalised_upper_cutoff = (fu + delta_f/2)/fs;

        filterCoeffCalculator(normalised_lower_cutoff,normalised_upper_cutoff);

        // Zero pad filter coefficients
        zeroPad(h_final, h_ZeroPad, sampleSize);

        #ifdef FASTCONV
            #ifdef UFFT
                /* FFT on filter coefficients */
                dataPrep(h_ZeroPad, filterFFT);
                fft(filterFFT, sampleSize, scratch );
                
                // TODO: DELETE
                #ifdef PLOTFILTERCOEFFS
                    // Gets the FFT magnitude of the filter coefficients 
                    cmplx filterFFTMag[sampleSize];
                    fftNormalize(filterFFT,filterFFTMag,sampleSize);

                    double filterFFTMagdouble[sampleSize];
                    for (int i=0;i<sampleSize;i++){
                        filterFFTMagdouble[i] = filterFFTMag[i].Re;
                    }
                    toTextFile(filterFFTMagdouble,0,"Filter Normal ufft"); 
                #endif
            #endif
            #ifdef FFTW
                /* Data Prep */
                fftw_complex h_ZeroPadfftw[sampleSize];
                for (int i=0; i<sampleSize; i++) {
                    h_ZeroPadfftw[i][REAL] = h_ZeroPad[i];
                    h_ZeroPadfftw[i][IMAG] = 0;
                }

                /* FFT */
                fftw_FFT(h_ZeroPadfftw, filterFFTfftw);
                #ifdef PLOTFILTERCOEFFS
                    // TODO: DELETE
                    // Gets the FFT magnitude of the adc data      
                    double filterFFTMagdouble[sampleSize];
                    for (int i=0;i<sampleSize;i++) {
                        filterFFTMagdouble[i] = 20*log10( ( sqrt ( pow( filterFFTfftw[i][REAL], 2) + pow(filterFFTfftw[i][IMAG], 2) ) ) );
                    }
                    toTextFile(filterFFTMagdouble,0,"Filter Normal fftw"); 
                #endif
            #endif
        #endif
    }

    double bpfFilter(int n_val, double f1, double f2) {
    double h_D_funcval=0;

    if (n_val == 0) {
        h_D_funcval = 2*(f2-f1);
    }
    else {
        h_D_funcval = (((2*f2)*sin(n_val*2*PI*f2))/(n_val*2*PI*f2))-(((2*f1)*sin(n_val*2*PI*f1))/(n_val*2*PI*f1));
    }
    return h_D_funcval;
    }

    void filterCoeffCalculator(double fp1,double fp2){
        double w[N];
        double h_D[N];
        
        int nstop = (N+1)/2;
        
        for (int i=0 ; i < nstop-1 ; i++){
            // Hanning 
            #ifdef HANNING
                w[i] = 0.5 + 0.5*cos((2*PI*i)/N);
            #endif
            // Blackman
            #ifdef BLACKMAN
                w[i]= 0.42 + 0.5*cos(2*PI*(i)/(N-1)) + 0.08*cos(4*PI*(i)/(N-1));
            #endif
            // Hamming
            #ifdef HAMMING
                w[i]= 0.54 + 0.46*cos(2*PI*(i)/N); 
            #endif
            h_D[i] = bpfFilter(i,fp1,fp2);                                          
            h[i] = h_D[i]*w[i];
        }

        // Reversing the array
        double h_beg[nstop];
        int idx = 0;
        for (int i = 0; i < nstop; i++){
            idx = nstop-1 - i;
            h_beg[i] = h[idx];  
        }
        
        int ifx = 0;
        // Comining the two arrays     
        for (int i = 0; i<N; i++) {
            if (i<(N+1)/2){
                h_final[i] = h_beg[i];
            }
            else{
                ifx = (i+1) - round((N+1)/2);
                h_final[i] = h[ifx];
            }
        }
    }

    #ifdef TRANSVERSAL
        /* ============================================================= */
        /*                     Transversal Function                      */
        /* ============================================================= */
        void transversalFilterFunction(double* result, double* h_t, double* x) {
            /* Have to take the zeropadded array of filtercoefficients */    
            for (int i=0; i<sampleSize; ++i) {
                for (int k=(N-1);k--;){    
                    if (i-k>0){
                        result[i] += h_t[k] * x[i-k];
                    }
                }
            }
        }
    #endif

    // =========================================================== //
    //                Data prep for input to FFT                   //
    // =========================================================== //
    //                 Zero Pad the array                          //
    // ------------------------------------------------------------//
    void zeroPad(double *inArr, double* result, int sSize) {
        for (int i=0; i<sSize; i++) {
            if (i < N){
                result[i] = inArr[i];
            }else{
                result[i] = 0;
            }
        }
    }
#endif

#ifdef UFFT 
    // =========================================================== //
    // =========================================================== //
    //                      FFT and IFFT code                      //
    // =========================================================== //
    // =========================================================== //
    /* Copyright (c) 2017 David Barina			                   */
    /* MIT License						                           */
    /* Link: https://www.fit.vut.cz/research/product/510/          */
    /* ___________________________________________________________ */

    // --- FFT --- //
    void fft( cmplx *v, int n_fft, cmplx *tmp )
    {
    if(n_fft>1) {			/* otherwise, do nothing and return */
        int k,m;    cmplx z, w, *vo, *ve;
        ve = tmp; vo = tmp+n_fft/2;
        for(k=0; k<n_fft/2; k++) {
        ve[k] = v[2*k];
        vo[k] = v[2*k+1];
        }
        fft( ve, n_fft/2, v );		/* FFT on even-indexed elements of v[] */
        fft( vo, n_fft/2, v );		/* FFT on odd-indexed elements of v[] */
        for(m=0; m<n_fft/2; m++) {
        w.Re = cos(2*PI*m/(double)n_fft);
        w.Im = -sin(2*PI*m/(double)n_fft);
        z.Re = w.Re*vo[m].Re - w.Im*vo[m].Im;	/* Re(w*vo[m]) */
        z.Im = w.Re*vo[m].Im + w.Im*vo[m].Re;	/* Im(w*vo[m]) */
        v[  m  ].Re = ve[m].Re + z.Re;
        v[  m  ].Im = ve[m].Im + z.Im;
        v[m+n_fft/2].Re = ve[m].Re - z.Re;
        v[m+n_fft/2].Im = ve[m].Im - z.Im;
        }
    }
    return;
    }

    // --- IFFT --- //
    void ifft( cmplx *v, int n_fft, cmplx *tmp )
    {
    if(n_fft>1) {			/* otherwise, do nothing and return */
        int k,m;    cmplx z, w, *vo, *ve;
        ve = tmp; vo = tmp+n_fft/2;
        for(k=0; k<n_fft/2; k++) {
        ve[k] = v[2*k];
        vo[k] = v[2*k+1];
        }
        ifft( ve, n_fft/2, v );		/* FFT on even-indexed elements of v[] */
        ifft( vo, n_fft/2, v );		/* FFT on odd-indexed elements of v[] */
        for(m=0; m<n_fft/2; m++) {
        w.Re = cos(2*PI*m/(double)n_fft);
        w.Im = sin(2*PI*m/(double)n_fft);
        z.Re = w.Re*vo[m].Re - w.Im*vo[m].Im;	/* Re(w*vo[m]) */
        z.Im = w.Re*vo[m].Im + w.Im*vo[m].Re;	/* Im(w*vo[m]) */
        v[  m  ].Re = (ve[m].Re + z.Re);
        v[  m  ].Im = (ve[m].Im + z.Im);
        v[m+n_fft/2].Re = (ve[m].Re - z.Re);
        v[m+n_fft/2].Im = (ve[m].Im - z.Im);
        }
    }
    return;
    }

    /* ============================================================= */
    /*                     Normalizing the FFT                       */
    /* ______________________________________________________________*/
    /* Obtaining the normalized logorithmic values from FFTAlgo      */
    /* ============================================================= */
    void fftNormalize(cmplx *fftArr, cmplx* result, int numN) {
        // magnitude 
        for (int i = 0; i<numN; i++ ) {
            result[i].Re = 20*log10( ( sqrt ( pow( fftArr[i].Re, 2) + pow(fftArr[i].Im, 2) ) ) );
            result[i].Im = 0;
        }
    }

    void dataPrep(double *inArr, cmplx* result) {
        for (int i=0;i<sampleSize;i++) { 
            result[i].Re = inArr[i];
            result[i].Im = 0;
        }
    }
#endif

#ifdef FFTW
    /* Computes the 1-D fast fourier transform */
    void fftw_FFT(fftw_complex *in, fftw_complex *out) {
        // create a DFT plan
        fftw_plan plan = fftw_plan_dft_1d(sampleSize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        // execure the plan
        fftw_execute(plan);
        // do some cleaning
        fftw_destroy_plan(plan);
        fftw_cleanup();
    }

    /* Computes the 1-D inverse fast fourier transform */
    void ifftw_FFT(fftw_complex *in, fftw_complex *out) {
        // create a IDFT plan
        fftw_plan plan = fftw_plan_dft_1d(sampleSize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        // execure the plan
        fftw_execute(plan);
        // do some cleaning
        fftw_destroy_plan(plan);
        fftw_cleanup();
        // scale the output to obtain the exact inverse
        for (int i=0; i<sampleSize; i++) {
            out[i][REAL] /=sampleSize;
            out[i][IMAG] /=sampleSize;
        }
    }
#endif

#ifdef IIR_PASCALS 
    /* ====================================================================================== */
    /* ====================================================================================== */
    /*                                                                                        */
    /*                                     IIR                                                */
    /*                                                                                        */
    /* ====================================================================================== */
    /* ====================================================================================== */

    /* ================================================================ */
    /*                         FILTER THE SIGNAL                        */
    /* ================================================================ */
    void iirImplementation(double* input_sig, double* z_filter_num, double* z_filter_den, double* y_n) {

        double firVal;
        double iirVal;
        
        int zNumLen = COL;
        int zDemLen = COL;


        for (int i=0;i<sampleSize;i++){
            y_n[i]=0.000;
        }

        // n represents the current output y(n)
        for (int i=0; i<sampleSize; i++){
            
            // Compute the FIR term:
            firVal = 0.0;
            for (int k=0; k<zNumLen; k++){
                // Make sure to avoid index errors 
                if ((i-k) >= 0) {
                    firVal += z_filter_num[k]*input_sig[i-k];
                } 
            }
        
            // Compute the IIR term:
            iirVal = 0.0;
            // printf("%lf \n",iirVal);
            for (int k=1; k<zDemLen; k++){
                // Make sure to avoid index errors 
                if ((i-k) >= 0) {
                    iirVal += z_filter_den[k]*y_n[i-k];
                } 
            }     

            // Put the two terms together  
            y_n[i] = (firVal - iirVal);
        }

        // Apply the gain to the output:
        for (int i=0; i<sampleSize; i++) {
            y_n[i] *= gain;
        }
    }

    /* ================================================================ */
    /*                FILTER COEFFICIENTS BANDPASS                      */
    /* ================================================================ */
    void iir_coef_bandpass(int ch_number, double* z_filter_num, double* z_filter_den, double* a_i, double* b_i) {

        double f_center = ((channel_bw/2)+ starting_channel_freq) + channel_bw*(ch_number-1);
        double fu = f_center + cuttoff_bw/2;
        double fl = f_center - cuttoff_bw/2;

        // please for the love of all that is good in this word, initialize the array
        double C[COL][COL];
        for(int i=0;i<(COL);i++) {
            for(int j=0;j<(COL);j++) {
                C[i][j]=0;
            }
        }
        /* Setting Up Pascals Matrix Outer Bounds */
        for(int i=0;i<(COL);i++) {
            for(int j=0;j<(COL);j++) {
                if (j == 0) { 
                    C[i][j] = pow((-1),i)*factorial(n)/(factorial(i)*factorial(n-i));
                }
                if (j == n) {
                    C[i][j] = factorial(n)/(factorial(i)*factorial(n-i));
                }
                if (i == 0) {
                    C[i][j] = 1;
                }
                if (i == n) {
                    C[i][j] = pow((-1),j);
                }
            }
        }
        
        /* Setting up the pascals matrix */
        //  Rows
        for (int j=1;j<(COL);j++) {
            //  Cols 
            for (int i=1;i<(COL);i++) {
                // Setup of Outer Bounds of Matrix
                C[i][j] = C[i-1][j]+C[i-1][j-1]+C[i][j-1]; 
            }
        }

        /* print */
        // Rows
        for (int j=0;j<(COL);j++) {
            //  Cols 
            for (int i=0;i<(COL);i++) {
            }
        }

        // Calculate the neccessary constants that are involved in prewarping the cut-off frequencies.
        double tl = tan((PI* fl)/fs);
        double cu = 1/(tan((PI* fu)/fs));

        double u = cu/(1-cu*tl);
        double l = tl/(1-cu*tl);

        /* ------------------------------------------------------------------------------ */
        /*      Calculate Matrix [Del_BP]. Set up the diagonals (U to U^n and L to L^n)   */
        /* ------------------------------------------------------------------------------ */
        // Set up where the middle coloumn of the matrix will be:
        int middle = (int)((COL)/2);

        /* ------------------------------------------------------------------------------ */
        /*               Calculating Del_Ai and Del_Bi Matrices:                          */
        /* ------------------------------------------------------------------------------ */
                    //Rows //Cols
        double del_i[ ROW ][ COL ];

        // Rows
        for(int i=0; i < ROW; i++) {
            // Cols
            for(int j=0; j < COL; j++) {
                del_i[i][j] = 0;
            }
        }
        
        /* Setting up of the triangle*/
        // Jumps over Rows
        for (int k=0; k<ROW; k++) {
            if (k==0) {
                del_i[k][middle] = 1;
            }
            if (k>0) {
                del_i[k][middle-k] = pow(u,k);
                del_i[k][middle+k] = pow(l,k);
            }
        }

        // rows
        for (int k=2; k<ROW; k++){
            // cols
            int start = (int)(n/2) - k + 2;
            int i=1;
            while (i<k) {
                del_i[k][start] = pow(u, (k-i)) * pow(l, i) * (factorial(k)/(factorial(i)*factorial(k-i)));
                start += 2;
                i += 1;
            }
        }

        // create transpose
        double trans[COL][ROW];
        // Rows
        for(int i = 0; i < ROW; i++){
            for(int j = 0; j < COL; j++){
                trans[j][i] = del_i[i][j];
            }
        }
        
        // Calculate the del_Ai and del_Bi Matrix by dot producting the Del_BP matrix with the given ai and bi coefficients of the analogue filter
        double delAi[COL];
        double delBi[COL];

        for (int i=0; i<COL; i++) {
            delAi[i] = dot(a_i,trans[i],ROW); 
            delBi[i] = dot(b_i,trans[i],ROW);
        }

        // Looking at each row in P matrix and multiplying it by the coloumn matrix (del_A or del_B)
        double sum = 0.00;
        for(int i =0; i<COL; i++){
            sum = 0.00;
            // Looking at each coloumn in row i of P matrix and multiplying by each row in
            // coloumn matrix
            for(int j=0; j<COL; j++){
                sum += C[i][j]*delAi[j];
            }
            z_filter_num[i] = sum;
        }

        // Looking at each row in P matrix and multiplying it by the coloumn matrix (del_A or del_B)
        sum = 0.00;
        for(int i =0; i<COL; i++){
            sum = 0.00;
            // Looking at each coloumn in row i of P matrix and multiplying by each row in
            // coloumn matrix
            for(int j=0; j<COL; j++){
                sum += C[i][j]*delBi[j];
            }
            z_filter_den[i] = sum;
        }
        
        // We dont want the gain to be fixed
        #ifndef PLOTFILTERCOEFFS
            // Fixing the gain
            gain = 1/z_filter_den[0];
            for (int i=0; i<COL; i++) z_filter_den[i] *= gain;
        #endif
    }

    /* ================================================================ */
    /*                FILTER COEFFICIENTS LOWPASS                       */
    /* ================================================================ */
    void iir_coef_highpass(int ch_number, double* z_filter_num, double* z_filter_den, double* a_i, double* b_i) {

        double f_center = ((channel_bw/2)+ starting_channel_freq) + channel_bw*(ch_number-1);
        double fc = f_center - cuttoff_bw/2;
        printf("Center Freqeuncy is %lf\n",f_center);
        
        // please for the love of all that is good in this word, initialize the array
        double C[COL][COL];
        for(int i=0;i<(COL);i++) {
            for(int j=0;j<(COL);j++) {
                C[i][j]=0;
            }
        }
        /* Setting Up Pascals Matrix Outer Bounds */
        for(int i=0;i<(COL);i++) {
            for(int j=0;j<(COL);j++) {
                if (j == 0) { 
                    C[i][j] = pow((-1),i)*factorial(n)/(factorial(i)*factorial(n-i));
                }
                if (j == n) {
                    C[i][j] = factorial(n)/(factorial(i)*factorial(n-i));
                }
                if (i == 0) {
                    C[i][j] = 1;
                }
                if (i == n) {
                    C[i][j] = pow((-1),j);
                }
            }
        }
        
        /* Setting up the pascals matrix */
        //  Rows
        for (int j=1;j<(COL-1);j++) {
            //  Cols 
            for (int i=1;i<(COL-1);i++) {
                // Setup of Outer Bounds of Matrix
                C[i][j] = C[i-1][j]+C[i-1][j-1]+C[i][j-1];
            }
        }

        // Calculate the neccessary constants that are involved in prewarping the cut-off frequencies.
        double t = tan(PI * fc/fs);
        double c = 1/(tan(PI* fc/fs));

        // Calculate the del_Ai and del_Bi Matrix by dot producting the Del_BP matrix with the given ai and bi coefficients of the analogue filter
        double delAi[COL];
        double delBi[COL];

        for (int i=0; i<COL; i++) {
            delAi[i] = a_i[i]*pow(t,i); 
            delBi[i] = b_i[i]*pow(t,i);
        }

        // Looking at each row in P matrix and multiplying it by the coloumn matrix (del_A or del_B)
        double sum = 0.00;
        for(int i =0; i<COL; i++){
            sum = 0.00;
            // Looking at each coloumn in row i of P matrix and multiplying by each row in
            // coloumn matrix
            for(int j=0; j<COL; j++){
                sum += C[i][j]*delAi[j];
            }
            z_filter_num[i] = sum;
        }

        // Looking at each row in P matrix and multiplying it by the coloumn matrix (del_A or del_B)
        sum = 0.00;
        for(int i =0; i<COL; i++){
            sum = 0.00;
            // Looking at each coloumn in row i of P matrix and multiplying by each row in
            // coloumn matrix
            for(int j=0; j<COL; j++){
                sum += C[i][j]*delBi[j];
            }
            z_filter_den[i] = sum;
        }

        // We dont want the gain to be fixed
        #ifndef PLOTFILTERCOEFFS
            // Fixing the gain
            gain = 1/z_filter_den[0];
            for (int i=0; i<COL; i++) z_filter_den[i] *= gain;
        #endif
    }

#endif

    /* ================================================================ */
    /*                FILTER COEFFICIENTS LOWPASS                       */
    /* ================================================================ */
#if defined IIR_PASCALS || defined DEMODSIG
    void iir_coef_lowpass(int ch_number, double* z_filter_num, double* z_filter_den, double* a_i, double* b_i) {
        int columnVar = 0;
        int new_n = 0;
        float fc = 0;
        if (ch_number != 1) {
            #ifdef IIR_PASCALS
                double f_center = ((channel_bw/2)+ starting_channel_freq) + channel_bw*(ch_number-1);
                fc = f_center + cuttoff_bw/2;
                columnVar = COL;
                new_n = n;
                printf("Column Var is %d\n",columnVar);    
            #endif
        }
        if (ch_number == -1) {
            #ifdef DEMODSIG
                // This is for low passing the output
                fc = f_c_sig_demod;
                columnVar = demodOrder + 1;
                new_n = demodOrder;
            #endif
        }

        // please for the love of all that is good in this word, initialize the array
        double C[columnVar][columnVar];
        for(int i=0;i<(columnVar);i++) {
            for(int j=0;j<(columnVar);j++) {
                C[i][j]=0;
            }
        }
        /* Setting Up Pascals Matrix Outer Bounds */
        for(int i=0;i<(columnVar);i++) {
            for(int j=0;j<(columnVar);j++) {
                if (j == 0) { 
                    C[i][j] = factorial(new_n)/(factorial(i)*factorial(new_n-i));
                }
                if (j == new_n) {
                    C[i][j] = pow((-1),i)*factorial(new_n)/(factorial(i)*factorial(new_n-i));
                }
                if (i == 0) {
                    C[i][j] = 1;
                }
                if (i == new_n) {
                    C[i][j] = pow((-1),j);
                }
            }
        }
        
        /* Setting up the pascals matrix */
        //  Rows
        for (int i=1;i<(columnVar);i++) {
            for (int j=1;j<(columnVar);j++) {
            //  columnVars 
                // Setup of Outer Bounds of Matrix
                C[i][j] = +C[i][j-1]-C[i-1][j]-C[i-1][j-1]; // TODO Check this
            }
        }

        // Calculate the neccessary constants that are involved in prewarping the cut-off frequencies.
        double t = tan((PI* fc)/fs);
        double c = 1/(tan((PI* fc)/fs));

        // Calculate the del_Ai and del_Bi Matrix by dot producting the Del_BP matrix with the given ai and bi coefficients of the analogue filter
        double delAi[columnVar];
        double delBi[columnVar];

        for (int i=0; i<columnVar; i++) {
            delAi[i] = a_i[i]*pow(c,i); 
            delBi[i] = b_i[i]*pow(c,i);
        }
        
        // Looking at each row in P matrix and multiplying it by the columnVaroumn matrix (del_A or del_B)
        double sum = 0.00;
        for(int i =0; i<columnVar; i++){
            sum = 0.00;
            // Looking at each columnVaroumn in row i of P matrix and multiplying by each row in
            // columnVaroumn matrix
            for(int j=0; j<columnVar; j++){
                sum += C[i][j]*delAi[j];
            }
            z_filter_num[i] = sum;
        }

        // Looking at each row in P matrix and multiplying it by the columnVaroumn matrix (del_A or del_B)
        sum = 0.00;
        for(int i =0; i<columnVar; i++){
            sum = 0.00;
            // Looking at each columnVaroumn in row i of P matrix and multiplying by each row in
            // columnVaroumn matrix
            for(int j=0; j<columnVar; j++){
                sum += C[i][j]*delBi[j];
            }
            z_filter_den[i] = sum;
        }

        // We dont want the gain to be fixed
        #ifndef PLOTFILTERCOEFFS
            // Fixing the gain
            gain = 1/z_filter_den[0];
            for (int i=0; i<columnVar; i++) z_filter_den[i] *= gain;
        #endif
    }

    // =============================================================//
    //                      DOT PRODUCT
    // =============================================================//
    double dot(double v[], double u[], int dotNum) {
        
        double result = 0.0;
        
        for (int i = 0; i < dotNum; i++)
            result += v[i]*u[i];

        return result;
    }

    /* ================================================================ */
    /*                    Factorial using recursions                    */
    /* ================================================================ */
    double factorial(int p) {  
    if (p == 0)  
        return 1;  
    else  
        return(p * factorial(p-1));  
    }  

    /* ================================================================ */
    /*                    Max value between two numbers                 */
    /* ================================================================ */
    double max(double num1, double num2) {
        return (num1 > num2 ) ? num1 : num2;
    }

#endif