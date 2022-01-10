#ifndef DORMAND_PRINCE_CPP_INCLUDED
#define DORMAND_PRINCE_CPP_INCLUDED

//I'm only storing here the weights of the method. Taking the numbers from the book 'Numerical Recipes'

struct DormandPrince
{
    r_number c2 = 1.0/5.0;
    r_number c3 = 3.0/10.0;
    r_number c4 = 4.0/5.0;
    r_number c5 = 8.0/9.0;
    r_number c6 = 1.0;
    
    r_number b1 = 35.0/384.0;
    r_number b2 = 0.0;
    r_number b3 = 500.0/1113.0;
    r_number b4 = 125.0/192.0;
    r_number b5 = -2187.0/6784.0;
    r_number b6 = 11.0/84.0;
    
    r_number B1 = 5179.0/57600.0;
    r_number B2 = 0.0;
    r_number B3 = 7571.0/16695.0;
    r_number B4 = 393.0/640.0;
    r_number B5 = -92097.0/339200.0;
    r_number B6 = 187.0/2100.0;
    
    r_number a21 = 1.0/5.0;
    
    r_number a31 = 3.0/40.0;
    r_number a32 = 9.0/40.0;
    
    r_number a41 = 44.0/45.0;
    r_number a42 = -56.0/15.0;
    r_number a43 = 32.0/9.0;
    
    r_number a51 = 19372.0/6561.0;
    r_number a52 = -25360.0/2187.0;
    r_number a53 = 64448.0/6561.0;
    r_number a54 = -212.0/729.0;
    
    r_number a61 = 9017.0/3168.0;
    r_number a62 = -355.0/33.0;
    r_number a63 = 46732.0/5247.0;
    r_number a64 = 49.0/176.0;
    r_number a65 = -5103.0/18656.0;
};



#endif
