/**
 * Mario Romero March 2021
**/

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <time.h>
#include "define.h"
using namespace std;

enum maps { Alpha, Jota, MeanIntensity, Flux, RadPressure, N_maps}; //alpha and j, I used the Spanish name of 'j' to avoid issues with the iterator variable 'j'
#include "tables.cpp"
#include "DormandPrince.cpp"
const DormandPrince DP; //This object only contains constant numbers.
#include "ConversorString.cpp"
#include "read_input.cpp"
#include "write_output.cpp"
#include "help.cpp"

//Functions
bool approx(r_number,r_number);
r_number ToRadians(r_number);
r_number s_to_z(r_number,r_number,r_number);
r_number s_to_R(r_number,r_number,r_number,r_number);

int main(int argc, char** argv)
{
    
    time_t real_tstart = time(NULL);
    
    if(argc<7){how_to_run_this_code();}
    //READ INPUT
    const int N_R = (int) (todouble(argv[1]));
    const int N_z = (int) (todouble(argv[2]));
    const int N_nu = (int) (todouble(argv[3]));
    
    
    const r_number* I0 = readExternalBackground(argv[4],N_nu);
    //const r_number I0 = 2.0; //PLACEHOLDER
    
    Table** magnitude = readFiles(argc,argv);
    
    /*MAIN LOOP*/
    //The structure of the main loop is a glorified Dormand-Prince on the line of sight (which will be the penultimate nested loop when I add frequencies)
    
    //Get the sphere whose radius is equal to the diagonal of the domain.
    #define MAP_ELEM Alpha //You can use any 'map' element, but to avoid writing 'Alpha' too much...
    r_number Rlast = magnitude[MAP_ELEM][0].X(N_R-1);
    r_number zlast = magnitude[MAP_ELEM][0].Y(N_z-1);
    r_number D2 = Rlast*Rlast + zlast*zlast; //D^2. Radius (Squared)
    
    //I'm declaring here the variables needed inside the loop.
    //r_number R_root[2],z_root[2]; //Points in which the line of sight crosses the sphere
    r_number s_root[2]; //Distance from the point (0,Ri,zj) to the sphere of radius D
    r_number theta,phi;//angles of the line in 3d
    r_number cotg,SinTheta,CosTheta; //cotangent (depreciated), Sin(theta) and Cos(theta). I wrote them with capital C and S to avoid confusion with the mathematical functions in cmath
    r_number SinPhi;//,CosPhi; //Sin(phi) and Cos(phi). I used variables instead of the mathematical functions to avoid being calculated again and again (which is slower)
    r_number a,b,c,sqrt_term; //Elements of a 2nd grade equation: a*x^2+b*x+c=0, and the squared term (sqrt(b^2-4ac)).
    
    r_number R_p,z_p,s_p,I_p; //Points at an arbitrary s. Same with I_p (needed for middle calculations)
    r_number s; //Line of sight coordinate
    r_number h; //stepsize
    r_number currI,nextI; //I_n and I_n+1, respectively
    r_number K[6]; //the 'k' of the runge-kutta methods, for this one, I need six of them
    //#if FIXED_STEP == 0
    //    r_number errorI; //Estimation of the error
    //    r_number control_param; //control parameter for the adaptive stepsize 
    //#endif
    //r_number FourPiJ[N_nu]; //4piJ. This is the result we are looking for.
    r_number dJ, dF, dP; //differential of J, F and P for each line of sight.
    r_number Ri,zj; //Points of the grid, R[i] and z[j]
    
    bool stopDP; //Condition to stop Dormand-Prince integrator.
    r_number n_steps = 2.0*max(N_R,N_z);
    
    //ofstream test_file;
    //test_file.open("Sphere.txt");
    
    //for(int i=30;i<N_R;i+=N_R)
    for(int i=0;i<N_R;++i)
    {
        Ri = magnitude[MAP_ELEM][0].X(i);
        
        //for(int j=30;j<N_z;j+=N_z)
        for(int j=0;j<N_z;++j)
        {
            zj = magnitude[MAP_ELEM][0].Y(j);
            
            //cout<<"iteration (Ri = "<<Ri<<", zj = "<<zj<<") "<<endl;
            
            //for(r_number theta_deg=90.0;theta_deg<180.0;theta_deg+=400.0) //That condition is done in purpose to only have one iteration for testing purposes
            for(r_number theta_deg=THETA_0; theta_deg<180.0; theta_deg+=THETA_STEP)
            {
                theta = ToRadians(theta_deg);
                CosTheta = cos(theta);
                SinTheta = sin(theta);
                
                //for(r_number phi_deg=90.0;phi_deg<180.0;phi_deg+=400.0) //Again, this condition is done in purpose to only have one iteration due to testing
                for(r_number phi_deg=PHI_0;phi_deg<180.0;phi_deg+=PHI_STEP)
                {
                    phi = ToRadians(phi_deg);
                    SinPhi = sin(phi);
                    //Find the two points that the line with angles theta and phi that crosses (0,Ri,zj) hit the sphere of radius D (see my notes)
                    a = 1.0;
                    b = 2.0*(Ri*SinTheta*SinPhi + zj*CosTheta);
                    c = Ri*Ri + zj*zj - D2;
                    sqrt_term = sqrt( b*b - 4.0*a*c ) ;
                    //Roots
                    s_root[0] = (-b + sqrt_term) / (2.0*a);
                    s_root[1] = (-b - sqrt_term) / (2.0*a);
                    
                    //cout<<Ri<<" "<<zj<<" "<<theta_deg<<" "<<phi_deg<<endl;
                    
                    for(int nu=0;nu<N_nu;++nu)
                    {
                        //for(int root=0;root<2;root+=5) //Only first root for testing
                        for(int root=0;root<2;++root) 
                        {
                            
                            nextI = I0[nu]; //For the situations that the if block is not run.
                            s = s_root[root];
                            //cout<<Ri<<" "<<zj<<" "<<s<<" "<<SinPhi<<" "<<SinTheta<<" "<<CosTheta<<endl;
                            if(!approx(s,0.0)) //Point is not at the border, or at the border but the ray is pointing inwards the domain
                            {
                                //sign = (root == 0) ? 1 : -1;
                                stopDP = false;
                                h = s/n_steps; //If adaptive step is on, this h is only used for the first iteration.
                                currI = I0[nu];
                                //cout<<currI<<" "<<CosTheta<<endl;
                                //DormandPrince iterator
                                //int n_iterations = 0; //Test
                                while(!stopDP)
                                {
                                    //n_iterations++;
                                    if(abs(s) < abs(h)) //Last iteration
                                    {
                                        h = s;
                                        stopDP = true;
                                    }
                                    //cout<<s<<" "<<currI<<" "<<stopDP<<endl;
                                    /*Euler method*/
                                    #if METHOD == 0
                                        R_p = s_to_R(Ri,s,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s,CosTheta);
                                    
                                    
                                        nextI = currI + abs(h)*( magnitude[Jota][nu].get_value( R_p , z_p ) - magnitude[Alpha][nu].get_value(R_p , z_p)*currI );
                                        currI = nextI;
                                    /*Midpoint method*/
                                    #elif METHOD == 1                                        
                                        R_p = s_to_R(Ri,s,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s,CosTheta);
                                        //First evaluation
                                        K[0] = abs(h) * ( magnitude[Jota].get_value( R_p , z_p ) - magnitude[Alpha].get_value(R_p , z_p)*currI );
                                        
                                        //Second evaluation
                                        R_p = s_to_R(Ri,s-0.5*h,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s-0.5*h,CosTheta);
                                        K[1] = abs(h) * ( magnitude[Jota].get_value( R_p , z_p ) - magnitude[Alpha].get_value(R_p , z_p)*(currI + 0.5*K[0]) );
                                        //Evaluation
                                        nextI = currI + K[1];
                                        currI = nextI;
                                    /*Dormand-Prince method*/
                                    #else
                                        //First evaluation
                                        s_p = s;
                                        R_p = s_to_R(Ri,s_p,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s_p,CosTheta);
                                        I_p = currI;
                                        K[0] = abs(h) * ( magnitude[Jota][nu].get_value( R_p , z_p ) - magnitude[Alpha][nu].get_value(R_p , z_p)*I_p );
                                        //Second evaluation
                                        s_p = s - DP.c2*h;
                                        R_p = s_to_R(Ri,s_p,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s_p,CosTheta);
                                        I_p = currI + DP.a21*K[0];
                                        K[1] = h * ( magnitude[Jota][nu].get_value( R_p , z_p ) - magnitude[Alpha][nu].get_value(R_p , z_p)*I_p );
                                        //Third evaluation
                                        s_p = s - DP.c3*h;
                                        R_p = s_to_R(Ri,s_p,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s_p,CosTheta);
                                        I_p = currI + DP.a31*K[0] + DP.a32*K[1];
                                        K[2] = abs(h) * ( magnitude[Jota][nu].get_value( R_p , z_p ) - magnitude[Alpha][nu].get_value(R_p , z_p)*I_p );
                                        //Fourth evaluation
                                        s_p = s - DP.c4*h;
                                        R_p = s_to_R(Ri,s_p,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s_p,CosTheta);
                                        I_p = currI + DP.a41*K[0] + DP.a42*K[1] + DP.a43*K[2];
                                        K[3] = abs(h) * ( magnitude[Jota][nu].get_value( R_p , z_p ) - magnitude[Alpha][nu].get_value(R_p , z_p)*I_p );
                                        //Fifth evaluation
                                        s_p = s - DP.c5*h;
                                        R_p = s_to_R(Ri,s_p,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s_p,CosTheta);
                                        I_p = currI + DP.a51*K[0] + DP.a52*K[1] + DP.a53*K[2] + DP.a54*K[3];
                                        K[4] = abs(h) * ( magnitude[Jota][nu].get_value( R_p , z_p ) - magnitude[Alpha][nu].get_value(R_p , z_p)*I_p );
                                        //Sixth evaluation
                                        s_p = s - DP.c6*h;
                                        R_p = s_to_R(Ri,s_p,SinTheta,SinPhi);
                                        z_p = s_to_z(zj,s_p,CosTheta);
                                        I_p = currI + DP.a61*K[0] + DP.a62*K[1] + DP.a63*K[2] + DP.a64*K[3] + DP.a65*K[4];
                                        K[5] = abs(h) * ( magnitude[Jota][nu].get_value( R_p , z_p ) - magnitude[Alpha][nu].get_value(R_p , z_p)*I_p );
                                        
                                        //Estimation
                                        nextI  = currI + DP.b1*K[0] + DP.b2*K[1] + DP.b3*K[2] + DP.b4*K[3] + DP.b5*K[4] + DP.b6*K[5];
                                    #endif
                                    
                                    //#if FIXED_STEP == 1
                                    
                                        currI = nextI;
                                        s -= h;
                                    
                                    /*
                                    #else
                                        errorI = sqrt( (DP.b1-DP.B1)*(DP.b1-DP.B1)*K[0]*K[0] + (DP.b2-DP.B2)*(DP.b2-DP.B2)*K[1]*K[1] + 
                                                (DP.b3-DP.B3)*(DP.b3-DP.B3)*K[2]*K[2] + (DP.b4-DP.B4)*(DP.b4-DP.B4)*K[3]*K[3] +
                                                (DP.b5-DP.B5)*(DP.b5-DP.B5)*K[4]*K[4] + (DP.b6-DP.B6)*(DP.b6-DP.B6)*K[5]*K[5]);
                                    
                                        //Adaptive stepsize
                                        if(!stopDP) //If I'm in the last iteration, stepsize doesn't need to be adapted!
                                        {
                                            control_param = DP_TOLERANCE*max(currI,nextI) / errorI;

                                            //cout<<s<<" "<<h<<" "<<control_param<<endl;
                                            //float dummy;
                                            //cin>>dummy;
                                            
                                            if(approx(control_param,1.0))
                                            //if(control_param >= 1.0)
                                            {//Neither too slow or too inaccurate. You're good to go
                                                currI = nextI;
                                                s -= h;  //Light has travel a distance h, so 's' is reduced
                                            }
                                            else if(control_param > 1.0)
                                            {//Too slow, but I'm accepting the result nonetheless
                                                s -= h;
                                                if(errorI != 0) //This is to avoid a division by 0
                                                {
                                                    h *= pow( abs(DP_TOLERANCE*max(currI,nextI)/errorI) , 0.2 ); //DP method is fifth-order, hence the 1/5 exponent
                                                }
                                                currI = nextI;
                                            }
                                            else
                                            {//inaccurate, repeat the step
                                                h *= pow( abs(DP_TOLERANCE*max(currI,nextI)/errorI) , 0.2 );
                                            }
                                        }
                                    #endif
                                    */
                                    
                                    
                                }//Closes the while iterator
                            //cout<<Ri<<" "<<zj<<" "<<n_iterations<<endl;;
                            }//Ends if not at border pointing outwards
                            //else{cout<<"Skip"<<endl;}
                            
                            //Add to the mean intensity integral
                            dJ = nextI * SinTheta * ToRadians(THETA_STEP) * ToRadians(PHI_STEP) / (4.0*PI);
                            //dJ = nextI;
                            //Add to the flux integral
                            //dF = dJ * CosTheta * 4.0*PI; //Cosine changes sign for the other root, careful!
                            //dF = (root == 0) ? dF : -dF;
                            
                            //Add to the radiation pressure integral
                            //dP = dJ * CosTheta * CosTheta * 4.0*PI * LIGHT_SPEED ; //Both roots have their cosine with different signs, but due to I'm squaring that value...
                            
                            //cout<<nextI<<endl;
                            //magnitude[MeanIntensity][nu].FillZ(nextI,i,j);
                            magnitude[MeanIntensity][nu].SumZ(dJ,i,j); //SumZ is the same as z[i][j] += dJ.
                            //magnitude[Flux][nu].SumZ(dF,i,j);
                            //magnitude[RadPressure][nu].SumZ(dP,i,j); //BUG: Results are displaced by some factor...
                            
                        }//Closes roots
                    }//Closes nu
                    
                }//Closes phi  
                //cout<<Ri<<" "<<zj<<" "<<nextI<<endl;
                //magnitude[MeanIntensity].FillZ(nextI,i,j);
            }//Closes theta
            //break;
        
        }//Closes z
    //break;
    }//Closes R
    
    //test_file.close();
    
    /*PRINT OUTPUT*/
    writeFile(magnitude[MeanIntensity],N_nu,"MeanIntensity.txt","# R | z | J => ");
    
    //Not properly tested
    //writeFile(magnitude[Flux],N_nu,"Flux.txt","# R | z | F => ");
    //writeFile(magnitude[RadPressure],N_nu,"RadPressure.txt","# R | z | P => ");
    
    //Testing...
    //magnitude[Alpha][0].write_table("alphaTest.txt"," R | z | alpha ");
    //magnitude[Jota][0].write_table("jotaTest.txt", " R | z | jota ");
    //magnitude[MeanIntensity].write_table("J.txt", "R | z | J ");
    
    /*DELETE THINGS */
    delete[] I0;
    deleteTables(magnitude);
    
    time_t real_tend = time(NULL);
    cout<<"time elapsed: "<<difftime(real_tend,real_tstart)<<" s"<<endl;
    
    return 0;
}

bool approx(r_number A, r_number B)
{
    //Checks if A is approximately B, using APPROX_TOLERANCE
    if( abs(A-B) <= APPROX_TOLERANCE*abs(A) ){return true;}
    else{ return false; }
}

r_number ToRadians(r_number angle)
{
    return angle*PI/180.0;
}

r_number s_to_z(r_number zj, r_number sp, r_number CosTheta)
{
    return zj + sp*CosTheta;
}

r_number s_to_R(r_number Ri, r_number sp, r_number SinTheta, r_number SinPhi)
{
    r_number radicant = Ri*Ri + 2.0*Ri*sp*SinTheta*SinPhi + sp*sp*SinTheta*SinTheta;
    //cout<<Ri<<" "<<sp<<" "<<radicant<<endl;
    if(radicant >= 0.0){ return sqrt(radicant); }
    else if( -radicant < SOFT_ZERO ){ return 0.0; }//Sounds stupid, but I found a case with Ri=sp and both Sin=1 where radicant gave -1e-13. The exact result is 0, and I don't want the code to explode because of that.
    else
    {
        cout<<"Arguments(Ri,sp,SinTheta,SinPhi):"<<endl;
        cout<<Ri<<" "<<sp<<" "<<SinTheta<<" "<<SinPhi<<endl;
        cout<<"Radicant = "<<radicant<<endl;
        throw runtime_error("Negative square root in s_to_R function.");
        //Better to stop this here rather than let it propagate and finish with a segment fault because of this.
    }
    //return sqrt( Ri*Ri + 2.0*Ri*sp*SinTheta*SinPhi + sp*sp*SinTheta*SinTheta ); 
}
