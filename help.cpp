#ifndef HELP_CPP_INCLUDED
#define HELP_CPP_INCLUDED

/*
 This function fires if you run the .exe without enough arguments (at the time of writing, less than 6, counting with the .exe).
 It will explain how to run this function and some shortcomings when running.
 
 Mario Romero       March 2021
 */

void how_to_run_this_code()
{
    cout<<endl<<"Syntax: ./code.exe N_R N_z N_nu IC.txt alpha-file.txt jota-file.txt"<<endl<<endl;
    
    cout<<"Last two arguments are the alpha (a) and jota(j) for the equation of radiative transfer:"<<endl;
    cout<<"\t dI/ds = j - a*I"<<endl<<endl;
    
    cout<<"The structure of both files are:"<<endl;
    cout<<" R(pc) | z(pc) | F(R,z) (units) => F_1 | F_2 | F_3 (...) F_N_nu "<<endl;
    cout<<" R0 z0   F_1,00   F_2,00   F_3,00         F_N_nu,00  "<<endl;
    cout<<" R0 z1   F_1,01   F_2,01   F_3,01   ...   F_N_nu,01  "<<endl;
    cout<<" R0 z2   F_1,02   F_2,02   F_3,02         F_N_nu,02  "<<endl;
    cout<<"                    .                     .          "<<endl;
    cout<<"                    .                     .          "<<endl;
    cout<<"                    .                     .          "<<endl;
    cout<<" R0 zN_z F_1,0N_z F_2,0N_z F_3,0N_z ...   F_N_nu,0N_z"<<endl;
    cout<<endl;
    cout<<" R1 z0   F_1,10   F_2,10   F_3,10         F_N_nu,10  "<<endl;
    cout<<" R1 z1   F_1,11   F_2,11   F_3,11   ...   F_N_nu,11  "<<endl;
    cout<<"(ECT)"<<endl<<endl;
    
    cout<<"When processing the input, the first line is ignored. Then"<<endl;
    cout<<"N_R  = Number of elements of R (i.e. R0, R1, R2, ..., RN_R)"<<endl;
    cout<<"N_z  = Number of elements of z (i.e. z0, z1, z2, ..., zN_z)"<<endl;
    cout<<"N_nu = Number of different F (alpha or jota)"<<endl<<endl;
    
    cout<<"N_nu is meant to be different frequencies, but you can assign each category to whatever you want (e.g.: different galaxies)"<<endl;
    cout<<"Nevertheless, you must make sure that each element of R and z are equispaced, or this code may not work."<<endl;
    cout<<"Furthermore, both files must have same dimensions of R and z, and have the same number of columns (that is, N_nu)."<<endl;
    cout<<"Ignoring the second requirement will lead to segfault"<<endl<<endl;
    
    cout<<"IC.txt is the external background mean intensity of your domain, and the initial condition when solving the equation of radiative transfer."<<endl;
    cout<<"Its file structure is a bit different:"<<endl;
    cout<<"I (units) => I_1 | I_2 | I_3 (...) I_N_nu "<<endl;
    cout<<"I_1 I_2 I_3 ... I_N_nu"<<endl;
    cout<<"First line is skipped and the second line includes the values for each nu, like the previous files"<<endl<<endl;
    
    cout<<"Output file will be the mean intensity, J_nu. The structure of this file is the same as alpha and jota files."<<endl<<endl<<endl;
    
    cout<<"Leaving without results..."<<endl;
    
    
    throw runtime_error("Bad syntax");
}

#endif
