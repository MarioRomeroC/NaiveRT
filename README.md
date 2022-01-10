# NaiveRT
Naive Radiative Transfer

This is a proof-of-concept of a 2d Radiative transfer code.
It solves the radiative transfer equation
dI/ds = j - a*I
and returns the mean intensity J = int(I*dOmega)/4pi for an axissymetric geometry
I uploaded to github in order to record this code.

This code was dropped in favour of Mixclask, located here: https://github.com/MarioRomeroC/Mixclask
If you want to run it, first read the define.h and change their options if you desire.
Next, compile it with
> g++ main.cpp -O2 -o NaiveRT.exe

Then, you can read the input instructions anytime by running the exe alone (./NaiveRT.exe).
In summary, you need three files, one for 'j', another for 'a' and a third one for the border initial condition.
First two files are 2D maps: First two columns are R and z, and following columns are the values for different wavelengths (one column = one particular wavelength)
Last file only contains a single row for the values of I(0) at the simulation-domain borders.

The output file is J, with the same structure as 'j' and 'a' files.
See the example folder to see the structure of these files.
