#ifndef WRITE_OUTPUT_CPP_INCLUDED
#define WRITE_OUTPUT_CPP_INCLUDED

//Here there are algorithms to write multiple tables in one single file (unlike Table::write_table)

/*string header(r_number* array, int N_array, string array_units="", string magnitude="", string magnitude_units="")
{//Ideally, array have all frequencies used in main
    string line = " R (pc) | z (pc) || ";
    line += magnitude+"("+magnitude_units+"); nu("+array_units+") => ";
    for(int nu=0;nu<N_array;++nu)
    {
        line += tostring(nu)+" | ";
    }
    return line;
    //Result should be, if magnitude = "J", magnitude_units "erg/cm2/s/Hz", array_units="Hz" and array = [1,3,5,7]
    // R (pc) | z (pc) || J(erg/cm2/s/Hz); nu(Hz) => 1 | 3 | 5 | 7
}*/

void writeFile(Table* allTables, int N_tables, string nameFile, string firstLine="")
{
    //Ideally, allTables should be an array of the same magnitude (alpha, meanIntensity...) where each element is a different frequency.
    //All tables must have same dimensions (that is, same number of X elements for each table, as well as with Y elements)
    
    cout<<"Creating file "<<nameFile<<endl;
    ofstream file;
    file.open(nameFile.c_str());
    
    int Nx = allTables[0].Xdim();
    int Ny = allTables[0].Ydim();
    
    file<<firstLine<<endl;
    for(int i=0;i<Nx;++i)
    {
        for(int j=0;j<Ny;++j)
        {
            file<<allTables[0].X(i)<<" "<<allTables[0].Y(j)<<" ";
            for(int nu=0;nu<N_tables;++nu)
            {
                file<<allTables[nu].Z(i,j)<<" ";
            }
            file<<endl;
        }
        file<<endl;
    }
    
    file.close();
}

#endif
