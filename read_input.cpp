#ifndef READ_INPUT_CPP_INCLUDED
#define READ_INPUT_CPP_INCLUDED

void skipline(ifstream& file)
{//What it says in the tin. Skips a whole line, or what it is left if you filled some data before.
    file.ignore(numeric_limits<streamsize>::max(), '\n');
}

const r_number* readExternalBackground(char* filename, const int N_nu)
{
    r_number* I0 = new r_number[N_nu];
    
    cout<<"Reading file "<<filename<<endl;
    ifstream file;
    file.open(filename);
    //Skip first line
    skipline(file);
    
    for(int nu=0;nu<N_nu;++nu)
    {
        file>>I0[nu];
    }
    
    file.close();
    
    return I0;
}

Table** readFiles(int argc, char** argv) 
{/*
   This routine reads the input you give to the exe. As
    > ./code.exe N_R N_z N_nu "ExternalBackground.txt" "file-Alpha.txt" "file-Jota.txt"
    First two parameters are the dimensions for R and z in BOTH files.
    Last two parameters are the files needed to fill each table, following the enum in main EXCEPT the last one
    So first file is 'Alpha' and second file is 'Jota'.
    ExternalBackground.txt is filled with the previous function.
    
    If the file does not respect N_R, N_z, code crashes with a segfault.
  */
  
  //Get dimensions
  const int Nx   = (int) (todouble(argv[1]));
  const int Ny   = (int) (todouble(argv[2]));
  const int N_nu = (int) (todouble(argv[3]));
  
  //Init data files
  Table** magnitude = new Table* [N_maps]; //This gets delete[] in main
  for(int variable=0;variable<N_maps;++variable)
  {
    magnitude[variable] = new Table[N_nu];
    for(int nu=0;nu<N_nu;++nu)
    {
        magnitude[variable][nu] = Table(Nx,Ny);
    }
  }
  
  //Get data
  for(int k=5;k<argc;++k) //Realize that the last magnitude table, corresponding to 'MeanIntensity' is not filled here (except X and Y), as it is the objective of the main function.
  {
      ifstream file;
      file.open(argv[k]);
      
      cout<<"Reading file "<<argv[k]<<endl;
      
      //Skip first line
      //file.ignore(numeric_limits<streamsize>::max(), '\n');
      skipline(file);
      
      r_number x,y;
      r_number prev_x;
      r_number z[N_nu];
      int i=0; //x index
      int j=0; //y index
      bool first_iteration = true;
      //float dummy;
      while(!file.eof())
      {
          
          //Read first two elements of the line
          file>>x;
          file>>y;
          
          //Get the indexes
          if(first_iteration){first_iteration=false;} //First iteration comes with i,j previously declared
          else if(x == prev_x){j++;}
          else{i++;j=0;}
          //Save the x for next iteration
          prev_x = x;
          
          for(int nu=0;nu<N_nu;++nu)
          {
              file>>z[nu];
              //Fill the corresponding table
              magnitude[k-5][nu].FillAll(x,y,z[nu],i,j);
              //Fill x and y of 'MeanIntensity' table, I'm using z=0.0 in order to know if this table is not filled properly later (also, it initializes J for the main loop).
              if(k==5)
              {
                magnitude[MeanIntensity][nu].FillAll(x,y,0.0,i,j);
                magnitude[Flux][nu].FillAll(x,y,0.0,i,j);
                magnitude[RadPressure][nu].FillAll(x,y,0.0,i,j);
              }
          }
          
          //Skip what remains of that line
          skipline(file);
          
          /*cout<<magnitude[k-4][0].X(i)<<" "<<magnitude[k-4][0].Y(j)<<" ";
          for(int nu=0;nu<N_nu;++nu)
          {
              cout<<magnitude[k-4][0].Z(i,j)<<" ";
          }
          cout<<" "<<i<<" "<<j<<endl;
          //cin>>dummy;
          */
      }
      
      file.close();
      //break;
  }
  //cout<<"Finished!"<<endl;
  
  return magnitude;
}

void deleteTables(Table** allTables)
{
  
  for(int variable=0;variable<N_maps;++variable)
  {
      delete[] allTables[variable];
  }
  delete[] allTables;
}

#endif
