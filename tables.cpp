#ifndef TABLES_CPP_INCLUDED
#define TABLES_CPP_INCLUDED
/*
 This is my bilinear class, BUT without any for loop except when you construct them.
 
 It is very important to have elements of x and y sorted and equispaced beforehand.
 Otherwise, this class is only useful to collect and print data in tables
 */

using namespace std;

class Table
{
    private:
        r_number *x; //x points
        int Nx;    //Number of x points

        r_number *y; //y points
        int Ny;     //Number of y points

        r_number **z; //z(x,y)

        r_number bilinear_interpolation(r_number X, r_number Y, int i1, int i2, int j1, int j2){
            //Create the needed values:
            r_number x1 = x[i1];
            r_number x2 = x[i2];

            r_number y1 = y[j1];
            r_number y2 = y[j2];

            r_number z11 = z[i1][j1];
            r_number z12 = z[i1][j2];
            r_number z21 = z[i2][j1];
            r_number z22 = z[i2][j2];

            //Vile copy of wikipedia!
            r_number den = (y2-y1)*(x2-x1);
            r_number num = (y2-Y)*(x2-X)*z11;
            num += (y2-Y)*(X-x1)*z21;
            num += (Y-y1)*(x2-X)*z12;
            num += (Y-y1)*(X-x1)*z22;

            //cout<<num/den<<endl;
            return num/den;
        }

        r_number bilog_interpolation(r_number X, r_number Y, int i1, int i2, int j1, int j2){
            //Create the needed values:
            r_number x1 = (x[i1]);
            r_number x2 = (x[i2]);

            r_number y1 = (y[j1]);
            r_number y2 = (y[j2]);
            
            r_number z11 = std::log10(z[i1][j1]);
            r_number z12 = std::log10(z[i1][j2]);
            r_number z21 = std::log10(z[i2][j1]);
            r_number z22 = std::log10(z[i2][j2]);

            //Vile copy of wikipedia!
            r_number den = (y2-y1)*(x2-x1);
            r_number num = (y2-Y)*(x2-X)*z11;
            num += (y2-Y)*(X-x1)*z21;
            num += (Y-y1)*(x2-X)*z12;
            num += (Y-y1)*(X-x1)*z22;

            return std::pow(10.,num/den);
            
        }
        
        void find_indexes(float k, int& k1, int& k2, bool& same)
        {
            k1 = floor(abs(k)); //Absolute value is to count mirror boundary conditions
            k2 = ceil(abs(k));
            if( k1 == k2 ){same = true;}
            
        }
        
    public:
        Table()
        {
            x = nullptr;
            y = nullptr;
            z = nullptr;
        }
        //Legacy constructor
        Table(r_number* X, r_number* Y, r_number** Z, int n_x, int n_y) //Remember, sorted arrays!
        {
            Nx = n_x;
            Ny = n_y;
            x = new r_number[Nx];
            y = new r_number[Ny];
            z = new r_number*[Nx];
            
            for(int i=0;i<Nx;++i)
            {
                x[i] = X[i];
                z[i] = new r_number[Ny];
            }
            for(int i=0;i<Nx;++i)
            {
                for(int j=0;j<Ny;++j)
                {
                    if(i==0){ y[j] = Y[j]; }
                    z[i][j] = Z[i][j];
                }
            }
        }
        
        //This constructor, combined with 'fill' methods, avoids you to do two loops (you don't need to have *x,*y and **z first)
        Table(int n_x, int n_y)
        {
            Nx = n_x;
            Ny = n_y;
            x = new r_number[Nx];
            y = new r_number[Ny];
            z = new r_number*[Nx];
            for(int i=0;i<Nx;++i){ z[i] = new r_number[Ny]; }
        }
        //Rule of 3
        ~Table(){
            if( x != nullptr){
                delete[] x;
                x = nullptr;
            }
            if( y != nullptr){
                delete[] y;
                y = nullptr;
            }
            if( z != nullptr){
                for(int i=0; i<Nx;i++){
                    delete[] z[i];
                }
                delete[] z;
                z = nullptr;
            }
        }
        
        Table(const Table& other){
        //C++ rule of 3: If you have pointer members to dynamic memory, you need a destructor, a copy constructor, and assigment operator.
        //This is to avoid shallow copy (copy of pointer members points to the same address as the original) and awful deletes

            //Copy x
            x = new r_number[other.Nx];
            Nx = other.Nx;
            for(int i=0;i<Nx;i++){
                x[i] = other.x[i];
            }

            //Copy y
            y = new r_number[other.Ny];
            Ny = other.Ny;
            for(int j=0;j<Ny;j++){
                y[j] = other.y[j];
            }

            //Copy z
            z = new r_number*[Nx];
            for(int i=0;i<Nx;i++){
                z[i] = new r_number[Ny];
            }
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    z[i][j] = other.z[i][j];
                }
            }
        }

        Table& operator=(const Table& other){
            //This represent the operation *this = other. WE ARE OVERRIDING *this, with potential memory leaks if we are not careful!

            //Override *this first (but not assign to null)
            if( x != nullptr){delete[] x;}
            if( y != nullptr){delete[] y;}
            if( z != nullptr){
                for(int i=0;i<Nx;i++){ delete[] z[i]; }
                delete[] z;
            }

            Nx = other.Nx;
            x = new r_number[Nx];
            for(int i=0;i<Nx;i++){
                x[i] = other.x[i];
            }

            Ny = other.Ny;
            y = new r_number[Ny];
            for(int j=0;j<Ny;j++){
                y[j] = other.y[j];
            }

            z = new r_number*[Nx];
            for(int i=0;i<Nx;i++){
                z[i] = new r_number[Ny];
            }
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny;j++){
                    z[i][j] = other.z[i][j];
                }
            }

            return *this;
        }
        //FILL METHODS
        void FillX(r_number X, int index)
        {
            x[index] = X;
        }
        void FillY(r_number Y, int index)
        {
            y[index] = Y;
        }
        void FillZ(r_number Z, int indexX, int indexY)
        {
            z[indexX][indexY] = Z;
        }
        void FillAll(r_number X,r_number Y, r_number Z, int indexX, int indexY)
        {
            FillX(X,indexX);
            FillY(Y,indexY);
            FillZ(Z,indexX,indexY);
        }
        void SumZ(r_number Z, int indexX, int indexY)
        {
            z[indexX][indexY] += Z;
        }
        
        //ACCESS METHODS
        r_number X(int index){ return x[index]; }
        r_number Y(int index){ return y[index]; }
        r_number Z(int indexX, int indexY){ return z[indexX][indexY]; }
        int Xdim(){return Nx;}
        int Ydim(){return Ny;}
        
        //GET VALUES
        r_number get_value(r_number X, r_number Y)
        {
            //Find the indexes with divisions
            float i = (Nx-1.0)*(X-x[0])/(x[Nx-1]-x[0]);
            float j = (Ny-1.0)*(Y-y[0])/(y[Ny-1]-y[0]);
            
            //cout<<X<<" "<<Y<<endl;
            
            if( abs(i) > Nx-1 || abs(j) > Ny-1 ) //If you are out of bounds (X > x[Nx-1], idem with Y), give 0.0. If you're negative, see find_indexes for the other boundary condition
            {
                //cout<<"Out of bounds "<<Nx-1<<" "<<Ny-1<<endl;
                return 0.0;
            }
            else
            {
                bool same = false;
                int i1,i2;
                find_indexes(i, i1, i2, same);
                if(same){ i2 = (i1 == Nx-1 ) ? i2-1 : i2+1; } //In paper, this does not matter. In code, it means segfault
                
                same = false;
                int j1,j2;
                find_indexes(j,j1,j2,same);
                if(same){ j2 = (j1 == Ny-1 ) ? j2-1 : j2+1; }
                
                //cout<<i1<<" "<<i2<<" , "<<j1<<" "<<j2<<endl;
            
                //Get the interpolation
                return bilinear_interpolation(X,Y,i1,i2,j1,j2);
            }
            
            
        }
        
        /*r_number get_log(r_number X, r_number Y)
        {
            //Find the indexes with divisions
            float i = Nx*(X-x[0])/(x[Nx-1]-x[0]);
            float j = Ny*(Y-y[0])/(y[Ny-1]-y[0]);
            if( abs(i) > Nx-1 && abs(j) > Ny-1 ) //If you are out of bounds (X > x[Nx-1], idem with Y), give 0.0. If you're negative, see find_indexes for the other boundary condition
            {
                return 0.0;
            }
            else
            {
                int i1,i2;
                find_indexes(i, i1, i2);
                
                int j1,j2;
                find_indexes(j,j1,j2);
            
                //Get the interpolation
                return bilog_interpolation(X,Y,i1,i2,j1,j2);
            }
        }*/
        
        //PRINT OUTPUTS
        //This function is used for testing purposes
        void write_table(string namefile, string firstline = "")
        {
            cout<<"Created file "<<namefile<<endl;
            
            ofstream file;
            file.open(namefile.c_str());
            
            file<<firstline<<endl;
            for(int i=0;i<Nx;++i)
            {
                for(int j=0;j<Ny;++j)
                {
                    file<<x[i]<<" "<<y[j]<<" "<<z[i][j]<<endl;
                }
                file<<endl; //This is to indicate a change in x. Helps gnuplot.
            }
            
            file.close();
        }
};


#endif
