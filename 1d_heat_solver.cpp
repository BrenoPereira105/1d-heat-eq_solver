#include<stdio.h>
#include<iostream>
#include<vector>
//this program uses cranck-nicolson method to solve 1d_heat_eq
//it generates implicit equations, what can be represented as a tridiagonal linear 
//system, and uses thomas algorithm to solve it. The Thomas algorithm 
//makes a LU decomposition of the tridiagonal original matrix by a generalized 
//form of gaussian elimination

void allocator(double ** &final_matrix, int xmax, int steps){

  final_matrix = new double*[xmax];
    for(int i = 0; i < xmax; ++i) {
    final_matrix[i] = new double[steps];
    }
}
void deallocator(double ** &final_matrix, int xmax, int steps){

    for(int i = 0; i < xmax; ++i) {
        delete [] final_matrix[i];
    }
    delete [] final_matrix;
    
}
void fill_boundary(double ** &final_matrix, int xmax, int steps, double u_0_t, double u_xmax_t){

    for(int steps_copy = 0; steps_copy < steps; steps_copy++){

            final_matrix[0][steps_copy] = u_0_t;
            final_matrix[xmax-1][steps_copy] = u_xmax_t;
        }
}
void fill_initial_condition(double ** &final_matrix, int xmax, double u_x_0){

    for(int xmax_copy = 0; xmax_copy < xmax; xmax_copy++){

        final_matrix[xmax_copy][0] = u_x_0;
    }
}
void fill_diagonals(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, int xmax, double r){

    a.push_back(0);
    for(int xmax_copy = 0; xmax_copy < xmax - 3; xmax_copy++){
        
        a.push_back(-r);
        b.push_back(2*(1+r));
        c.push_back(-r);
    }
    c.push_back(0);
    b.push_back(2*(1+r));
}

void solve(std::vector<double> a,std::vector<double> b,std::vector<double> c, std::vector<double> d, int n,
std::vector<double> &final_value) {
    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
    Finally we have a triangular matrix:
    |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report results
    Written by Keivan Moradi, 2014
    */
    std::vector<double> c_copy = c;
    std::vector<double> d_copy = d;
    //copy(c.begin(),c.begin()+c.size(),c_copy);
    //copy(d.begin(),d.begin()+d.size(),d_copy);
    
    n--; // since we start from x0 (not x1)
    c_copy[0] /= b[0];
    d_copy[0] /= b[0];

   
    for (int i = 1; i < n; i++) {
        c_copy[i] /= b[i] - a[i]*c_copy[i-1];
        d_copy[i] = (d_copy[i] - a[i]*d_copy[i-1]) / (b[i] - a[i]*c_copy[i-1]);
    }

    d_copy[n] = (d_copy[n] - a[n]*d_copy[n-1]) / (b[n] - a[n]*c_copy[n-1]);

    for (int i = n; i-- > 0;) {
        d_copy[i] -= c_copy[i]*d_copy[i+1];
    }

   
    final_value = d_copy;
    
}
void fill_d(int k,std::vector<double> &d, int xmax, double ** u, double r){

    if(d.size()!=0){
        d.clear();
    }

    for(int xmax_copy = 0; xmax_copy < xmax-2; xmax_copy++){

        if(xmax_copy == 0){

            d.push_back(2*u[xmax_copy+1][k]+r*(u[xmax_copy][k]
            -2*u[xmax_copy+1][k]+u[xmax_copy+2][k])+r*u[0][1]);
        }
        else if(xmax_copy == xmax-3){

            d.push_back(2*u[xmax_copy+1][k]+r*(u[xmax_copy][k]
            -2*u[xmax_copy+1][k]+u[xmax_copy+2][k])+r*u[xmax-1][1]);
        }
        else{
  
            d.push_back(2*u[xmax_copy+1][k]+r*(u[xmax_copy][k]
            -2*u[xmax_copy+1][k]+u[xmax_copy+2][k]));
        }
    }
}
void calculate_and_fill_final(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c,
std::vector<double> d, double r,int xmax, int steps, double ** &final_matrix){

    std::vector<double> new_u;
    int xmax_copy = 1;
    //int steps_copy = 0;
    
    for(int k = 0; k < steps-1; k++){

        fill_d(k,d,xmax,final_matrix,r);
        
    
        if(new_u.size()!=0){
            new_u.clear();
        }
        solve(a, b, c, d, xmax-2, new_u);
        
       
        //new_u.clear();
        //now filling the new values...
        while(xmax_copy < xmax-1){

            //std::cout<<new_u[xmax_copy-1]; 
            
            final_matrix[xmax_copy][k+1] = new_u[xmax_copy-1];
            xmax_copy++;
        }
        xmax_copy = 1;
    }
}
void printer(double ** u, int xmax, int steps){

    for(int steps_copy = 0; steps_copy<steps; steps_copy++){
        for(int xmax_copy = 0; xmax_copy<xmax; xmax_copy++){

            std::cout<<u[xmax_copy][steps_copy]<<" ";
        }
        std::cout<<"\n";
    }
}
int main(){

    //OBS: I will consider 'delta_x = delta_t = 1'

    //boundary conditions
    double u_0_t, u_xmax_t;
    //initial condition
    double u_x_0;
    //lentgh
    int xmax;
    //termal diffusion constant
    double alpha;
    //time steps
    int steps;
    

    std::cout<<"--set boundary conditions--\n";
    std::cout<<"\nu[0,t]:";
    std::cin>>u_0_t;
    std::cout<<"\nu[xmax,t]:";
    std::cin>>u_xmax_t;
    //lentgh
    std::cout<<"\n--set lentgh--\n";
    std::cout<<"\nxmax:";
    std::cin>>xmax;
    std::cout<<"\nsteps:";
    std::cin>>steps;
    std::cout<<"\n--set initial condition--\n";
    std::cout<<"\nu[x,0]:";
    std::cin>>u_x_0;
    std::cout<<"\n--set termal diffusion alpha--\n";
    std::cout<<"\nalpha:";
    std::cin>>alpha;
    
    //we will have to allocate memory for the final matrix (2d_array) representing time steps

    double ** final_matrix;
    
    allocator(final_matrix, xmax, steps);


    //first we fill boundary
    fill_boundary(final_matrix, xmax, steps, u_0_t, u_xmax_t);

    //now we fill initial condition
    fill_initial_condition(final_matrix, xmax, u_x_0);

    //first diagonal
    std::vector <double> a;
    //second diagonal
    std::vector <double> b;
    //third diagonal
    std::vector <double> c;
    //right side column
    std::vector <double> d;
    //r value from tridiagonal system
    double r = alpha/2;

    fill_diagonals(a,b,c, xmax, r);

    calculate_and_fill_final(a,b,c,d,r,xmax,steps,final_matrix);

    printer(final_matrix, xmax, steps);

    deallocator(final_matrix, xmax, steps);

    return 0;
}
