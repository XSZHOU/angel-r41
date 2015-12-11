#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>

#include <flens/flens.h>

//------------------------------------------------------------------------------

using namespace flens;
using namespace std;


//------------------------------------------------------------------------------

enum CartDir { North=0, East=1, South=2, West=3};

struct MpiCart
{
    MpiCart()
    {}
    
    MpiCart(int _numRows, int _numCols)
        : numRows(_numRows), numCols(_numCols)
    {
        int  dims[2];
        bool periods[2];

        // create mpi_cart:  2-dimensional, 2 x 2, non-periodic
        dims[0]    = numRows;
        dims[1]    = numCols;
        periods[0] = false;
        periods[1] = false;
        comm = MPI::COMM_WORLD.Create_cart(2, dims, periods, true);
        
        // find your position in the cart
        int rank = MPI::COMM_WORLD.Get_rank();
        int coords[2];
        comm.Get_coords(rank, 2, coords);
        row = coords[0];
        col = coords[1];
    
        // find your neighbors
        comm.Shift(0, 1, neighbors[North], neighbors[South]);
        comm.Shift(1, 1, neighbors[West],  neighbors[East]);
    }
    
    MpiCart(const MpiCart &rhs)
        : comm(rhs.comm),
          numRows(rhs.numRows), numCols(rhs.numCols),
          row(rhs.row), col(rhs.col)
    {
        for (int k=0; k<8; ++k) {
            neighbors[k] = rhs.neighbors[k];
        }
    }

    int
    cartRank(int row, int col)
    {
        if ((row<0) || (col<0) || (row==numRows) || (col==numCols)) {
            return MPI_PROC_NULL;
        }
        int coords[2];
        coords[0] = row;
        coords[1] = col;
        return comm.Get_cart_rank(coords);
    }
    
    MPI::Cartcomm comm;
    int numRows, numCols;
    int row, col;
    int neighbors[4];
};

struct DistributedVector
{
    typedef GeMatrix<FullStorage<double, RowMajor> > Grid;

    DistributedVector()
    {
    }

    DistributedVector(const MpiCart &_mpiCart, int rh)
        : mpiCart(_mpiCart), N(rh-1)
    {
        // global matrix has size (N+1) x (N+1)
        
        // compute m and n, such that local matrix has size (-1,m+1) x (-1,n+1)
        m = (N+1) / mpiCart.numRows - 1;
        n = (N+1) / mpiCart.numCols - 1;
        grid.resize(_(-1,m+1),_(-1,n+1));
        
        // m0, n0: row and col offset
        m0 = (m+1)*mpiCart.row;
        n0 = (n+1)*mpiCart.col;
        
        // firstRow, firstCol:
        // -> non-zero value of local grid is stored in (firstRow, firstCol)
        firstRow = (mpiCart.row==0) ? 1 : 0;
        firstCol = (mpiCart.col==0) ? 1 : 0;
                     
        // types to send rows and cols
        // MpiRow = MPI::DOUBLE.Create_contiguous(n+1); 
        MpiRow = MPI::DOUBLE.Create_contiguous(n+2); 
        MpiRow.Commit();
            
        MpiCol = MPI::DOUBLE.Create_vector(m+1, 1, grid.leadingDimension()); 
        MpiCol.Commit();
    }
    
    void
    updateBoundary()
    {
        MPI::Request req[8];
        MPI::Status  stat[8];
    
        req[ 0] = mpiCart.comm.Isend(&grid(0,n), 1, MpiCol, mpiCart.neighbors[East],  0);  
        req[ 1] = mpiCart.comm.Isend(&grid(m,0), 1, MpiRow, mpiCart.neighbors[South], 0);    
        req[ 2] = mpiCart.comm.Isend(&grid(0,0), 1, MpiCol, mpiCart.neighbors[West],  0);  
        req[ 3] = mpiCart.comm.Isend(&grid(0,0), 1, MpiRow, mpiCart.neighbors[North], 0);    

        req[ 4] = mpiCart.comm.Irecv(&grid(  0,n+1), 1, MpiCol, mpiCart.neighbors[East],  0);
        req[ 5] = mpiCart.comm.Irecv(&grid(m+1,  0), 1, MpiRow, mpiCart.neighbors[South], 0);
        req[ 6] = mpiCart.comm.Irecv(&grid(  0, -1), 1, MpiCol, mpiCart.neighbors[West],  0);
        req[ 7] = mpiCart.comm.Irecv(&grid( -1,  0), 1, MpiRow, mpiCart.neighbors[North], 0);
        
        MPI::Request::Waitall(8, req, stat);
    }
    
    MpiCart         mpiCart;
    int             N;
    int             m, n;
    int             m0, n0;
    int             firstRow, firstCol;
    Grid            grid;
    MPI::Datatype   MpiRow, MpiCol;
};

double
norm(const DistributedVector &x)
{
    double h = 1./(x.N+1);
    
    // compute (frobenius) norm
    // -> local norm
    double r = 0;
    for (int i=0; i<=x.m; ++i) {
        for (int j=0; j<=x.n; ++j) {
            r += x.grid(i,j)*x.grid(i,j);
        }
    }

    // -> global norm
    double R;
    x.mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    return sqrt(h*R);
}

double
sum(const DistributedVector &x)
{
    // -> local sum
    double r = 0;
    for (int i=0; i<=x.m; ++i) {
        for (int j=0; j<=x.n; ++j) {
            r += x.grid(i,j);
        }
    }

    // -> global sum
    double R;
    x.mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    return R;
}

void
setRHS(DistributedVector &f)
{
    double h = 1./(f.N+1);
    
    for (int i=0; i<=f.m; ++i) {
        for (int j=0; j<=f.n; ++j) {
            double x = (i+f.m0+0.5)*h;
            double y = (j+f.n0+0.5)*h;
            
            f.grid(i,j) =  2*x+2*y-2;
        }
    }
}

void
setExactSolution(DistributedVector &e)
{
    double h = 1./(e.N+1);
    
    for (int i=0; i<=e.m; ++i) {
        for (int j=0; j<=e.n; ++j) {
            double x = (i+e.m0+0.5)*h;
            double y = (j+e.n0+0.5)*h;
            
            e.grid(i,j) =  x*x/2 - x*x*x/3 + y*y/2 - y*y*y/3 - 1./6;
        }
    }
}

void
gaussSeidelRedBlack(double omega, const DistributedVector &f, DistributedVector &u)
{
    const DistributedVector::Grid &F = f.grid;
    DistributedVector::Grid       &U = u.grid;
    
    double c         = omega/4.;
    double hh        = 1./((u.N+1)*(u.N+1));

    if (u.mpiCart.row==0) {
        U(-1,_) = U(0,_);
    }
    if (u.mpiCart.col==0) {
        U(_,-1) = U(_,0);
    }
    if (u.mpiCart.row==u.mpiCart.numRows-1) {
        U(u.m+1,_) = U(u.m,_);
    }
    if (u.mpiCart.col==u.mpiCart.numCols-1) {
        U(_,u.n+1) = U(_,u.n);
    }

    // red nodes
    for (int i=0; i<=u.m; ++i) {
        for (int j=0; j<=u.n; ++j) {
            
            if ((i+j)%2==1) {
                continue;
            }
             
            U(i,j) = (1-omega)*U(i,j)
                   + c*(F(i,j)*hh + U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1));        
        }
    }

    u.updateBoundary();

    // black nodes

    for (int i=0; i<=u.m; ++i) {
        for (int j=0; j<=u.n; ++j) {
            
            if ((i+j)%2==0) {
                continue;
            }
             
            U(i,j) = (1-omega)*U(i,j)
                   + c*(F(i,j)*hh + U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1));        
        }
    }

    u.updateBoundary();
    
    if (u.mpiCart.row==0) {
        U(-1,_) = U(0,_);
    }
    if (u.mpiCart.col==0) {
        U(_,-1) = U(_,0);
    }
    if (u.mpiCart.row==u.mpiCart.numRows-1) {
        U(u.m+1,_) = U(u.m,_);
    }
    if (u.mpiCart.col==u.mpiCart.numCols-1) {
        U(_,u.n+1) = U(_,u.n);
    }

    int rhh = (u.N+1)*(u.N+1);
    double s = sum(u)/rhh;
    
    for (int i=0; i<=u.m; ++i) {
        for (int j=0; j<=u.n; ++j) {
            U(i,j) -= s;
        }
    }
}

void
residualPoisson(const DistributedVector &f,
                DistributedVector &u,
                DistributedVector &r)
{
    int rhh = (f.N+1)*(f.N+1);
 
    const DistributedVector::Grid &F = f.grid;
    DistributedVector::Grid &U = u.grid;
    DistributedVector::Grid       &R = r.grid;
    
    if (u.mpiCart.row==0) {
        U(-1,_) = U(0,_);
    }
    if (u.mpiCart.col==0) {
        U(_,-1) = U(_,0);
    }
    if (u.mpiCart.row==u.mpiCart.numRows-1) {
        U(u.m+1,_) = U(u.m,_);
    }
    if (u.mpiCart.col==u.mpiCart.numCols-1) {
        U(_,u.n+1) = U(_,u.n);
    }
    
    for (int i=0; i<=f.m; ++i) {
        for (int j=0; j<=f.n; ++j) {
            R(i,j) = F(i,j) - rhh*(4*U(i,j) -U(i+1,j) -U(i-1,j) -U(i,j+1) -U(i,j-1));
        }
    }
}

double
errorNorm(const DistributedVector &e, const DistributedVector &u)
{
    double h = 1./(e.N+1);
    
    double r = 0;
    for (int i=0; i<=e.m; ++i) {
        for (int j=0; j<=e.n; ++j) {
            r += (e.grid(i,j)-u.grid(i,j))*(e.grid(i,j)-u.grid(i,j));
        }
    }
    
    double R;
    u.mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    return sqrt(h*h*R);
}

void
printStat(int k,
          const DistributedVector &r,
          const DistributedVector &u)
{
    double rNorm = norm(r);

    if ((r.mpiCart.row==0) && (r.mpiCart.col==0)) {
        cout.width(3);
        cout << k << ") | ";
    
        cout.precision(12);
        cout.setf(std::ios::fixed);
        cout.width(17);
        cout << rNorm << endl;
    }
}

// f_c = R*r
void
restriction(const DistributedVector &r, DistributedVector &f_c)
{
    const DistributedVector::Grid &R   = r.grid;
    DistributedVector::Grid       &F_C = f_c.grid;
    
    for (int i=0, I=0; i<=f_c.m; ++i, I+=2) {
        for (int j=0, J=0; j<=f_c.n; ++j, J+=2) {
            F_C(i,j) = (R(I,J)+R(I,J+1)+R(I+1,J)+R(I+1,J+1))/4;
        }
    }
}

// u = u + P*u_c
void
prolongation(const DistributedVector &u_c, DistributedVector &u)
{
      
    const DistributedVector::Grid &U_C   = u_c.grid;
    DistributedVector::Grid       &U     = u.grid;

    for (int i=0, I=0; i<=u_c.m+1; ++i, I+=2) {
        for (int j=0, J=0; j<=u_c.n+1; ++j, J+=2) {
            U(  I,  J) += 0.75*(0.75*U_C(i,j)+0.25*U_C(i-1,j)) + 0.25*(0.75*U_C(i,j-1)+0.25*U_C(i-1,j-1));
            U(I-1,  J) += 0.75*(0.25*U_C(i,j)+0.75*U_C(i-1,j)) + 0.25*(0.25*U_C(i,j-1)+0.75*U_C(i-1,j-1));
            U(  I,J-1) += 0.25*(0.75*U_C(i,j)+0.25*U_C(i-1,j)) + 0.75*(0.75*U_C(i,j-1)+0.25*U_C(i-1,j-1));
            U(I-1,J-1) += 0.25*(0.25*U_C(i,j)+0.75*U_C(i-1,j)) + 0.75*(0.25*U_C(i,j-1)+0.75*U_C(i-1,j-1));
        }
    }
    u.updateBoundary();
}

void
write(const DistributedVector &v, const string &name)
{
    ostringstream s;
    s << name << "_local_" << v.mpiCart.row << "_" << v.mpiCart.col << ".dat";
    ofstream out(s.str().c_str());
    out << v.N << " " << v.N << endl;
    out << v.m-v.firstRow+1 << " " << v.n-v.firstCol+1<< endl; 
    out << v.grid(_(v.firstRow, v.m),_(v.firstCol, v.n));
}

void
read(DistributedVector::Grid &grid, const string &name, int numProcX, int numProcY)
{
    int row = 1;
    for (int px=0; px<numProcX; ++px) {
        int col = 1;
	    for (int py=0; py<numProcY; ++py) {
	        ostringstream s;
	        s << name << "_local_" << px << "_" << py << ".dat";
	        ifstream in(s.str().c_str());
	        int N, M;
            in >> N >> M;
	        if ((px==0) && (py==0)) {
		        grid.resize(N+2,M+2,0,0);
	        }
            int r; int c;
            in >> r; in >> c;
	        for (int i=0; i<r; ++i) {
                for (int j=0; j<c; ++j) {
                    in >> grid(row+i,col+j);
		        }
	        }
            col += c;
	        if (py==numProcY-1) {
	            row += r;
            }
	    }
    }
}

//------------------------------------------------------------------------------

const long l = 8;
const int  n = (1<<l) - 1;  // n = 2^l -1

DistributedVector  u[l+1], f[l+1], r[l+1];

void
initMultigrid(const MpiCart &mpiCart)
{
    typedef DistributedVector Vec;

    int N = n;
    for (int i=l; i>=1; --i, N/=2) {
        u[i] = Vec(mpiCart, N+1);
        f[i] = Vec(mpiCart, N+1);
        r[i] = Vec(mpiCart, N+1);
    }
    for (int i=l; i>=3; --i) {
        restriction(f[i], f[i-1]);
    }
}

void
multigrid(int l, int v1, int v2)
{
    if (l<=2) {
        for (int k=1; k<=10; ++k) {
            gaussSeidelRedBlack(1., f[l], u[l]);
        }
    } else {
        for (int v=1; v<=v1; ++v) {
            gaussSeidelRedBlack(1., f[l], u[l]);
        }
        residualPoisson(f[l], u[l], r[l]);

        restriction(r[l], f[l-1]);


        u[l-1].grid = 0;
        multigrid(l-1,v1,v2);
        
        prolongation(u[l-1], u[l]);
        
        
        for (int v=1; v<=v2; ++v) {
            gaussSeidelRedBlack(1., f[l], u[l]);
        }
    }
}

void
fullMultigrid()
{
    for (int k=1; k<=1000; ++k) {
        gaussSeidelRedBlack(1., f[2], u[2]);
    }
    
    for (int i=2; i<l; ++i) {
        prolongation(u[i], u[i+1]);
        multigrid(i+1, 1, 1);
    }
}

//------------------------------------------------------------------------------

typedef GeMatrix<FullStorage<double, RowMajor> > Mat;

// physics
double h = 1./n,
       re = 1000,
       pr = 7,
       alpha = 1./(re*pr),
       Gamma = 0.9,
       gx = 0.,
       gy = -0.00981;


// for time steps
double delta = 0.5,
       uMax,
       vMax;

void
setBoundaryCondition(Mat &U, Mat &V)
{
	U(  0,_) = 0.;
	U(n+1,_) = 0.;

	V(_,  0) = 0.;
	V(_,n+1) = 0.; 

	V( -1,_) = -V(0,_);
	V(n+1,_) = -V(n,_);

	double u0 = 1;

	U(_,-1) = -U(_,0);
	
	U(_,n+1)  = 2*u0;
	U(_,n+1) -= U(_,n);
}

double
computeTimeStep(const Mat &U, const Mat &V)
{
	double f1, f2, f3;

    uMax = 0.01;
    vMax = 0.01;

	for (int i=0; i<=n+1; ++i) {
		for (int j=-1; j<=n+1; ++j) {
			if (fabs(U(i,j)) > uMax) {
				uMax = fabs(U(i,j));
			}
        }
    }
	for (int i=-1; i<=n+1; ++i) {
		for (int j=0; j<=n+1; ++j) {
			if (fabs(V(i,j)) > vMax) {
				vMax = fabs(V(i,j));
			}
		}
	}
	f1 = re*pr*h*h/4;
	f2 = h/uMax;
	f3 = h/vMax;
	
	return delta*std::min(std::min(f1, f2), f3);
}

double
sqr(double x)
{
	return x*x;
}

void
computeForce(const Mat &U, const Mat &V, double dt, Mat &Fx, Mat &Fy)
{
    double rh = (1+n),
           rhh = rh*rh;
	
	for (int i=1; i<=n; ++i) {
		for (int j=0; j<=n; ++j) {
			
			double duu = sqr(U(i,j)   + U(i+1,j))
			           - sqr(U(i-1,j) + U(i,j))
			    +Gamma*( fabs(  U(i,j)+U(i+1,j))*(  U(i,j)-U(i+1,j))
			            -fabs(U(i-1,j)+  U(i,j))*(U(i-1,j)-U(i,j)));
			duu *= rh/4;

			double duv = (V(i,j+1)+V(i-1,j+1))*(  U(i,j)+U(i,j+1))
                       - (V(i,  j)+V(i-1,  j))*(U(i,j-1)+U(i,j))
			    +Gamma*( fabs(V(i,j+1)+V(i-1,j+1))*(U(i,j)  -U(i,j+1))
			            -fabs(V(i,j)  +V(i-1,j))  *(U(i,j-1)-U(i,j)));
			duv *= rh/4;
			
			double ddux = (U(i+1,j) - 2*U(i,j) + U(i-1,j))*rhh;
            double dduy = (U(i,j+1) - 2*U(i,j) + U(i,j-1))*rhh;               		
			
			Fx(i,j) = U(i,j) + dt*((ddux+dduy)/re - duu -duv +gx);
		}
	}

	for (int i=0; i<=n; ++i) {
		for (int j=1; j<=n; ++j) {
			double dvv = sqr(V(i,j)   + V(i,j+1))
			           - sqr(V(i,j-1) + V(i,j))
			    +Gamma*( fabs(V(i,j)+V(i,j+1))*(V(i,j)-V(i,j+1))
			            -fabs(V(i,j-1)+V(i,j))*(V(i,j-1)-V(i,j)));
			dvv *= rh/4;

			double duv = (U(i+1,j)+U(i+1,j-1))*(V(  i,j)+V(i+1,j))
			           - (U(  i,j)+U(  i,j-1))*(V(i-1,j)+V(  i,j))
			    +Gamma*( fabs(U(i+1,j)+U(i+1,j-1))*(V(  i,j)-V(i+1,j))
			            -fabs(U(  i,j)+U(  i,j-1))*(V(i-1,j)-V(  i,j)));
			duv *= rh/4;
			
			double ddvx = (V(i+1,j)-2*V(i,j)+V(i-1,j))*rhh;
			
			double ddvy = (V(i,j+1)-2*V(i,j)+V(i,j-1))*rhh;
			
			Fy(i,j) = V(i,j) + dt*((ddvx + ddvy)/re - duv -dvv);
		}
	}
	
	Fx(  0,_) = U(  0,_);
	Fx(n+1,_) = U(n+1,_);
	
	Fy(_,  0) = V(_,  0);
	Fy(_,n+1) = V(_,n+1);
}

void
computePoissonRhs(const Mat &Fx, const Mat &Fy, double dt, Mat &Rhs)
{
    for (int i=0; i<=n; ++i) {
	    for (int j=0; j<=n; ++j) {
            Rhs(i,j) =-((Fx(i+1,j)-Fx(i,j))/h+(Fy(i,j+1)-Fy(i,j))/h)/dt;
        }
    }
}

void
solvePoisson(const Mat &Rhs, Mat &P)
{
    f[l].grid = Rhs;
    
    
    residualPoisson(f[l], u[l], r[l]);  
    //printStat(0, r[l], u[l]);

    fullMultigrid();
    residualPoisson(f[l], u[l], r[l]);  
    //printStat(1, r[l], u[l]);
    
    int it=0;
    while (norm(r[l])>0.01) {
        multigrid(l,1,1);
        
        residualPoisson(f[l], u[l], r[l]); 
        printStat(1, r[l], u[l]);
        
        ++it;
        if (it>20) {
            break;
        }
    }
    printStat(1, r[l], u[l]);
    P = u[l].grid;
    
    P(_,-1) = P(_,0);
    P(-1,_) = P(0,_);

    P(_,n+1) = P(_,n);
    P(n+1,_) = P(n,_);

}

void
computeVelocity(const Mat &Fx, const Mat &Fy, const Mat &P, double dt,
                Mat &U, Mat &V)
{
    for (int i=1; i<=n; ++i) {
		for (int j=0; j<=n; ++j) {
			U(i,j) = Fx(i,j) -  dt*(P(i,j)-P(i-1,j))/h;
		}
	}
	for (int i=0; i<=n; ++i) {
		for (int j=1; j<=n; ++j) {
			V(i,j) = Fy(i,j) - dt*(P(i,j)-P(i,j-1))/h;
		}
	}
}

void
writeVelocity(std::ostream &out, const Mat &U, const Mat &V)
{
	int xStep = 3;
	int yStep = 3;
	
	for (int j=0; j<n; j+=yStep) {
		for (int i=0; i<n; i+=xStep) {
			float vx = (U(i+1,  j) + U(i,j))*0.5;
			float vy = (V(  i,j+1) + V(i,j))*0.5;
			
			float nrm = sqrt(vx*vx+vy*vy);
			
			vx /= nrm;
			vy /= nrm;
			
			out << setw(4) << (i+0.5)*h << " " 
			    << setw(4) << (j+0.5)*h << " "
			    << setw(15) << vx << " " 
			    << setw(15) << vy << std::endl;
		} 
	}
}

int
main(int argc, char *argv[])
{
    MPI::Init(argc, argv);
    
    MpiCart mpiCart(1,1);

    initMultigrid(mpiCart);
    
    
    DistributedVector Rhs(mpiCart, n+1);
    
    Mat U(_(-1,n+1), _(-1,n+1)),
        V(_(-1,n+1), _(-1,n+1)),
        P(_(-1,n+1), _(-1,n+1)),
        Fx(_(-1,n+1), _(-1,n+1)),
        Fy(_(-1,n+1), _(-1,n+1));
        
    double dt, t, storeTime=0;
    
    int dataSet = 0;
        
/*
    setBoundaryCondition(U, V);
    
    cout << "U = " << U << endl;
    cout << "V = " << V << endl;
    
    t = 0;
    dt = computeTimeStep(U, V);
    
    cout << "dt = " << dt << endl;
    
   
	computeForce(U, V, dt, Fx, Fy);
	
	cout << "Fx = " << Fx << endl;
	cout << "Fy = " << Fy << endl;
	
	computePoissonRhs(Fx, Fy, dt, Rhs.grid);
			
	cout << "Rhs = " << Rhs.grid << endl;
	
	solvePoisson(Rhs.grid, P);

	cout << "P = " << P << endl;

    computeVelocity(Fx, Fy, P, dt, U, V);
	cout << "U = " << U << endl;
    cout << "V = " << V << endl;

    cout << "----------------" << endl;
    cout << "dt = " << dt << endl;

    setBoundaryCondition(U, V);

    computeForce(U, V, dt, Fx, Fy);
	
	cout << "Fx = " << Fx << endl;
	cout << "Fy = " << Fy << endl;
	
	Rhs.grid = 0;
	
	computePoissonRhs(Fx, Fy, dt, Rhs.grid);
			
	cout << "Rhs = " << Rhs.grid << endl;
	
	solvePoisson(Rhs.grid, P);

	cout << "P = " << P << endl;

    computeVelocity(Fx, Fy, P, dt, U, V);
	cout << "U = " << U << endl;
    cout << "V = " << V << endl;

    
    ostringstream s;
	s << "data/velocity" << setw(8) << setfill('0') << dataSet++ << ".dat";
	ofstream out(s.str().c_str());
	writeVelocity(out, U, V);
*/

    setBoundaryCondition(U, V);
    t = 0;
    dt = computeTimeStep(U, V);
    while (t<10) {
        cout << "t = " << t << endl;
        computeForce(U, V, dt, Fx, Fy);
        computePoissonRhs(Fx, Fy, dt, Rhs.grid);
        solvePoisson(Rhs.grid, P);
        computeVelocity(Fx, Fy, P, dt, U, V);
		t += dt;
        setBoundaryCondition(U, V);
        dt = computeTimeStep(U, V);
        
        if (t>storeTime) {
			storeTime += .1;

			ostringstream s;
			s << "data/velocity" << setw(8) << setfill('0') << dataSet++ << ".dat";
			ofstream out(s.str().c_str());
			writeVelocity(out, U, V);
		}
    }
    
    
    MPI::Finalize();
}
