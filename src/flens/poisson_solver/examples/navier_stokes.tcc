namespace flens {

//------------------------------------------------------------------------------

template <typename U, typename V>
void
setBoundaryCondition(int rh, GeMatrix<U> &u, GeMatrix<V> &v, double u0)
{
    int N = rh - 1;

    u(0,_) = 0; u(N+1,_) = 0;
    v(_,0) = 0; v(_,N+1) = 0;

    u(_,-1)  = -u(_,0);
    u(_,N+1) = 2*u0;
    u(_,N+1) -= u(_,N);

    v( -1,_) = -v(0,_);
    v(N+1,_) = -v(N,_);
}

//------------------------------------------------------------------------------

template <typename TMP>
void
setTemperatureBoundaryCondition(GeMatrix<TMP> &T, double t0)
{
    int i0 = T.firstRow()+1;
    int i1 = T.lastRow()-1;
    int j0 = T.firstCol()+1;
    int j1 = T.lastCol()-1;

    for (int i=T.firstRow()+1; i<=T.lastRow()-1; ++i) {
        double x = double(i)/T.numRows();
        T(i,j1+1) = t0*(1-x);
    }

    T(_,j1+1) -= T(_,j1);

    T(_,j0-1) = T(_,j0);
    T(i0-1,_) = T(i0,_);
    T(i1+1,_) = T(i1,_);

}

//------------------------------------------------------------------------------

template <typename U, typename V>
double
computeTimeStep(int rh, const GeMatrix<U> &u, const GeMatrix<V> &v,
                double re, double pr, double delta)
{
    double uMax = 0;
    for (int i=u.firstRow(); i<=u.lastRow(); ++i) {
        for (int j=u.firstCol(); j<=u.lastCol(); ++j) {
            if (std::abs(u(i,j)) > uMax) {
                uMax = std::abs(u(i,j));
            }
        }
    }

    double vMax = 0;
    for (int i=v.firstRow(); i<=v.lastRow(); ++i) {
        for (int j=v.firstCol(); j<=v.lastCol(); ++j) {
            if (std::abs(v(i,j)) > vMax) {
                vMax = std::abs(v(i,j));
            }
        }
    }
    std::cout << "uMax = " << uMax << ", vMax = " << vMax << std::endl;

    uMax = std::max(uMax,0.01);
    vMax = std::max(vMax,0.01);

    double h = 1./rh;
    double f1, f2, f3;

    f1 = 0.25*re*pr*h*h/4;
    f2 = h/uMax;
    f3 = h/vMax;

    return delta*std::min(std::min(f1, f2), f3);
}

//------------------------------------------------------------------------------

template <typename T>
T
sqr(T x)
{
    return x*x;
}

template <typename U, typename V, typename FX, typename FY>
void
computeForce(int rh, const GeMatrix<U> &u, const GeMatrix<V> &v,
             double dt, double gamma, double re, double gx, double gy,
             GeMatrix<FX> &fx, GeMatrix<FY> &fy)
{
    double rhh = rh*rh;

    for (int i=u.firstRow()+1; i<=u.lastRow()-1; ++i) {
        for (int j=u.firstCol()+1; j<=u.lastCol()-1; ++j) {
            double duu = sqr(u(  i,j) + u(i+1,j))
                       - sqr(u(i-1,j) + u(  i,j))
                +gamma*( std::abs(u(  i,j)+u(i+1,j))*(u(  i,j)-u(i+1,j))
                        -std::abs(u(i-1,j)+u(  i,j))*(u(i-1,j)-u(  i,j)));
            duu *= rh/4;

            double duv = (v(i,j+1)+v(i-1,j+1))*(u(i,  j)+u(i,j+1))
                       - (v(i,  j)+v(i-1,  j))*(u(i,j-1)+u(i,  j))
                +gamma*( std::abs(v(i,j+1)+v(i-1,j+1))*(u(i,  j)-u(i,j+1))
                        -std::abs(v(i,  j)+v(i-1,  j))*(u(i,j-1)-u(i,  j)));
            duv *= rh/4;

            double ddux = (u(i+1,  j) - 2*u(i,j) + u(i-1,j))*rhh;
            double dduy = (u(  i,j+1) - 2*u(i,j) + u(  i,j-1))*rhh;

            fx(i,j) = u(i,j) + dt*((ddux+dduy)/re - duu -duv +gx);
        }
    }
    fx(u.firstRow(),_) = u(u.firstRow(),_);
    fx(u.lastRow(),_) = u(u.lastRow(),_);

    for (int i=v.firstRow()+1; i<=v.lastRow()-1; ++i) {
        for (int j=v.firstCol()+1; j<=v.lastCol()-1; ++j) {
            double dvv = sqr(v(i,j)   + v(i,j+1))
                       - sqr(v(i,j-1) + v(i,j))
                +gamma*( std::abs(v(i,  j)+v(i,j+1))*(v(i,  j)-v(i,j+1))
                        -std::abs(v(i,j-1)+v(i,  j))*(v(i,j-1)-v(i,  j)));
            dvv *= rh/4;

            double duv = (u(i+1,j)+u(i+1,j-1))*(v(  i,j)+v(i+1,j))
                       - (u(  i,j)+u(  i,j-1))*(v(i-1,j)+v(  i,j))
                +gamma*( std::abs(u(i+1,j)+u(i+1,j-1))*(v(  i,j)-v(i+1,j))
                        -std::abs(u(  i,j)+u(  i,j-1))*(v(i-1,j)-v(  i,j)));
            duv *= rh/4;

            double ddvx = (v(i+1,j)-2*v(i,j)+v(i-1,j))*rhh;

            double ddvy = (v(i,j+1)-2*v(i,j)+v(i,j-1))*rhh;

            fy(i,j) = v(i,j) + dt*((ddvx + ddvy)/re - duv -dvv +gy);
        }
    }
    fy(_,v.firstCol()) = v(_,v.firstCol());
    fy(_,v.lastCol()) = v(_,v.lastCol());
}

//------------------------------------------------------------------------------

template <typename TMP, typename FY>
void
computeBouyantForce(int rh, const GeMatrix<TMP> &T, double beta, double dt,
                    double gy, GeMatrix<FY> &fy)
{
    for (int i=fy.firstRow()+1; i<=fy.lastRow()-1; ++i) {
        for (int j=fy.firstCol()+1; j<=fy.lastCol()-1; ++j) {
            fy(i,j) -= beta*(T(i,j)+T(i,j-1))*0.5*gy*dt;
        }
    }
}

//------------------------------------------------------------------------------

template <typename FX, typename FY, typename RHS>
void
computePoissonRhs(int rh, const GeMatrix<FX> &fx, const GeMatrix<FY> &fy,
                  double dt, GeMatrix<RHS> &rhs)
{
    for (int i=rhs.firstRow()+1; i<=rhs.lastRow()-1; ++i) {
        for (int j=rhs.firstCol()+1; j<=rhs.lastCol()-1; ++j) {
            rhs(i,j) =-rh*((fx(i+1,j)-fx(i,j))+(fy(i,j+1)-fy(i,j)))/dt;
        }
    }
}

//------------------------------------------------------------------------------

template <typename FX, typename FY, typename P, typename U, typename V>
void
computeVelocity(int rh, const GeMatrix<FX> &fx, const GeMatrix<FY> &fy,
                const GeMatrix<P> &p, double dt,
                GeMatrix<U> &u, GeMatrix<V> &v)
{
    for (int i=u.firstRow()+1; i<=u.lastRow()-1; ++i) {
        for (int j=u.firstCol()+1; j<=u.lastCol()-1; ++j) {
            u(i,j) = fx(i,j) - rh*dt*(p(i,j)-p(i-1,j));
        }
    }
    for (int i=v.firstRow()+1; i<=v.lastRow()-1; ++i) {
        for (int j=v.firstCol()+1; j<=v.lastCol()-1; ++j) {
            v(i,j) = fy(i,j) - rh*dt*(p(i,j)-p(i,j-1));
        }
    }
}

//------------------------------------------------------------------------------

template <typename U, typename V, typename TMP>
void
computeTemperature(int rh, const GeMatrix<U> &u, const GeMatrix<V> &v,
                   double dt, double re, double pr, const GeMatrix<TMP> &T,
                   GeMatrix<TMP> &T2)
{
    for (int i = T.firstRow()+1; i<=T.lastRow()-1; ++i) {
        for (int j = T.firstCol()+1; j<=T.lastCol()-1; ++j) {
            double T_xx = (T(i-1,  j)-2*T(i,j)+T(i+1,  j))*rh*rh;
            double T_yy = (T(  i,j-1)-2*T(i,j)+T(  i,j+1))*rh*rh;

            double uc = (u(i+1,  j)+u(i,j))/2;
            double vc = (v(  i,j+1)+v(i,j))/2;

            double uT_x = (uc>0) ? uc*(T(  i,j)-T(i-1,j))*rh
                                 : uc*(T(i+1,j)-T(  i,j))*rh;

            double vT_y = (vc>0) ? vc*(T(i,  j)-T(i,j-1))*rh
                                 : vc*(T(i,j+1)-T(i,  j))*rh;

            T2(i,j) = T(i,j) + dt*((T_xx + T_yy)/(re*pr) - uT_x - vT_y);
        }
    }
}

//------------------------------------------------------------------------------

template <typename U, typename V>
void
writeVelocity(const char *filename, int counter, int rh,
              const GeMatrix<U> &u, const GeMatrix<V> &v)
{
    double h = 1./rh;
    int xStep = 1;
    int yStep = 1;

    std::ostringstream s;
    s << filename << "_"
      << std::setw(4) << std::setfill('0') << counter
      << ".dat";

    std::ofstream out(s.str().c_str());

    for (int j=v.firstCol(); j<=v.lastCol()-1; j+=yStep) {
        for (int i=u.firstRow(); i<=u.lastRow()-1; i+=xStep) {
            float vx = (u(i+1,  j) + u(i,j))*0.5;
            float vy = (v(  i,j+1) + v(i,j))*0.5;

            float nrm = sqrt(vx*vx+vy*vy);

            if (nrm>1e-15) {
                vx /= nrm;
                vy /= nrm;
            }

            out << std::setw(4) << (i+0.5)*h << " "
                << std::setw(4) << (j+0.5)*h << " "
                << std::setw(15) << vx << " "
                << std::setw(15) << vy << std::endl;
        }
    }
}

//------------------------------------------------------------------------------

template <typename TMP>
void
writeTemperature(const char *filename, int counter, int rh,
                 const GeMatrix<TMP> &T)
{
    double h = 1./rh;

    std::ostringstream s;
    s << filename << "_"
      << std::setw(4) << std::setfill('0') << counter
      << ".dat";

    std::ofstream out(s.str().c_str());

    for (int j=T.firstCol()+1; j<=T.lastCol()-1; ++j) {
        for (int i=T.firstRow()+1; i<=T.lastRow()-1; ++i) {
            out << std::setw(4) << (i+0.5)*h << " "
                << std::setw(4) << (j+0.5)*h << " "
                << std::setw(15) << T(i,j) << std::endl;
        }
        out << std::endl;
    }
}

//------------------------------------------------------------------------------

template <typename T>
ParticleField<T>::ParticleField()
    : numParticles(21*21)
{
    double h = 0.98/20;
    int count=0;

    particle = static_cast<Particle *>(calloc(numParticles, sizeof(Particle)));

    for (int i=0; i<=20; ++i) {
        for (int j=0; j<=20; ++j, ++count) {
            particle[count].x = i*h + 0.01;
            particle[count].y = j*h + 0.01;
            particle[count].u = particle[count].v = 0;
        }
    }
}

template <typename T>
ParticleField<T>::~ParticleField()
{
    free(particle);
}

template <typename T>
template <typename U, typename V>
void
ParticleField<T>::update(int rh, const GeMatrix<U> &u, const GeMatrix<V> &v, double dt)
{
    for (int p=0; p<numParticles; ++p) {
        double gridPosX = rh*particle[p].x;
        double gridPosY = rh*particle[p].y-0.5;

        int i = int(gridPosX);
        int j = int(gridPosY);

        double offsetX = gridPosX - i;
        double offsetY = gridPosY - j;

        /*
        std::cerr << "  gridPosX = " << gridPosX
                  << ", gridPosY = " << gridPosY
                  << ", i = " << i
                  << ", j = " << j
                  << ", offsetX = " << offsetX
                  << ", offsetY = " << offsetY
                  << ", u = " << u(i,j)
                  << std::endl;
        */
        particle[p].u = (1-offsetY)*((1-offsetX)*u(i,  j) + offsetX*u(i+1,  j))
                          + offsetY*((1-offsetX)*u(i,j+1) + offsetX*u(i+1,j+1));
    }

    for (int p=0; p<numParticles; ++p) {
        double gridPosX = rh*particle[p].x-0.5;
        double gridPosY = rh*particle[p].y;

        int i = int(gridPosX);
        int j = int(gridPosY);

        double offsetX = gridPosX - i;
        double offsetY = gridPosY - j;

        particle[p].v = (1-offsetX)*((1-offsetY)*v(i,  j) + offsetY*v(i,  j+1))
                          + offsetX*((1-offsetY)*v(i+1,j) + offsetY*v(i+1,j+1));
    }

    for (int p=0; p<numParticles; ++p) {
        particle[p].x += dt*particle[p].u;
        particle[p].y += dt*particle[p].v;

        if (particle[p].x<=0) {
            particle[p].x = 0.0001;
        }
        if (particle[p].y<=0) {
            particle[p].y = 0.0001;
        }
        if (particle[p].x>=1) {
            particle[p].x = 1 - 0.0001;
        }
        if (particle[p].y>=1) {
            particle[p].y = 1 - 0.0001;
        }
    }
}

template <typename T>
void
ParticleField<T>::writeToFile(const char *filename, int counter)
{
    std::ostringstream s;
    s << filename << "_"
      << std::setw(4) << std::setfill('0') << counter
      << ".dat";

    std::ofstream out(s.str().c_str());

    for (int p=0; p<numParticles; ++p) {
        out << std::setw(4) << particle[p].x << " "
            << std::setw(4) << particle[p].y << " "
            << std::setw(4) << particle[p].u << " "
            << std::setw(4) << particle[p].v << std::endl;
    }
}

} // namespace flens
