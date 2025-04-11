#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>    // for popen, pclose
#include <limits>
using namespace std;

// ------------------------------------------
// (A) Thomas Algorithm (Tridiagonal Matrix Solver)
// ------------------------------------------
vector<double> thomas(
    const vector<double>& alpha, // lower diagonal
    const vector<double>& beta,  // main diagonal
    const vector<double>& gamma, // upper diagonal
    const vector<double>& f      // right-hand side
){
    int n = (int)f.size();
    vector<double> v(n, 0.0);

    vector<double> a = alpha;
    vector<double> b = beta;
    vector<double> c = gamma;
    vector<double> r = f;

    // forward elimination
    for(int i=1; i<n; i++){
        double mult = a[i] / b[i-1];
        b[i]   -= mult * c[i-1];
        r[i]   -= mult * r[i-1];
    }
    // back substitution
    v[n-1] = r[n-1]/b[n-1];
    for(int i=n-2; i>=0; i--){
        v[i] = (r[i] - c[i]*v[i+1]) / b[i];
    }
    return v;
}

// ------------------------------------------
// (B) Plot 3D Surface to JPEG using gnuplot
//     Reverse x-axis: "set xreverse"
// ------------------------------------------
void plotSurfaceToJpeg(
    const vector<double>& x,
    const vector<double>& y,
    const vector<double>& M,  // Nx * Ny
    int Nx,
    int Ny,
    const char* outFilename
){
    FILE* gp = popen("gnuplot","w");
    if(!gp){
        cerr << "[Error] Cannot open gnuplot.\n";
        return;
    }
    // JPEG settings
    fprintf(gp, "set terminal jpeg size 800,600\n");
    fprintf(gp, "set output '%s'\n", outFilename);
    fprintf(gp, "set view 30, -132\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set zlabel 'Price'\n");

    // Reverse x-axis
    fprintf(gp, "set xreverse\n");

    // pm3d
    fprintf(gp, "set pm3d\n");

    // splot
    fprintf(gp, "splot '-' using 1:2:3 with pm3d notitle\n");

    // Nx x Ny data
    for(int j=0; j<Ny; j++){
        for(int i=0; i<Nx; i++){
            double xx = x[i];
            double yy = y[j];
            double zz = M[i + Nx*j];
            fprintf(gp, "%f %f %f\n", xx, yy, zz);
        }
        fprintf(gp, "\n"); // empty line between rows
    }
    fprintf(gp, "e\n");
    fprintf(gp, "quit\n");
    pclose(gp);

    cout << "Saved plot to " << outFilename << endl;
}

// ------------------------------------------
// (C) Main Function: 3D ELS PDE + Knock-In/Early Redemption + X-axis Reversed Chart
// ------------------------------------------
int main(){

    // ===========
    // 1) Initial Parameters
    // ===========
    double facevalue = 10000.0;
    double x_vol=0.2662, y_vol=0.2105, z_vol=0.2111;
    double rho_xy=0.279, rho_yz=0.5256, rho_xz=0.2895;
    double r=0.0165;
    double T=3.0;
    int Nt=180*3;    // 540
    double dt = T/(double)Nt;

    double x0=100.0, y0=100.0, z0=100.0;

    // Coupons/Strike Prices
    vector<double> coupon_rate={0.28,0.24,0.192,0.144,0.096,0.048};
    vector<double> strike_price={0.85,0.85,0.90,0.90,0.95,0.95};

    // Early redemption points
    vector<int> step;
    for(int i=1; i<=6; i++){
        step.push_back( (int)(Nt*(double)i/6.0) );
    }

    double dummy=0.06;  // Coupon for first period
    double kib=0.50;    // Knock-In Barrier

    // ===========
    // 2) Non-uniform Grid
    // A=[0], B=65..130 (2.5 increments), C=[160,180,200,220]
    // ===========
    vector<double> A={0.0}, B, C={160,180,200,220};
    for(double v=65.0; v<=132.4999; v+=2.5){
        B.push_back(v);
    }
    vector<double> x; 
    x.insert(x.end(), A.begin(), A.end());
    x.insert(x.end(), B.begin(), B.end());
    x.insert(x.end(), C.begin(), C.end());
    vector<double> y=x, z=x;

    int Nx=(int)x.size();
    int Ny=(int)y.size();
    int Nz=(int)z.size();

    vector<double> hx(Nx-1), hy(Ny-1), hz(Nz-1);
    for(int i=0;i<Nx-1;i++){
        hx[i] = x[i+1]- x[i];
    }
    for(int j=0;j<Ny-1;j++){
        hy[j] = y[j+1]- y[j];
    }
    for(int k=0;k<Nz-1;k++){
        hz[k] = z[k+1]- z[k];
    }

    // 3D index function
    auto idx3=[&](int i,int j,int k){
        return i + Nx*( j + Ny*k );
    };
    int Ntotal= Nx*Ny*Nz;

    // ===========
    // 3) u, ku arrays
    // ===========
    vector<double> u(Ntotal,0.0), ku(Ntotal,0.0);

    // Initial conditions
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                int idx= idx3(i,j,k);
                double xi= x[i], yj= y[j], zk= z[k];

                // Knock-in region => principal loss
                if(xi<= kib*x0 || yj<= kib*y0 || zk<= kib*z0){
                    double payoff= min({xi,yj,zk})/x0 * facevalue;
                    u[idx]  = payoff;
                    ku[idx] = payoff;
                }
                // Below first period strike => dummy coupon
                else if( xi< strike_price[0]*x0
                      || yj< strike_price[0]*y0
                      || zk< strike_price[0]*z0 ){
                    u[idx] = facevalue*(1.0+ dummy);
                    double payoff= min({xi,yj,zk})/x0 * facevalue;
                    ku[idx]= payoff;
                }
                else {
                    // First coupon
                    double val= facevalue*(1.0+ coupon_rate[0]);
                    u[idx]  = val;
                    ku[idx] = val;
                }
            }
        }
    }

    // ===========
    // 4) Tridiagonal matrix coefficients (ax, dx, cx, etc.)
    // ===========
    vector<double> ax(Nx-2), dx(Nx-2), cx(Nx-2);
    vector<double> ay(Ny-2), dy(Ny-2), cy(Ny-2);
    vector<double> az(Nz-2), dz(Nz-2), cz(Nz-2);

    // x-coefficients
    for(int i=0;i<Nx-2;i++){
        double Xmid= x[i+1];
        double hx_i= hx[i], hx_ip1= hx[i+1];
        double sig2= (x_vol*Xmid)*(x_vol*Xmid);

        ax[i] = ( -sig2 + r*Xmid*hx_ip1 )
                / ( hx_i*(hx_i+hx_ip1) );
        dx[i] = ( 1.0/dt + sig2/(hx_i*hx_ip1)
                  - r*Xmid*(hx_ip1-hx_i)/(hx_i*hx_ip1)
                  + r/3.0 );
        cx[i] = -( sig2 + r*Xmid*hx_i )
                / ( hx_ip1*(hx_i+hx_ip1) );
    }
    ax[Nx-3] -= cx[Nx-3];
    dx[Nx-3] += 2.0*cx[Nx-3];

    // y-coefficients
    for(int j=0;j<Ny-2;j++){
        double Ymid= y[j+1];
        double hy_j= hy[j], hy_jp1= hy[j+1];
        double sig2= (y_vol*Ymid)*(y_vol*Ymid);

        ay[j] = ( -sig2 + r*Ymid*hy_jp1 )
                / ( hy_j*(hy_j+hy_jp1) );
        dy[j] = ( 1.0/dt + sig2/(hy_j*hy_jp1)
                  - r*Ymid*(hy_jp1-hy_j)/(hy_j*hy_jp1)
                  + r/3.0 );
        cy[j] = -( sig2 + r*Ymid*hy_j )
                / ( hy_jp1*(hy_j+hy_jp1) );
    }
    ay[Ny-3] -= cy[Ny-3];
    dy[Ny-3] += 2.0*cy[Ny-3];

    // z-coefficients
    for(int k=0;k<Nz-2;k++){
        double Zmid= z[k+1];
        double hz_k= hz[k], hz_kp1= hz[k+1];
        double sig2= (z_vol*Zmid)*(z_vol*Zmid);

        az[k] = ( -sig2 + r*Zmid*hz_kp1 )
                / ( hz_k*(hz_k+hz_kp1) );
        dz[k] = ( 1.0/dt + sig2/(hz_k*hz_kp1)
                  - r*Zmid*(hz_kp1-hz_k)/(hz_k*hz_kp1)
                  + r/3.0 );
        cz[k] = -( sig2 + r*Zmid*hz_k )
                / ( hz_kp1*(hz_k+hz_kp1) );
    }
    az[Nz-3] -= cz[Nz-3];
    dz[Nz-3] += 2.0*cz[Nz-3];

    // Auxiliary arrays
    vector<double> old_u(u), old_ku(ku);
    vector<double> fx(Nx-2), fy(Ny-2), fz(Nz-2);

    int tag=0; // Early redemption index

    // ===========
    // 5) Time Loop (ADI)
    // ===========
    for(int iter_t=0; iter_t<Nt; iter_t++){

        // =====================
        // (a) Early Redemption Check
        // =====================
        if(tag<(int)step.size() && iter_t== step[tag]){
            double sp= strike_price[tag];
            double cr= coupon_rate[tag];
            // gx, gy, gz
            int gx=0, gy=0, gz=0;
            for(int i=0;i<Nx;i++){
                if(x[i]>= x0*sp){ gx=i; break;}
            }
            for(int j=0;j<Ny;j++){
                if(y[j]>= y0*sp){ gy=j; break;}
            }
            for(int k=0;k<Nz;k++){
                if(z[k]>= z0*sp){ gz=k; break;}
            }
            // Knock-out => facevalue*(1+ cr)
            for(int i=gx; i<Nx; i++){
                for(int j=gy; j<Ny; j++){
                    for(int k=gz; k<Nz; k++){
                        u[idx3(i,j,k)]  = facevalue*(1.0+ cr);
                        ku[idx3(i,j,k)] = facevalue*(1.0+ cr);
                    }
                }
            }
            // Knock-in => u=ku
            int kx=0, ky=0, kz=0;
            for(int i=0;i<Nx;i++){
                if(x[i]>= x0*kib){ kx=i; break;}
            }
            for(int j=0;j<Ny;j++){
                if(y[j]>= y0*kib){ ky=j; break;}
            }
            for(int k=0;k<Nz;k++){
                if(z[k]>= z0*kib){ kz=k; break;}
            }

            for(int i=0;i<=kx;i++){
                for(int j=0;j<Ny;j++){
                    for(int k=0;k<Nz;k++){
                        u[idx3(i,j,k)] = ku[idx3(i,j,k)];
                    }
                }
            }
            for(int j=0;j<=ky;j++){
                for(int i=0;i<Nx;i++){
                    for(int k=0;k<Nz;k++){
                        u[idx3(i,j,k)] = ku[idx3(i,j,k)];
                    }
                }
            }
            for(int k=0;k<=kz;k++){
                for(int i=0;i<Nx;i++){
                    for(int j=0;j<Ny;j++){
                        u[idx3(i,j,k)] = ku[idx3(i,j,k)];
                    }
                }
            }

            tag++;
        }

        // =====================
        // (b) X-Sweep (u)
        // =====================
        old_u= u;
        for(int j=1; j<Ny-1; j++){
            for(int k=1; k<Nz-1; k++){
                for(int i=1; i<(Nx-1); i++){
                    int idxFx= i-1;
                    double cxy= (1.0/3.0)*rho_xy*x_vol*y_vol* x[i]*y[j];
                    double cxz= (1.0/3.0)*rho_xz*x_vol*z_vol* x[i]*z[k];
                    double cyz= (1.0/3.0)*rho_yz*y_vol*z_vol* y[j]*z[k];

                    double txy=(
                        old_u[(i+1)+ Nx*((j+1)+ Ny*k)]
                      - old_u[(i+1)+ Nx*((j-1)+ Ny*k)]
                      - old_u[(i-1)+ Nx*((j+1)+ Ny*k)]
                      + old_u[(i-1)+ Nx*((j-1)+ Ny*k)]
                    );
                    double dxy=(
                        hx[i-1]*hy[j] + hx[i]*hy[j]
                      + hx[i]*hy[j-1] + hx[i-1]*hy[j-1]
                    );

                    double txz=(
                        old_u[(i+1)+ Nx*( j+ Ny*(k+1))]
                      - old_u[(i+1)+ Nx*( j+ Ny*(k-1))]
                      - old_u[(i-1)+ Nx*( j+ Ny*(k+1))]
                      + old_u[(i-1)+ Nx*( j+ Ny*(k-1))]
                    );
                    double dxz=(
                        hx[i-1]*hz[k] + hx[i]*hz[k]
                      + hx[i]*hz[k-1] + hx[i-1]*hz[k-1]
                    );

                    double tyz=(
                        old_u[i+ Nx*((j+1)+ Ny*(k+1))]
                      - old_u[i+ Nx*((j+1)+ Ny*(k-1))]
                      - old_u[i+ Nx*((j-1)+ Ny*(k+1))]
                      + old_u[i+ Nx*((j-1)+ Ny*(k-1))]
                    );
                    double dyz=(
                        hy[j-1]*hz[k] + hy[j]*hz[k]
                      + hy[j]*hz[k-1] + hy[j-1]*hz[k-1]
                    );

                    double crossPart= cxy*(txy/dxy)
                                    + cxz*(txz/dxz)
                                    + cyz*(tyz/dyz);

                    double mainPart= old_u[i+ Nx*( j+ Ny*k)] / dt;
                    fx[idxFx]= crossPart + mainPart;
                }
                // Thomas solver
                vector<double> sol= thomas(ax, dx, cx, fx);
                for(int i=1; i<(Nx-1); i++){
                    u[idx3(i,j,k)] = sol[i-1];
                }
            }
        }
        // Boundary correction (x-direction)
        for(int j=1;j<Ny-1;j++){
            for(int k=1;k<Nz-1;k++){
                u[(Nx-1)+ Nx*( j+ Ny*k)] = 2.0*u[(Nx-2)+ Nx*( j+ Ny*k)]
                                             - u[(Nx-3)+ Nx*( j+ Ny*k)];
            }
        }
        for(int i=1;i<(Nx-1);i++){
            for(int k=1;k<(Nz-1);k++){
                u[i+ Nx*((Ny-1)+ Ny*k)] = 2.0*u[i+ Nx*((Ny-2)+Ny*k)]
                                             - u[i+ Nx*((Ny-3)+Ny*k)];
            }
        }
        for(int i=1;i<(Nx-1);i++){
            for(int j=1;j<(Ny-1);j++){
                u[i+ Nx*( j+ Ny*(Nz-1))] = 2.0*u[i+ Nx*( j+Ny*(Nz-2))]
                                             - u[i+ Nx*( j+Ny*(Nz-3))];
            }
        }

        // =====================
        // (c) Y-Sweep (u)
        // =====================
        old_u= u;
        for(int k=1;k<Nz-1;k++){
            for(int i=1;i<Nx-1;i++){
                for(int j=1;j<Ny-1;j++){
                    int idxFy= j-1;
                    double cxy= (1.0/3.0)*rho_xy*x_vol*y_vol * x[i]*y[j];
                    double cxz= (1.0/3.0)*rho_xz*x_vol*z_vol * x[i]*z[k];
                    double cyz= (1.0/3.0)*rho_yz*y_vol*z_vol * y[j]*z[k];

                    double txy=(
                       old_u[(i+1)+ Nx*((j+1)+Ny*k)]
                     - old_u[(i+1)+ Nx*((j-1)+Ny*k)]
                     - old_u[(i-1)+ Nx*((j+1)+Ny*k)]
                     + old_u[(i-1)+ Nx*((j-1)+Ny*k)]
                    );
                    double dxy=(
                       hx[i-1]*hy[j] + hx[i]*hy[j]
                     + hx[i]*hy[j-1] + hx[i-1]*hy[j-1]
                    );

                    double txz=(
                       old_u[(i+1)+ Nx*( j+ Ny*(k+1))]
                     - old_u[(i+1)+ Nx*( j+ Ny*(k-1))]
                     - old_u[(i-1)+ Nx*( j+ Ny*(k+1))]
                     + old_u[(i-1)+ Nx*( j+ Ny*(k-1))]
                    );
                    double dxz=(
                       hx[i-1]*hz[k] + hx[i]*hz[k]
                     + hx[i]*hz[k-1] + hx[i-1]*hz[k-1]
                    );

                    double tyz=(
                       old_u[i+ Nx*((j+1)+ Ny*(k+1))]
                     - old_u[i+ Nx*((j+1)+ Ny*(k-1))]
                     - old_u[i+ Nx*((j-1)+ Ny*(k+1))]
                     + old_u[i+ Nx*((j-1)+ Ny*(k-1))]
                    );
                    double dyz=(
                       hy[j-1]*hz[k] + hy[j]*hz[k]
                     + hy[j]*hz[k-1]+ hy[j-1]*hz[k-1]
                    );

                    double crossPart= cxy*(txy/dxy)
                                    + cxz*(txz/dxz)
                                    + cyz*(tyz/dyz);

                    double mainPart= old_u[i+ Nx*( j+Ny*k)]/ dt;
                    fy[idxFy]= crossPart + mainPart;
                }
                vector<double> sol= thomas(ay, dy, cy, fy);
                for(int j=1;j<Ny-1;j++){
                    u[idx3(i,j,k)] = sol[j-1];
                }
            }
        }
        // Boundary correction (y)
        for(int i=1;i<(Nx-1);i++){
            for(int k=1;k<(Nz-1);k++){
                u[i+ Nx*((Ny-1)+ Ny*k)] = 2.0*u[i+ Nx*((Ny-2)+Ny*k)]
                                             - u[i+ Nx*((Ny-3)+Ny*k)];
            }
        }
        for(int j=1;j<(Ny-1);j++){
            for(int k=1;k<(Nz-1);k++){
                u[(Nx-1)+ Nx*( j+Ny*k)] = 2.0*u[(Nx-2)+ Nx*( j+Ny*k)]
                                             - u[(Nx-3)+ Nx*( j+Ny*k)];
            }
        }
        for(int i=1;i<(Nx-1);i++){
            for(int j=1;j<(Ny-1);j++){
                u[i+ Nx*( j+Ny*(Nz-1))] = 2.0*u[i+ Nx*( j+Ny*(Nz-2))]
                                             - u[i+ Nx*( j+Ny*(Nz-3))];
            }
        }

        // =====================
        // (d) Z-Sweep (u)
        // =====================
        old_u= u;
        for(int j=1;j<Ny-1;j++){
            for(int i=1;i<Nx-1;i++){
                for(int k=1;k<(Nz-1);k++){
                    int idxFz= k-1;
                    double cxy= (1.0/3.0)*rho_xy*x_vol*y_vol* x[i]*y[j];
                    double cxz= (1.0/3.0)*rho_xz*x_vol*z_vol* x[i]*z[k];
                    double cyz= (1.0/3.0)*rho_yz*y_vol*z_vol* y[j]*z[k];

                    double txy=(
                       old_u[(i+1)+ Nx*((j+1)+Ny*k)]
                     - old_u[(i+1)+ Nx*((j-1)+Ny*k)]
                     - old_u[(i-1)+ Nx*((j+1)+Ny*k)]
                     + old_u[(i-1)+ Nx*((j-1)+Ny*k)]
                    );
                    double dxy=(
                       hx[i-1]*hy[j]+ hx[i]*hy[j]
                     + hx[i]*hy[j-1]+ hx[i-1]*hy[j-1]
                    );

                    double txz=(
                       old_u[(i+1)+ Nx*( j+Ny*(k+1))]
                     - old_u[(i+1)+ Nx*( j+Ny*(k-1))]
                     - old_u[(i-1)+ Nx*( j+Ny*(k+1))]
                     + old_u[(i-1)+ Nx*( j+Ny*(k-1))]
                    );
                    double dxz=(
                       hx[i-1]*hz[k] + hx[i]*hz[k]
                     + hx[i]*hz[k-1]+ hx[i-1]*hz[k-1]
                    );

                    double tyz=(
                       old_u[i+ Nx*((j+1)+ Ny*(k+1))]
                     - old_u[i+ Nx*((j+1)+ Ny*(k-1))]
                     - old_u[i+ Nx*((j-1)+ Ny*(k+1))]
                     + old_u[i+ Nx*((j-1)+ Ny*(k-1))]
                    );
                    double dyz=(
                       hy[j-1]*hz[k] + hy[j]*hz[k]
                     + hy[j]*hz[k-1] + hy[j-1]*hz[k-1]
                    );

                    double crossPart= cxy*(txy/dxy)
                                    + cxz*(txz/dxz)
                                    + cyz*(tyz/dyz);

                    double mainPart= old_u[i+ Nx*( j+Ny*k)]/ dt;
                    fz[idxFz]= crossPart + mainPart;
                }
                vector<double> sol= thomas(az,dz,cz, fz);
                for(int k=1;k<(Nz-1);k++){
                    u[idx3(i,j,k)] = sol[k-1];
                }
            }
        }
        // Boundary correction (z)
        for(int i=1;i<(Nx-1);i++){
            for(int j=1;j<(Ny-1);j++){
                u[i+ Nx*( j+Ny*(Nz-1))] = 2.0*u[i+ Nx*( j+Ny*(Nz-2))]
                                             - u[i+ Nx*( j+Ny*(Nz-3))];
            }
        }
        for(int j=1;j<(Ny-1);j++){
            for(int k=1;k<(Nz-1);k++){
                u[(Nx-1)+ Nx*( j+Ny*k)] = 2.0*u[(Nx-2)+ Nx*( j+Ny*k)]
                                             - u[(Nx-3)+ Nx*( j+Ny*k)];
            }
        }
        for(int i=1;i<(Nx-1);i++){
            for(int k=1;k<(Nz-1);k++){
                u[i+ Nx*((Ny-1)+Ny*k)] = 2.0*u[i+ Nx*((Ny-2)+Ny*k)]
                                             - u[i+ Nx*((Ny-3)+Ny*k)];
            }
        }

        // =====================
        // (e) ku (Knock-in) X→Y→Z Sweeps
        // =====================
        // X-Sweep
        old_ku= ku;
        for(int j=1;j<Ny-1;j++){
            for(int k=1;k<Nz-1;k++){
                for(int i=1;i<(Nx-1);i++){
                    int idxFx= i-1;
                    double cxy= (1.0/3.0)*rho_xy*x_vol*y_vol* x[i]*y[j];
                    double cxz= (1.0/3.0)*rho_xz*x_vol*z_vol* x[i]*z[k];
                    double cyz= (1.0/3.0)*rho_yz*y_vol*z_vol* y[j]*z[k];

                    double txy=(
                       old_ku[(i+1)+ Nx*((j+1)+Ny*k)]
                     - old_ku[(i+1)+ Nx*((j-1)+Ny*k)]
                     - old_ku[(i-1)+ Nx*((j+1)+Ny*k)]
                     + old_ku[(i-1)+ Nx*((j-1)+Ny*k)]
                    );
                    double dxy=(
                       hx[i-1]*hy[j]+ hx[i]*hy[j]
                     + hx[i]*hy[j-1]+ hx[i-1]*hy[j-1]
                    );

                    double txz=(
                       old_ku[(i+1)+ Nx*( j+Ny*(k+1))]
                     - old_ku[(i+1)+ Nx*( j+Ny*(k-1))]
                     - old_ku[(i-1)+ Nx*( j+Ny*(k+1))]
                     + old_ku[(i-1)+ Nx*( j+Ny*(k-1))]
                    );
                    double dxz=(
                       hx[i-1]*hz[k]+ hx[i]*hz[k]
                     + hx[i]*hz[k-1]+ hx[i-1]*hz[k-1]
                    );

                    double tyz=(
                       old_ku[i+ Nx*((j+1)+ Ny*(k+1))]
                     - old_ku[i+ Nx*((j+1)+ Ny*(k-1))]
                     - old_ku[i+ Nx*((j-1)+ Ny*(k+1))]
                     + old_ku[i+ Nx*((j-1)+ Ny*(k-1))]
                    );
                    double dyz=(
                       hy[j-1]*hz[k]+ hy[j]*hz[k]
                     + hy[j]*hz[k-1]+ hy[j-1]*hz[k-1]
                    );

                    double crossPart= cxy*(txy/dxy)
                                    + cxz*(txz/dxz)
                                    + cyz*(tyz/dyz);
                    double mainPart= old_ku[i+ Nx*( j+Ny*k)]/dt;
                    fx[idxFx]= crossPart + mainPart;
                }
                vector<double> sol= thomas(ax,dx,cx, fx);
                for(int i=1;i<(Nx-1);i++){
                    ku[idx3(i,j,k)] = sol[i-1];
                }
            }
        }
        // Boundary correction (ku-x)
        for(int j=1;j<Ny-1;j++){
            for(int k=1;k<Nz-1;k++){
                ku[(Nx-1)+ Nx*( j+Ny*k)] =
                    2.0*ku[(Nx-2)+ Nx*( j+Ny*k)]
                  -   ku[(Nx-3)+ Nx*( j+Ny*k)];
            }
        }
        for(int i=1;i<(Nx-1);i++){
            for(int k=1;k<(Nz-1);k++){
                ku[i+ Nx*((Ny-1)+Ny*k)] =
                    2.0*ku[i+ Nx*((Ny-2)+Ny*k)]
                  -   ku[i+ Nx*((Ny-3)+Ny*k)];
            }
        }
        for(int i=1;i<(Nx-1);i++){
            for(int j=1;j<(Ny-1);j++){
                ku[i+ Nx*( j+Ny*(Nz-1))] =
                    2.0*ku[i+ Nx*( j+Ny*(Nz-2))]
                  -   ku[i+ Nx*( j+Ny*(Nz-3))];
            }
        }

        // Y-Sweep (ku)
        old_ku= ku;
        for(int k=1;k<Nz-1;k++){
            for(int i=1;i<Nx-1;i++){
                for(int j=1;j<Ny-1;j++){
                    int idxFy= j-1;
                    double cxy= (1.0/3.0)*rho_xy*x_vol*y_vol* x[i]*y[j];
                    double cxz= (1.0/3.0)*rho_xz*x_vol*z_vol* x[i]*z[k];
                    double cyz= (1.0/3.0)*rho_yz*y_vol*z_vol* y[j]*z[k];

                    double txy=(
                       old_ku[(i+1)+ Nx*((j+1)+Ny*k)]
                     - old_ku[(i+1)+ Nx*((j-1)+Ny*k)]
                     - old_ku[(i-1)+ Nx*((j+1)+Ny*k)]
                     + old_ku[(i-1)+ Nx*((j-1)+Ny*k)]
                    );
                    double dxy=(
                       hx[i-1]*hy[j]+ hx[i]*hy[j]
                     + hx[i]*hy[j-1]+ hx[i-1]*hy[j-1]
                    );

                    double txz=(
                       old_ku[(i+1)+ Nx*( j+ Ny*(k+1))]
                     - old_ku[(i+1)+ Nx*( j+ Ny*(k-1))]
                     - old_ku[(i-1)+ Nx*( j+ Ny*(k+1))]
                     + old_ku[(i-1)+ Nx*( j+ Ny*(k-1))]
                    );
                    double dxz=(
                       hx[i-1]*hz[k]+ hx[i]*hz[k]
                     + hx[i]*hz[k-1]+ hx[i-1]*hz[k-1]
                    );

                    double tyz=(
                       old_ku[i+ Nx*((j+1)+Ny*(k+1))]
                     - old_ku[i+ Nx*((j+1)+Ny*(k-1))]
                     - old_ku[i+ Nx*((j-1)+Ny*(k+1))]
                     + old_ku[i+ Nx*((j-1)+Ny*(k-1))]
                    );
                    double dyz=(
                       hy[j-1]*hz[k]+ hy[j]*hz[k]
                     + hy[j]*hz[k-1]+ hy[j-1]*hz[k-1]
                    );

                    double crossPart= cxy*(txy/dxy)
                                    + cxz*(txz/dxz)
                                    + cyz*(tyz/dyz);
                    double mainPart= old_ku[i+ Nx*( j+Ny*k)]/dt;
                    fy[idxFy]= crossPart + mainPart;
                }
                vector<double> sol= thomas(ay, dy, cy, fy);
                for(int j=1;j<Ny-1;j++){
                    ku[idx3(i,j,k)] = sol[j-1];
                }
            }
        }
        // Boundary correction (ku-y)
        for(int i=1;i<(Nx-1);i++){
            for(int k=1;k<(Nz-1);k++){
                ku[i+ Nx*((Ny-1)+Ny*k)] =
                    2.0*ku[i+ Nx*((Ny-2)+Ny*k)]
                  -   ku[i+ Nx*((Ny-3)+Ny*k)];
            }
        }
        for(int j=1;j<(Ny-1);j++){
            for(int k=1;k<(Nz-1);k++){
                ku[(Nx-1)+ Nx*( j+Ny*k)] =
                    2.0*ku[(Nx-2)+ Nx*( j+Ny*k)]
                  -   ku[(Nx-3)+ Nx*( j+Ny*k)];
            }
        }
        for(int i=1;i<(Nx-1);i++){
            for(int j=1;j<(Ny-1);j++){
                ku[i+ Nx*( j+Ny*(Nz-1))] =
                    2.0*ku[i+ Nx*( j+Ny*(Nz-2))]
                  -   ku[i+ Nx*( j+Ny*(Nz-3))];
            }
        }

        // Z-Sweep (ku)
        old_ku= ku;
        for(int j=1;j<Ny-1;j++){
            for(int i=1;i<Nx-1;i++){
                for(int k=1;k<(Nz-1);k++){
                    int idxFz= k-1;
                    double cxy= (1.0/3.0)*rho_xy*x_vol*y_vol* x[i]*y[j];
                    double cxz= (1.0/3.0)*rho_xz*x_vol*z_vol* x[i]*z[k];
                    double cyz= (1.0/3.0)*rho_yz*y_vol*z_vol* y[j]*z[k];

                    double txy=(
                       old_ku[(i+1)+ Nx*((j+1)+Ny*k)]
                     - old_ku[(i+1)+ Nx*((j-1)+Ny*k)]
                     - old_ku[(i-1)+ Nx*((j+1)+Ny*k)]
                     + old_ku[(i-1)+ Nx*((j-1)+Ny*k)]
                    );
                    double dxy=(
                       hx[i-1]*hy[j]+ hx[i]*hy[j]
                     + hx[i]*hy[j-1]+ hx[i-1]*hy[j-1]
                    );

                    double txz=(
                       old_ku[(i+1)+ Nx*( j+Ny*(k+1))]
                     - old_ku[(i+1)+ Nx*( j+Ny*(k-1))]
                     - old_ku[(i-1)+ Nx*( j+Ny*(k+1))]
                     + old_ku[(i-1)+ Nx*( j+Ny*(k-1))]
                    );
                    double dxz=(
                       hx[i-1]*hz[k]+ hx[i]*hz[k]
                     + hx[i]*hz[k-1]+ hx[i-1]*hz[k-1]
                    );

                    double tyz=(
                       old_ku[i+ Nx*((j+1)+ Ny*(k+1))]
                     - old_ku[i+ Nx*((j+1)+ Ny*(k-1))]
                     - old_ku[i+ Nx*((j-1)+ Ny*(k+1))]
                     + old_ku[i+ Nx*((j-1)+ Ny*(k-1))]
                    );
                    double dyz=(
                       hy[j-1]*hz[k]+ hy[j]*hz[k]
                     + hy[j]*hz[k-1]+ hy[j-1]*hz[k-1]
                    );

                    double crossPart= cxy*(txy/dxy)
                                    + cxz*(txz/dxz)
                                    + cyz*(tyz/dyz);
                    double mainPart= old_ku[i+ Nx*( j+Ny*k)]/dt;
                    fz[idxFz]= crossPart + mainPart;
                }
                vector<double> sol= thomas(az,dz,cz, fz);
                for(int k=1;k<(Nz-1);k++){
                    ku[idx3(i,j,k)] = sol[k-1];
                }
            }
        }
        // Boundary correction (ku-z)
        for(int i=1;i<(Nx-1);i++){
            for(int j=1;j<(Ny-1);j++){
                ku[i+ Nx*( j+Ny*(Nz-1))] =
                    2.0*ku[i+ Nx*( j+Ny*(Nz-2))]
                  -   ku[i+ Nx*( j+Ny*(Nz-3))];
            }
        }
        for(int j=1;j<(Ny-1);j++){
            for(int k=1;k<(Nz-1);k++){
                ku[(Nx-1)+ Nx*( j+Ny*k)] =
                    2.0*ku[(Nx-2)+ Nx*( j+Ny*k)]
                  -   ku[(Nx-3)+ Nx*( j+Ny*k)];
            }
        }
        for(int i=1;i<(Nx-1);i++){
            for(int k=1;k<(Nz-1);k++){
                ku[i+ Nx*((Ny-1)+Ny*k)] =
                    2.0*ku[i+ Nx*((Ny-2)+Ny*k)]
                  -   ku[i+ Nx*((Ny-3)+Ny*k)];
            }
        }

    } // end time loop

    // ===========
    // 6) Result at (x=100, y=100, z=100)
    // ===========
    int i100=-1, j100=-1, k100=-1;
    for(int i=0;i<Nx;i++){
        if(fabs(x[i]-100.0)<1e-12){ i100=i; break;}
    }
    for(int j=0;j<Ny;j++){
        if(fabs(y[j]-100.0)<1e-12){ j100=j; break;}
    }
    for(int k=0;k<Nz;k++){
        if(fabs(z[k]-100.0)<1e-12){ k100=k; break;}
    }
    if(i100<0 || j100<0 || k100<0){
        cout<<"Warning: No grid point at (x=100, y=100, z=100).\n";
    } else {
        double val= u[i100 + Nx*( j100+ Ny*k100 )];
        cout<<"Price at (100,100,100)= "<< val <<endl;
    }

    // ===========
    // 7) Plot z=100 slice to JPEG using gnuplot
    //     (x-axis reversed)
    // ===========
    if(k100>=0){
        vector<double> uSlice(Nx*Ny, 0.0);
        vector<double> kuSlice(Nx*Ny, 0.0);

        for(int j=0;j<Ny;j++){
            for(int i=0;i<Nx;i++){
                int idx2= i + Nx*j;
                uSlice[idx2]  = u[idx3(i,j,k100)];
                kuSlice[idx2] = ku[idx3(i,j,k100)];
            }
        }
        // u_slice
        plotSurfaceToJpeg(x, y, uSlice, Nx, Ny, "u_slice.jpg");
        // ku_slice
        plotSurfaceToJpeg(x, y, kuSlice, Nx, Ny, "ku_slice.jpg");
    }

    return 0;
}
