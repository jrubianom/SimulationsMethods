#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "vector.h"

struct PEFRL{
    double const1 = 0.1786178958448091;   //Zeta
    double const4 = -0.2123418310626054;  //Lambda
    double const3 = -0.06626458266981849; //Chi
    double const2 = (1-2*const4)/2;         //(1-2*Lambda)/2
    double const5 = 1-2*(const3+const1);  //1+2*(Xi+Ji)
};

struct Config{

    double L = 60.;
    double H = 120.;
    double K = 1e4;
    double E_min = 1.;
    double r0 = 10.;
    double F_mag = 12*E_min*pow(r0, -2);
    int seed = 0;

    double m = 1.;
    double r = 2.5;
    double kT = 0.1;
    double V0 = sqrt(2*kT/m);

    double tmax = 200.;
    double dt = 1e-3;
    double vis_max = 1000;
};

class Crandom{
	unsigned long long u,v,w;
	
	public:
		Crandom(unsigned long long j){
            v=4101842887655102017LL; w=1;
            u = j ^ v; int64();
            v = u; int64();
            w = v; int64();
        }
		unsigned long long int64(){
            u = u * 2862933555777941757LL + 7046029254386353087LL;
            v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
            w = 4294957665U*(w & 0xffffffff) + (w >> 32);
            unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
            return (x + v) ^ w;
        }
		double r(){
            return 5.42101086242752217E-20 * int64();
        }
		unsigned int int32(){
            return (unsigned int) int64();
        };
		double exponencial(float tau){
            return -tau*log(r());
        }
		double gauss(float mu,float sigma){
            return sigma*sqrt(-2*log(r()))*cos(2*M_PI*r())+mu;
        }
};

class Body{
    public:
        Body():
            m(1.), r(1.), R(0.,0.,0.), V(0.,0.,0.), F(0.,0.,0.)
        {};
        void Init(double m_new, double r_new,
                  double Rx, double Ry, double Rz,
                  double Vx, double Vy, double Vz){
            m = m_new;
            r = r_new;
            R.load(Rx, Ry, Rz);
            V.load(Vx, Vy, Vz);
            F.load(0., 0., 0.);
        };
        void Move_R(double dt){R +=V*dt;};
        void Move_V(double dt){V +=F*(dt/m);};
        void Add_F(vector3D Fext){F += Fext;};
        void Erase_F(){F.load(0,0,0);};
        void Print(int id, int vis_print){
            std::string fname = "data/data_" + std::to_string(vis_print) + ".txt";  
            std::ofstream fout(fname, std::ios::app);
            fout << id << "," << m << "," << r << "," 
                 << R.x() << "," << R.y() << "," << R.z() << ","
                 << V.x() << "," << V.y() << "," << V.z() << ","
                 << F.x() << "," << F.y() << "," << F.z() << "\n";
            fout.close();
        };
        friend class Collider;
    private:
        double m, r;
        vector3D R, V, F;
};

class Collider{
    public:
        Collider(Config config):
            config(config)
        {};
        void Collide(std::vector<Body> &bodies){
            //Global forces
            for (auto& b: bodies){
                //Erase previous history
                b.Erase_F();

                //Walls     
                if (double d = b.r - b.R.y(); d > 0){                   //[0,-1,0]
                    vector3D F(0., config.K*pow(d, 1.5), 0.);
                    b.Add_F(F);
                }
                else if (double d = b.r - (config.H-b.R.y()); d > 0){   //[0,1,0]
                    vector3D F(0., -config.K*pow(d, 1.5), 0.);
                    b.Add_F(F);
                }
                if (double d = b.r - b.R.x(); d > 0){                   //[-1,0,0]
                    vector3D F(config.K*pow(d, 1.5), 0., 0.);
                    b.Add_F(F);
                }
                else if (double d = b.r - (config.L-b.R.x()); d > 0){   //[1,0,0]
                    vector3D F(-config.K*pow(d, 1.5), 0., 0.);
                    b.Add_F(F);
                }
            }

            //Interaction forces
            for (int ii = 0; ii < bodies.size(); ii++){
                for (int jj = ii+1; jj < bodies.size(); jj++){
                    vector3D F(bodies[ii].R - bodies[jj].R);
                    double r = F.norm();
                    double s = r/config.r0;

                    //Fast calculation
                    F *= config.F_mag*(pow(s, -14.)-pow(s, -8.));
                    //
                    
                    //Safe calculation
                    /*if (s < 1)
                        F *= config.F_mag*pow(s, -14.)*std::abs(1 - pow(s, 6.));
                    else
                        F *= -config.F_mag*pow(s, -8.)*std::abs(1 - pow(s, -6.));*/
                    //

                    bodies[ii].Add_F(F);
                    bodies[jj].Add_F(F*(-1.));
                }
            }
        };
    private:
        Config config;
};


int main(){

    //Create global variables
    Config config;
    Crandom random(config.seed);
    Collider newton(config);
    std::vector<Body> particles(25);

    //Init system 
    for (int ii = 0; ii < 5; ii++){
        for (int jj = 0; jj < 5; jj++){
                double theta = 2*M_PI*random.r();
                particles[ii+5*jj].Init(
                config.m, config.r, 
                config.L*(ii+1)/6, config.L*(jj+1)/6, 0., 
                config.V0*cos(theta), config.V0*sin(theta), 0.);
        }
    }

    newton.Collide(particles);

    const PEFRL pefrl;
    double t = 0.;
    int iteration = 0;
    int vis_iteration = 0;
    int vis_print = 0;
    bool last;

    //Init print envioroment
    int id = 0;
    for (auto& p: particles){
        p.Print(id, vis_print);
        id++;
    }

    //Start program checks
    std::cout.precision(4);
    std::cout << "--------------------------------------------------------------\n"
              << std::left << std::setw(12)
              << "Step" << std::setw(12)
              << "Printed" << std::setw(12)
              << "Time(s)" << std::setw(12)
              << "Progress" << "\n"
              << std::left << std::setw(12)
              << "--------------------------------------------------------------\n"
              << std::left << std::setw(12)
              << 0 << std::setw(12)
              << vis_print << std::setw(12)
              << t << std::setw(12)
              << "0%" << "\n";

    for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++){

        //Check end of simulation
        last = (t >= config.tmax - 1e-8*config.dt);
        config.dt = std::min(config.dt, config.tmax - t);

        //Evolve system
        for (auto& p: particles){
            p.Move_R(config.dt*pefrl.const1);
        }
        newton.Collide(particles);
        for (auto& p: particles){
            p.Move_V(config.dt*pefrl.const2);
            p.Move_R(config.dt*pefrl.const3);
        }
        newton.Collide(particles);
        for (auto& p: particles){
            p.Move_V(config.dt*pefrl.const4);
            p.Move_R(config.dt*pefrl.const5);
        }
        newton.Collide(particles);
        for (auto& p: particles){
            p.Move_V(config.dt*pefrl.const4);
            p.Move_R(config.dt*pefrl.const3);
        }
        newton.Collide(particles);
        for (auto& p: particles){
            p.Move_V(config.dt*pefrl.const2);
            p.Move_R(config.dt*pefrl.const1);
        }

        t += config.dt;

        //Print system
        if (vis_iteration >= config.vis_max || last){
            vis_iteration = 0;
            vis_print++;

            int id = 0;
            for (auto& p: particles){
                p.Print(id, vis_print);
                id++;
            }
        }

        //Check system
        double percentage = 100*t/config.tmax;
        std::string progress = std::to_string((int)percentage)+"%";
        std::cout << std::left << std::setw(12)
                  << iteration << std::setw(12)
                  << vis_print << std::setw(12)
                  << t << std::setw(12)
                  << progress << "\n";
    }

    return 0;
}
