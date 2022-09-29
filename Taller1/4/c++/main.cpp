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

    double g = 980.;
    double K = 1e9;

    double m = 100.;
    double l = 12.;
    double r = 1.5;

    double tmax = 10.;
    double dt = 1e-3;
    double vis_max = 10;
};

class Body{
    public:
        Body():
            m(1.), r(1.), l(1.), theta(0.), omega(0.), tau(0.), R0(0.,1.,0.), R(0.,0.,0.)
        {};
        void Init(double m_new, double r_new, double l_new,
                  double theta_new, double omega_new, double tau_new,
                  vector3D R0_new){
            m = m_new;
            r = r_new;
            l = l_new;
            theta = theta_new;
            omega = omega_new;
            tau = tau_new;
            R0 = R0_new;
            vector3D dR(l*sin(theta), -l*cos(theta), 0.);
            R = R0 + dR;
        };
        void Move_R(double dt){
            theta += omega*dt;
            vector3D dR(l*sin(theta), -l*cos(theta), 0.);
            R = R0 + dR;
        };
        void Move_V(double dt){omega += tau*(dt/(m*pow(l, 2)));};
        void Add_F(double Fext){tau += Fext;};
        void Erase_F(){tau = 0.;};
        void Print(int id, int vis_print){
            std::string fname = "data/data_" + std::to_string(vis_print) + ".txt";  
            std::ofstream fout(fname, std::ios::app);
            fout << id << "," << m << "," << r << "," << l << ","
                 << theta << "," << omega << "," << tau << ","
                 << R.x() << "," << R.y() << "," << R.z() << ","
                 << R0.x() << "," << R0.y() << "," << R0.z() << "\n";
            fout.close();
        };
        friend class Collider;
    private:
        double m, r, l;
        double theta, omega, tau;
        vector3D R, R0;
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
                b.Add_F(-config.g*b.l*b.m*sin(b.theta));
            }

            //Interaction forces
            for (int ii = 0; ii < bodies.size(); ii++){
                for (int jj = ii+1; jj < bodies.size(); jj++){
                    vector3D F(bodies[ii].R - bodies[jj].R);
                    double d = F.norm();
                    double s = bodies[ii].r + bodies[jj].r - d;

                    if (s > 1){
                        F *= config.K*pow(s, 1.5)/d;
                        bodies[ii].Add_F(bodies[ii].l*(F.x()*cos(bodies[ii].theta)+F.y()*sin(bodies[ii].theta)));
                        bodies[jj].Add_F(-bodies[jj].l*(F.x()*cos(bodies[jj].theta)+F.y()*sin(bodies[jj].theta)));
                    }
                }
            }
        };
    private:
        Config config;
};


int main(){

    //Create global variables
    Config config;
    Collider newton(config);
    std::vector<Body> particles(3);

    //Init system 
    vector3D R0;

    R0.load(-config.r, config.l, 0.);
    particles[0].Init(config.m, config.r, config.l, -M_PI/12, 0., 0., R0);

    R0.load(0., config.l, 0.);
    particles[1].Init(config.m, config.r, config.l, 0., 0., 0., R0);

    R0.load(config.r, config.l, 0.);
    particles[2].Init(config.m, config.r, config.l, 0., 0., 0., R0);

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
