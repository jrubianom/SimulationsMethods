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
    double m = 1.;
    double r = 1.;
    double dt = 0.01;
    double T = 1.;
    double vis_max = 10;
};

class Body{
    public:
        Body():
            m(1.), r(1.)
        {
            R.load(0., 0., 0.);
            V.load(0., 0., 0.);
            F.load(0., 0., 0.);
        };
        void Init(double m_new, double r_new){
            m = m_new;
            r = r_new;
        };
        void Set(double Rx, double Ry, double Rz,
                 double Vx, double Vy, double Vz){
            R.load(Rx, Ry, Rz);
            V.load(Vx, Vy, Vz);
            F.load(0, 0, 0);
        };
        void Move_R(double dt){R +=V*dt;};
        void Move_V(double dt){V +=F*(dt/m);};
        void Add_F(vector3D Fext){F += Fext;};
        void Erase_F(){F.load(0,0,0);};
        void Print(int id, int vis_print){
            std::string fname = "data/p_" + std::to_string(vis_print) + ".dat";  
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
        Collider(double G = 1.):
            G(G)
        {};
        void Collide(std::vector<Body> &bodies){
            //Erase forces
            for (auto& b: bodies)
                b.Erase_F();

            //Calculate forces
            for (int ii = 0; ii < bodies.size(); ii++){
                for (int jj = ii+1; jj < bodies.size(); jj++){
                    Contact(bodies[ii], bodies[jj]);
                }
            }
        };
        void Contact(Body &body_1, Body &body_2){ 
            vector3D R12 = body_2.R - body_1.R;
            double d = R12.norm();
            vector3D F12 = R12*(G*body_1.m*body_2.m*std::pow(d, -3.));
            body_1.Add_F(F12);
            body_2.Add_F(F12*(-1.));
        };
    private:
        double G;
};


int main(){

    //Create global variables
    Config config;
    Collider newton;
    std::vector<Body> planets(2);

    //Init system 
    planets[0].Set(1., 0., 0., 0., 0., 0.);
    planets[1].Set(-1., 0., 0., 0., 0., 0.);
    newton.Collide(planets);

    const PEFRL pefrl;
    double t = 0.;
    int iteration = 0;
    int vis_iteration = 0;
    int vis_print = 0;
    bool last;

    //Init print envioroment
    int id = 0;
    for (auto& p: planets){
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
        last = (t >= config.T - 1e-8*config.dt);
        config.dt = std::min(config.dt, config.T - t);

        //Evolve system
        for (auto& p: planets){
            p.Move_R(config.dt*pefrl.const1);
        }
        newton.Collide(planets);
        for (auto& p: planets){
            p.Move_V(config.dt*pefrl.const2);
            p.Move_R(config.dt*pefrl.const3);
        }
        newton.Collide(planets);
        for (auto& p: planets){
            p.Move_V(config.dt*pefrl.const4);
            p.Move_R(config.dt*pefrl.const5);
        }
        newton.Collide(planets);
        for (auto& p: planets){
            p.Move_V(config.dt*pefrl.const4);
            p.Move_R(config.dt*pefrl.const3);
        }
        newton.Collide(planets);
        for (auto& p: planets){
            p.Move_V(config.dt*pefrl.const2);
            p.Move_R(config.dt*pefrl.const1);
        }

        t += config.dt;

        //Print system
        if (vis_iteration >= config.vis_max || last){
            vis_iteration = 0;
            vis_print++;

            int id = 0;
            for (auto& p: planets){
                p.Print(id, vis_print);
                id++;
            }
        }

        //Check system
        double percentage = 100*t/config.T;
        std::string progress = std::to_string((int)percentage)+"%";
        std::cout << std::left << std::setw(12)
                  << iteration << std::setw(12)
                  << vis_print << std::setw(12)
                  << t << std::setw(12)
                  << progress << "\n";
    }

    return 0;
}
