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
    double m = 1047.;       //Mass ratio
    double m_per = 0.005;   //Mass ratio of perturbation planet
    double r = 1000.;       //Distance between planets
    double G = 1.;          //Gravitational constant

    double w = std::sqrt(G*(1.+m)/std::pow(r,3));
    double T = (2*M_PI/w)*40;
    double dt = (2*M_PI/w)/100;

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
        double Get_Theta(){return std::atan2(R.y(), R.x());};
        void Move_R(double dt){R +=V*dt;};
        void Move_V(double dt){V +=F*(dt/m);};
        void Add_F(vector3D Fext){F += Fext;};
        void Erase_F(){F.load(0,0,0);};
        void Print(int id, int vis_print, double theta){
            std::string fname = "data/p_" + std::to_string(vis_print) + ".txt";  
            std::ofstream fout(fname, std::ios::app);
            fout << id << "," << m << "," << r << "," 
                 << std::cos(theta)*R.x() + std::sin(theta)*R.y() << "," << std::cos(theta)*R.y() - std::sin(theta)*R.x() << "," << R.z() << ","
                 << std::cos(theta)*V.x() + std::sin(theta)*V.y() << "," << std::cos(theta)*V.y() - std::sin(theta)*V.x() << "," << V.z() << ","
                 << std::cos(theta)*F.x() + std::sin(theta)*F.y() << "," << std::cos(theta)*F.y() - std::sin(theta)*F.x() << "," << F.z() << "\n";
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
            //Erase forces
            for (auto& b: bodies)
                b.Erase_F();

            //Calculate forces
            //Contact(bodies[0], bodies[1]);
            //Disturb(bodies[0], bodies[2]);
            //Disturb(bodies[1], bodies[2]);

            for (int ii = 0; ii < bodies.size(); ii++){
                for (int jj = ii+1; jj < bodies.size(); jj++){
                    Contact(bodies[ii], bodies[jj]);
                }
            }
        };
        void Contact(Body &body_1, Body &body_2){ 
            vector3D R12 = body_2.R - body_1.R;
            double d = R12.norm();
            vector3D F12 = R12*(config.G*body_1.m*body_2.m*std::pow(d, -3.));
            body_1.Add_F(F12);
            body_2.Add_F(F12*(-1.));
        };
        void Disturb(Body &body_1, Body &body_2){
            vector3D R12 = body_2.R - body_1.R;
            double d = R12.norm();
            vector3D F12 = R12*(config.G*body_1.m*body_2.m*std::pow(d, -3.));
            body_2.Add_F(F12*(-1.));
        }
    private:
        Config config;
        double G;
};


int main(){

    //Create global variables
    Config config;
    Collider newton(config);
    std::vector<Body> planets(3);

    //Init system 
    planets[0].Init(config.m, 1.);
    planets[1].Init(1., 1.);
    planets[2].Init(config.m_per, 1.);

    double r1 = config.r*1./(1.+config.m);
    double r2 = config.r*config.m/(1.+config.m);
    planets[0].Set(-r1, 0., 0., 0., -r1*config.w, 0.);
    planets[1].Set(r2, 0., 0., 0., r2*config.w, 0.);
    planets[2].Set(r2/2, r2*std::sqrt(3)/2, 0., -r2*config.w*std::sqrt(3)/2, r2*config.w/2, 0.);

    newton.Collide(planets);

    const PEFRL pefrl;
    double t = 0.;
    int iteration = 0;
    int vis_iteration = 0;
    int vis_print = 0;
    bool last;

    //Init print envioroment
    int id = 0;
    double theta = planets[1].Get_Theta();
    for (auto& p: planets){
        p.Print(id, vis_print, theta);
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
            double theta = planets[1].Get_Theta();
            for (auto& p: planets){
                p.Print(id, vis_print, theta);
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
