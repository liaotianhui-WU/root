/************************************************************************
 * Simple DEM script: polydisperse spheres packing by free-fall         *
 ************************************************************************/

// Std
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>

using namespace DEM;
using Util::PI;


int main(int argc, char** argv)
try
{
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    //number of threads
    size_t Nproc = 1; 
    if (argc>=3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    double verlet;  
    double rho;  
    double Rmax;         
    double Rmin;       
    double Lx;          
    double Ly;       
    double Lz;       
    double tf;
    double dt;
    double dtOut;
    double fraction;
    double Rref_mul; // Rref = Rref_mul * Rmax
    double planeR;
    double Cf;
    unsigned seed;
    double Kn;
    double Kt;
    double Gn;
    double Gt;
    double Mu;
    {
    infile >> verlet;       infile.ignore(200,'\n');
    infile >> rho;          infile.ignore(200,'\n');
    infile >> Rmax;         infile.ignore(200,'\n');
    infile >> Rmin;         infile.ignore(200,'\n');
    infile >> Lx;           infile.ignore(200,'\n');
    infile >> Ly;           infile.ignore(200,'\n');
    infile >> Lz;           infile.ignore(200,'\n');
    infile >> tf;           infile.ignore(200,'\n');
    infile >> dt;           infile.ignore(200,'\n');
    infile >> dtOut;        infile.ignore(200,'\n');
    infile >> fraction;     infile.ignore(200,'\n');
    infile >> Rref_mul;     infile.ignore(200,'\n');
    infile >> planeR;       infile.ignore(200,'\n');
    infile >> Cf;           infile.ignore(200,'\n');
    infile >> seed;         infile.ignore(200,'\n');
    infile >> Kn;           infile.ignore(200,'\n');
    infile >> Kt;           infile.ignore(200,'\n');
    infile >> Gn;           infile.ignore(200,'\n');
    infile >> Gt;           infile.ignore(200,'\n');
    infile >> Mu;           infile.ignore(200,'\n');
    }

    DEM::Domain dom;
    dom.Alpha = verlet;

    // generate a lattice of particles inside a box using GenSpheresBox (tag 1)
    // X0/X1 are the full box extents 
    Vec3_t X0(-Lx/2, -Ly/2, -Lz/2);
    Vec3_t X1( Lx/2,  Ly/2,  Lz/2);
    double Rref = Rref_mul * Rmax;          // reference radius / lattice spacing
    double RminFraction = Rmin / Rref;      // minimum radius as fraction of Rref

    dom.GenSpheresBox(1, X0, X1, Rref, rho, "HCP", seed, fraction, RminFraction);

    // set interaction properties for tag 1 and container walls (set walls later with tags -1..-5)
    Dict prps;
    prps.Set(1, "Kn Kt Gn Gt Mu Eta Beta", Kn, Kt, Gn, Gt, Mu, 0.0, 0.0);
    dom.SetProps(prps);

    // build container walls with GenOpenBox (tags -1..-5) (centered in the function)
    dom.GenOpenBox(-1, Lx+planeR*2, Ly+planeR*2, Lz+planeR*2, planeR, Cf);
    std::cout << "Total particles: " << dom.Particles.Size() << "\n";
    for (size_t i=0;i<dom.Particles.Size();++i)
    {
        Particle* p = dom.Particles[i];
        if (!p) continue;
        if (p->Tag < 0)
        {
            std::cout << "particle idx="<<i<<" Tag="<<p->Tag<<" nfaces="<<p->Faces.Size()
                    <<" pos="<<p->x<<" Dmax="<<p->Dmax<<"\n";
        }
    }

    // set same props for walls
    for (int tag=-1;tag>=-5;tag--){
        prps.Set(tag,"Kn Kt Gn Gt Mu", Kn, Kt, Gn, Gt, Mu);
    }
    dom.SetProps(prps);

    // fix container walls
    for (int tag=-1;tag>=-5;tag--){
        dom.GetParticle(tag)->FixVeloc();
    }

    // apply gravity to particles
    for (size_t i=0;i<dom.Particles.Size();++i)
    {
        Particle* p = dom.Particles[i];
        if (p->Tag == 1)
        {
            p->Ff = p->Props.m * Vec3_t(0.0, 0.0, -9.8); // gravity
        }
    }

    dom.Initialize(dt);
    dom.Solve(tf, dt, dtOut, /*setup*/ NULL, /*report*/ NULL, "soil_packing", false, Nproc);
    dom.Save("soil_dom");

    // --- dump wall-particle contacts
    // writes: idx p1_index p1_tag p2_index p2_tag cx cy cz Fx Fy Fz
    /*
    std::ofstream wf("wall_contacts_at_tf.txt");
    wf << "# idx p1_index p1_tag p2_index p2_tag cx cy cz Fx Fy Fz\n";
    for (size_t ic=0; ic<dom.CInteractons.Size(); ++ic)
    {
        DEM::CInteracton * ci = dom.CInteractons[ic];
        int tag1 = ci->P1->Tag;
        int tag2 = ci->P2->Tag;
        if (tag1 < 0 || tag2 < 0) // one side is a wall
        {
            double R1 = ci->P1->Props.R;
            double R2 = ci->P2->Props.R;
            Vec3_t contact_point = (R2*ci->P1->x + R1*ci->P2->x)/(R1 + R2);            Vec3_t force_total   = ci->Fnet + ci->Ftnet;          // normal + tangential
            wf << ic << " "
               << ci->P1->Index << " " << tag1 << " "
               << ci->P2->Index << " " << tag2 << " "
               << contact_point(0) << " " << contact_point(1) << " " << contact_point(2) << " "
               << force_total(0) << " " << force_total(1) << " " << force_total(2) << "\n";
        }
    }
    wf.close();

    // save particle list as ASCII (id, x,y,z, r, vx,vy,vz)
    std::ofstream of("polydisp_spheres.dat");
    of << "# id  xc   yc   zc   R   vx   vy   vz\n";
    size_t id = 0;
    for (size_t i=0;i<dom.Particles.Size();++i)
    {
        Particle* p = dom.Particles[i];
        if (p->Tag != 1) continue;
        of << id++ << " "
           << p->x(0) << " " << p->x(1) << " " << p->x(2) << " "
           << p->Props.R << " "
           << p->v(0) << " " << p->v(1) << " " << p->v(2) << "\n";
    }
    of.close();
    std::cout << "Written polydisp_spheres.dat with " << id << " particles\n";
*/
    return 0;
}
MECHSYS_CATCH
