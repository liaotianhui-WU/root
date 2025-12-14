#include <mechsys/dem/domain.h>

using namespace DEM;

struct RootSegment
{
    int pid;              // particle id in DOM
    Vec3_t alpha0;        // equilibrium direction
    double Kb;            // bending stiffness to previous
    double L;             // segment length
};

std::vector<RootSegment> Root;

// parameters
const double Rroot = 0.002;     // root radius
const double Lg    = 0.01;      // max segment length
const double ug    = 1e-4;      // growth rate
const double Ke    = 1e4;       // axial stiffness
const double Kb_tip    = 1e-2;  // soft bending
const double Kb_mature = 1.0;   // stiff bending

Vec3_t Axis(const Particle &P)
{
    // rice axis is local z
    return P.Q.RotateVec(Vec3_t(0,0,1));  // Use RotateVec instead of Rotate
}

Vec3_t EndPoint(const Particle &P, int sign)
{
    return P.x + sign * 0.5 * Lg * Axis(P);
}

void UserForce(Domain &dom, void *UD)
{
    double dt = dom.Dt;

    // --- growth of the tip ---
    RootSegment &tip = Root.back();
    Particle &Ptip = dom.Particles[tip.pid];
    Ptip.L += ug * dt;

    if (Ptip.L >= Lg)
    {
        Ptip.L = Lg;

        // New segment
        Particle Pnew;
        Pnew.Shape = Shape::Rice;
        Pnew.R = Rroot;
        Pnew.L = 0.0;
        Pnew.x = EndPoint(Ptip, -1) - 1e-6 * Axis(Ptip);
        Pnew.Q = Quaternion_t(1, 0, 0, 0);  // Correct Quaternion_t constructor
        int pid = dom.Particles.Push(&Pnew);  // Pass particle pointer

        RootSegment RS;
        RS.pid = pid;
        RS.alpha0 = Axis(Ptip);
        RS.Kb = Kb_tip;
        RS.L = 0.0;
        Root.push_back(RS);

        // Remodel older segments
        for (size_t i = 0; i + 1 < Root.size(); i++)
        {
            Root[i].alpha0 = Axis(dom.Particles[Root[i].pid]);
            Root[i].Kb = Kb_mature;
        }
    }

    // --- springs between segments ---
    for (size_t i = 0; i + 1 < Root.size(); i++)
    {
        Particle &Pi = dom.Particles[Root[i].pid];
        Particle &Pj = dom.Particles[Root[i + 1].pid];

        // axial spring
        Vec3_t xi = EndPoint(Pi, -1);
        Vec3_t xj = EndPoint(Pj, +1);
        Vec3_t F = Ke * (xj - xi);
        Pi.F += F;
        Pj.F -= F;

        // bending spring
        Vec3_t ei = Axis(Pi);
        Vec3_t ej = Axis(Pj);
        Vec3_t theta = cross(ei, ej);
        Vec3_t theta0 = cross(ei, Root[i].alpha0);
        Vec3_t M = Root[i].Kb * (theta - theta0);
        Pi.T -= M;
        Pj.T += M;
    }
}

int main(int argc, char **argv)
{
    Domain dom;
    dom.Dt = 1e-5;
    dom.Gravity = Vec3_t(0,0,0);

    // --- fixed spherical obstacle ---
    Particle Obs;
    Obs.Shape = Shape::Sphere;
    Obs.R = 0.01;
    Obs.x = Vec3_t(0.0, 0.0, -0.05);
    Obs.v = Vec3_t(0,0,0);
    Obs.w = Vec3_t(0,0,0);
    Obs.FixVeloc();
    dom.Particles.Push(&Obs);

    // --- initial root seed ---
    Particle P;
    P.Shape = Shape::Rice;
    P.R = Rroot;
    P.L = 0.0;
    P.x = Vec3_t(0,0,0);
    P.Q = Quaternion_t(1,0,0,0);

    int pid = dom.Particles.Push(&P);

    RootSegment RS;
    RS.pid = pid;
    RS.alpha0 = Vec3_t(0,0,1);
    RS.Kb = Kb_tip;
    Root.push_back(RS);

    // contact parameters
    dom.SetProps(0.0, 0.0, 0.0, 0.0); // frictionless for clarity

    dom.Solve(5.0, dom.Dt, 0.0, nullptr, nullptr, "root_output.vtk");
}
