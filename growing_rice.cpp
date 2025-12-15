#include <mechsys/dem/domain.h>

using namespace DEM;

void grow(Domain &dom){
    const double dx = 0.01;
    Particle *p = dom.GetParticle(2);
    Quaternion_t q = p->Q;
    Vec3_t elongation;
    Rotation(Vec3_t(0, 0, -dx), q, elongation);
    p->Verts[1] = p->Verts[1] + elongation;
    p->Edges[0] = DEM::Edge(p->Verts[0], p->Verts[1])
}

int main(int argc, char **argv) try
{
    Domain dom;
    const double dt = 1;
    const double dtOut = 1;
    const double tf = 100;
    const double Rroot = 0.04;
    const double Rbead = 0.2;
    const double Kb_tip = 1e-3;
    const double Kb_bead = 1;
    const double rho_root = 1;

    dom.AddPlane(/*Tag*/-1, /*position, origin*/OrthoSys::O, /*R*/0.1, /*Lx*/3.0, /*Ly*/3.0, /*rho*/1.0); // by default lies on x-y plane
    dom.AddSphere(/*Tag*/1, /*position, move one radius away to avoid head on contact*/Vec3_t(0.0,Rbead,0.5), Rbead, 1.0);
    dom.AddRice(/*Tag*/2, /*position*/Vec3_t(0.0,0.0,1.0), Rroot, /*L*/0.1, rho_root, /*angle*/0.0, /*axis*/&OrthoSys::e2);

    DEM::Particle *p = dom.GetParticle(1);
    p->FixVeloc();
    p = dom.GetParticle(-1);
    p->FixVeloc();

    dom.Initialize(dt);
    dom.Solve(tf, dt, dtOut, &grow, NULL, "rice");
}
MECHSYS_CATCH