#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <tuple>
#include <random>

#include <boost/format.hpp>

#define NATMAX 5000
#define NPAIRMAX 40000
#define NATP 20
#define FACTE 4186.0

/// This offset is to compensate for C's indexes starting at 0 in comparison to Fortran's 1
#define OFFSET 1

struct Input
{

    Input(std::vector<std::string> &lines_input)
    {
        ///Por ahora se asignaran los valores mediante hardcoding
        file9 = "structurecg.pdb";
        file7="topcg.dat";
        file12="energy.dat";
        file20="snapshots.pdb";
        file16="beads.dat";
        file17="potentials.dat";
        tsnap=1E-12;
        tene=1E-12;
        nbloc=10000;
        temp=300;
        kkk=2857;
        fvdw=8;
        fsolv=12;
        asolv=10;
        eps=16.5;
    }

    double tsnap=1E-11;
    double tene=1E-12;
    //tterm;
    double tact=2E-14;
    double sigma=0.05;
    double temp=300.0;
    int kkk=2381;
    std::string  file7 = "topologia.dat";
    std::string  file9 = "nativain.pdb";
	std::string file10 = "res";
	std::string file11 = "distancia.dat";
	std::string file12 = "energia.dat";
    std::string file15 = "res";
    std::string file16 = "atomtypes.dat";
	std::string file17 = "potentials.dat";
    std::string file19 = "output.pdb";
	std::string file20 = "snapcg.pdb";
	std::string file21 = "snapca.pdb";
    int nbloc = 10000;
    double rcutgo=8.0;
    double tmin=1E-30;
    double dcut=10.0;
    //ego;
    int isolv=1; //Booleano
    //fpot;
    double ehb=3.0;
    double ehbc=4.0;
    int idab=0;
    //igoab
    int irig=0;
    double rshake=50.0;
    double rpot=50.0;
    //beta;
    double ebond=1000.0;
    double dstep=1E-4;
    double dijmin=1E-4;
    int isec=0; /// Boolean
    double tpush=5E-4;
    int iterm=1;
    double factm=1.0;
    //ihbr;
    double fvdw=8.0;
    double fsolv=15.0;
    double eps=16.5;
    double rbox=0.0;
    int iwr=0;
    double fterm=4.0;
    double facthc=0.8;
    double factr=0.9;
    int icons=1;
    int iprint=1;   /// Boolean
    double rsolv=3.5;
    double asolv=10.0;
    double bsolv=0.5;
    double dwat=6.0;

    int icm = 0;

    double a = 1E-10;
    double xlamb = 3.5;

};

struct Distancies
{

    Distancies(std::vector<std::string> &lines_input)
    {
        ///Por ahora se asignaran los valores mediante hardcoding

    }

    double rohmin  = 1.75;
    double rohmax  = 2.50;
    double rnomin  = 2.75;
    double rnomax  = 3.50;
    double rchmin  = 2.90;
    double rchmax  = 3.75;

    double roha = 2.15;
    double rohb = 2.34;
    double rnoa = 3.1;
	double rnob = 3.2;
	double rcha = 3.27;
	double rchb = 3.4;
};

/// Get the program parameters from input redirection ("dmdcg.dat")
auto getInput(std::istream &cin) -> std::vector<std::string>
{

    std::vector<std::string> lines;

    while (true)
    {
        std::string line;
        cin >> line;

        if (cin.eof())
            break;

        line.erase(0, 1);

        lines.push_back(line);
    }

    return lines;
}

auto getUniformRandom() -> double
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);

    static double st_rnd = 0.0;

    st_rnd += 0.05;

    if (st_rnd > 0.99)
        st_rnd = 0.05;

    return st_rnd;

    //return dist(gen);
}

struct Xoc
{
    Xoc(int natom) : r(natom), v(natom), xm(natom) { }

    std::vector<std::array<double, 3>> r;   /// Value. Size of natom x 3. Holds the coordinates of each atom
    std::vector<std::array<double, 3>> v;   /// Value. Size of natom x 3. Holds the velocity of each atom
    std::vector<double> xm; ///Value. Size of natom ....
    double rbox;
    int ierr; /// Apparently unused
};

struct Pous
{
    Pous(int natom) : rstep(natom, std::vector<std::array<double, 2>>(natom)), estep(natom, std::vector<std::array<double, 2>>(natom)) { }

    std::vector<std::vector<std::array<double, 2>>> rstep;
    std::vector<std::vector<std::array<double, 2>>> estep;
};

/// Note: even if the sizes of this class members are nres x nres, some parts of the program ask for values in the range natom x natom
struct Intr
{
    Intr(int natom) : nstep(natom, std::vector<int>(natom)), istruct(natom, std::vector<int>(natom)), inter(natom, std::vector<int>(natom)) { }

    std::vector<std::vector<int>> nstep;    /// Size. Size of nres x nres (???). Each element can have the values of 0 or 2;
                                            /// goes hand in hand with pous.rstep and pous.estep, and the value indicates if we use the pair of values of rstep/estep
                                            /// to be used to sum some energies (value=2) or not use any values to sum (value=0)
    std::vector<std::vector<int>> istruct;  /// Boolean. Size of nres x nres. Represents if a pair of atoms i,j form a hydrogen bond or not
    std::vector<std::vector<int>> inter;    /// Boolean. Size of nres x nres. Represents if a pair of atoms i,j potencials are defined (???)
};

struct Cov
{
    Cov(int natom, int npairmax) : icov(natom, std::vector<int>(natom)), rbound(npairmax), ibound(npairmax), rhc(natom) { }

    std::vector<std::vector<int>> icov;
    std::vector<double> rbound;
    std::vector<std::array<int, 2>> ibound;
    std::vector<double> rhc;
    double sigma;
};

struct Pdb
{
    Pdb(int natom) : atom(natom), res(natom), ind2(natom), nat(natom), imol(natom) { }

    std::vector<std::string> atom;  /// Atom name of each atom
    std::vector<std::string> res;   /// Residue name of each atom
    std::vector<int> ind2;  /// Index. Size of natom. Holds the residue index of each atom.
    std::vector<int> nat;   /// Number of amino groups in an atom. Each atom is specified as an atom or as a amino group, thus the number.
                            /// The values tipically aren't very large, maybe less than 10 or 15
    std::vector<int> imol;
};

struct Atpres
{
    Atpres(int natom) : ihb(natom, 0), ica(natom, 0 - OFFSET), io(natom, 0 - OFFSET), ih(natom, 0 - OFFSET), ico(natom, 0 - OFFSET), in(natom, 0 - OFFSET) { }

    std::vector<int> ihb; /// Boolean. Size of nres?. Represents whether or not the residue has a hydrogen bond
    std::vector<int> ica; /// Index. Size of nres. Each element contains the atom index of the current residue (-1 if the residue doesn't have a Ca (Calcium))
    std::vector<int> io;  /// Index. Size of nres. Each element contains the atom index of the current residue (-1 if the residue doesn't have an O (Oxygen))
    std::vector<int> ih;  /// Index. Size of nres. Each element contains the atom index of the current residue (-1 if the residue doesn't have a H (Hydrogen))
    std::vector<int> ico; /// Index. Size of nres. Each element contains the atom index of the current residue (-1 if the residue doesn't have a C (Carbon))
    std::vector<int> in;  /// Index. Size of nres. Each element contains the atom index of the current residue (-1 if the residue doesn't have a N (Nitrogen))
    //std::vector<int> icb;
};

struct Shake
{
    Shake(int natom) : ishk(natom), nshk(natom, std::vector<int>(natom, 0 - OFFSET)) { }

    std::vector<int> ishk;  /// Size. Size of natom. Represents number of possible overlaps between atoms respect to rshake2 value
    std::vector<std::vector<int>> nshk; /// Index. Size of natom x size of possible overlaps (varies depending on the atoms positions).
                                        /// Dimension 1 represents the atom, dimension 2 represents the n-overlap (if they overlap) and the element
                                        /// contains the index of the other atom that is overlapping
};

struct Npt
{
    Npt(int natom) : ipot(natom), npot(natom, std::vector<int>(natom, 0 - OFFSET)) { }

    std::vector<int> ipot;  /// Size. Size of natom. Represents number of possible overlaps between atoms respect to rpot2; the name implies some potential
    std::vector<std::vector<int>> npot; /// Index. Size of natom x size of possible overlaps (varies depending on the atoms positions).
                                        /// Dimension 1 represents the atom, dimension 2 represents the n-overlap (if they overlap) and the element
                                        /// contains the index of the other atom that is overlapping
};

/// This class members holds multiple energies and potentials of each atom
struct Fisic
{
    Fisic(int natom) : evdw(natom), rvdw(natom), qq(natom), gfree(natom), vol(natom) { }

    std::vector<double> evdw;
    std::vector<double> rvdw;
    std::vector<double> qq;
    std::vector<double> gfree;
    std::vector<double> vol;
};

/// This constants are used mainly (only) in the potencials function
struct Param
{
    double fvdw;
    double fsolv;
    double eps;
    double xlamb;
};

struct Parmsolv
{
    Parmsolv(int natom) : icont(natom), fcont(natom) { }

    double rsolv;   /// Value. Threshold of distance, with it we decide if a particle face is truncated by another particle
    double asolv;   /// Value. Its the alpha factor specified in PACSAB Appendix A
    double bsolv;   /// Value. Its the beta factor specified in PACSAB Appendix A
    double dwat;    /// Value. An offset that describes the size of a cube (from cube center to any face center); this implies that there are 6 faces,
                    /// and if we take the distance to a vertex (dwar/sqrt(3)) we have 8 vertices, adding up to 14 offsetted distances to calculate
                    /// if some particles are inside of a particle's cube
    std::vector<int> icont;     /// Size. Size of natom. Represents index of packaging; 14 distances are calculated with xoc.r (go to potencial function) and for each distance
                                /// that is less than a threshold it increments by one, so the value is never greater than 14
    std::vector<double> fcont;  /// Value. Size of natom. Its the gamma factor specified in PACSAB Appendix A; gamma≈1 -> exposed particle, gamma≈0 -> buried particle
};

struct Other
{
    Other(int natom, int natp) : qa(natom, std::vector<double>(natp)), gfreea(natom, std::vector<double>(natp)), va(natom, std::vector<double>(natp)),
                                 evdwa(natom, std::vector<double>(natp)), rvdwa(natom, std::vector<double>(natp)), rhca(natom, std::vector<double>(natp)), xma(natom, std::vector<double>(natp)),
                                 rant(natom),
                                 tpart(natom),
                                 ipart(natom), ibeta(natom), ind1(natom),
                                 inb1(natom), inb2(natom),
                                 nblist1(natom, std::vector<int>(natom)), nblist2(natom, std::vector<int>(natom)),
                                 ireg(natom, std::vector<int>(natom)),
                                 timp(natom, std::vector<double>(natom)),
                                 cad(natom),
                                 atp(natom, std::vector<std::string>(natp))
                                 { }

    std::vector<std::vector<double>> qa, gfreea, va;
    std::vector<std::vector<double>> evdwa, rvdwa, rhca, xma;
    std::vector<std::array<double, 3>> rant;
    std::vector<double> tpart;
    std::vector<int> ipart, ibeta, ind1;
    std::array<double, 3> rcm, vcm, vd, vm1, vm2;
    std::vector<int> inb1, inb2;
    std::vector<std::vector<int>> nblist1, nblist2;
    std::vector<std::vector<int>> ireg;
    std::vector<std::vector<double>> timp;
    //std::array<double. 20> v1, v2;
    std::vector<std::string> cad;
    std::vector<std::vector<std::string>> atp;  /// Value. Size of natom x pdb.nat. Holds every amino group name contained in an atom;
                                                /// the second dimension can hold elements with no string due to some atoms being only one atom (H, C, Ca, N, ...)

};

auto dbox(int n1, int n2, int k, const Xoc &xoc) -> double
{

    double rbox2 = 0.5 * xoc.rbox;
    double r12 = xoc.r[n2][k] - xoc.r[n1][k];

    if (r12 > rbox2)
        r12 -= xoc.rbox;
    else if (r12 < -rbox2)
        r12 += xoc.rbox;

    return r12;

}

auto potencial(int natom, Xoc &xoc, Pous &pous, Intr &intr, Cov &cov, Pdb &pdb, Fisic &fisic, Param &param, Parmsolv &parmsolv) -> void
{

    double rsolv2 = parmsolv.rsolv * parmsolv.rsolv;
    double rbox2 = 0.5 * xoc.rbox;
    int nres = pdb.ind2[natom - 1] + OFFSET;
    double dwatd = parmsolv.dwat / std::sqrt(3.0);

    for (int i = 0; i < natom; ++i)
    {
        std::array<std::array<double, 3>, 6> r =
        {{
            {xoc.r[i][0] + parmsolv.dwat, xoc.r[i][1]                , xoc.r[i][2]                },
            {xoc.r[i][0] - parmsolv.dwat, xoc.r[i][1]                , xoc.r[i][2]                },
            {xoc.r[i][0]                , xoc.r[i][1] + parmsolv.dwat, xoc.r[i][2]                },
            {xoc.r[i][0]                , xoc.r[i][1] - parmsolv.dwat, xoc.r[i][2]                },
            {xoc.r[i][0]                , xoc.r[i][1]                , xoc.r[i][2] + parmsolv.dwat},
            {xoc.r[i][0]                , xoc.r[i][1]                , xoc.r[i][2] - parmsolv.dwat}
        }};

        std::array<std::array<double, 3>, 8> rv =
        {{
            {xoc.r[i][0] + dwatd,   xoc.r[i][1] + dwatd,    xoc.r[i][2] + dwatd},
            {xoc.r[i][0] + dwatd,   xoc.r[i][1] + dwatd,    xoc.r[i][2] - dwatd},
            {xoc.r[i][0] + dwatd,   xoc.r[i][1] - dwatd,    xoc.r[i][2] + dwatd},
            {xoc.r[i][0] + dwatd,   xoc.r[i][1] - dwatd,    xoc.r[i][2] - dwatd},
            {xoc.r[i][0] - dwatd,   xoc.r[i][1] + dwatd,    xoc.r[i][2] + dwatd},
            {xoc.r[i][0] - dwatd,   xoc.r[i][1] + dwatd,    xoc.r[i][2] - dwatd},
            {xoc.r[i][0] - dwatd,   xoc.r[i][1] - dwatd,    xoc.r[i][2] + dwatd},
            {xoc.r[i][0] - dwatd,   xoc.r[i][1] - dwatd,    xoc.r[i][2] - dwatd}
        }};

        parmsolv.icont[i] = 0;

        /// for each element in r
        for (int r_ind = 0; r_ind < r.size(); ++r_ind)
        {
            for (int j = 0; j < natom; ++j)
            {
                if (pdb.imol[i] != pdb.imol[j])
                    continue;

                if (i == j)
                    continue;

                double rij1 = r[r_ind][0] - xoc.r[j][0];

                if (rij1 > rbox2)
                    rij1 -= xoc.rbox;
                else if (rij1 < -rbox2)
                    rij1 += xoc.rbox;

                double rij2 = r[r_ind][1] - xoc.r[j][1];

                if (rij2 > rbox2)
                    rij2 -= xoc.rbox;
                else if (rij2 < -rbox2)
                    rij2 += xoc.rbox;

                double rij3 = r[r_ind][2] - xoc.r[j][2];

                if (rij3 > rbox2)
                    rij3 -= xoc.rbox;
                else if (rij3 < -rbox2)
                    rij3 += xoc.rbox;

                double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                if (rmod2 < rsolv2)
                {
                    ++parmsolv.icont[i];
                    break;
                }
            }
        }

        /// for each element in rv (vertexs)
        for (int rv_ind = 0; rv_ind < rv.size(); ++rv_ind)
        {
            for (int j = 0; j < natom; ++j)
            {
                if (pdb.imol[i] != pdb.imol[j])
                    continue;

                if (i == j)
                    continue;

                double rij1 = rv[rv_ind][0] - xoc.r[j][0];

                if (rij1 > rbox2)
                    rij1 -= xoc.rbox;
                else if (rij1 < -rbox2)
                    rij1 += xoc.rbox;

                double rij2 = rv[rv_ind][1] - xoc.r[j][1];

                if (rij2 > rbox2)
                    rij2 -= xoc.rbox;
                else if (rij2 < -rbox2)
                    rij2 += xoc.rbox;

                double rij3 = rv[rv_ind][2] - xoc.r[j][2];

                if (rij3 > rbox2)
                    rij3 -= xoc.rbox;
                else if (rij3 < -rbox2)
                    rij3 += xoc.rbox;

                double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                if (rmod2 < rsolv2)
                {
                    ++parmsolv.icont[i];
                    break;
                }
            }
        }
    }

    for (int i = 0; i < natom; ++i)
        parmsolv.fcont[i] = 1.0 / (1.0 + std::exp((parmsolv.icont[i] - parmsolv.asolv) / parmsolv.bsolv));

    /// TODO: the is a small difference in the calculations of eij in the order of 1E-4,
    /// so we should inspect what is causing it
    for (int i = 0; i < natom - 1; ++i)
    {
        int ii = pdb.ind2[i];

        for (int j = i + 1; j < natom; ++j)
        {
            int jj = pdb.ind2[j];

            if (cov.icov[i][j] == 0)
            {
                if (intr.istruct[i][j] == 0)
                {
                    double rvdwij = fisic.rvdw[i] + fisic.rvdw[j];
                    double sto = std::pow(2.0 / (std::pow(pdb.nat[i], 0.33) + std::pow(pdb.nat[j], 0.33)), 6.0);
                    double potvdw = std::sqrt(fisic.evdw[i] * fisic.evdw[j]) * sto * (sto - 2.0);
                    double potlk = -0.09 / param.xlamb * (fisic.gfree[i] * fisic.vol[j] + fisic.gfree[j] * fisic.vol[i]) / (rvdwij * rvdwij * std::exp(std::pow(rvdwij / param.xlamb, 2.0)));
                    double eij = param.fvdw * potvdw + param.fsolv * potlk * parmsolv.fcont[i] * parmsolv.fcont[j] + param.eps * fisic.qq[i] * fisic.qq[j] / rvdwij;

                    intr.nstep[i][j] = 2 - OFFSET; //Originaly was 2, but it is an index so I decremented it
                    pous.rstep[i][j][0] = 0.9 * rvdwij;
                    pous.rstep[i][j][1] = 1.1 * rvdwij;

                    if (eij < 0.0)
                    {
                        pous.estep[i][j][0] = 3.0 * eij; //3.d0==3.0?
                        pous.estep[i][j][1] = -eij;
                    }
                    else
                    {
                        pous.estep[i][j][0] = -eij;
                        pous.estep[i][j][1] = -eij;
                    }
                }
            }
        }
    }

    return;

}

auto enchufa(int natom, double dcut, Xoc &xoc, Intr &intr, Cov &cov, Atpres &atpres, Pdb &pdb, Npt &npt) -> void
{

    double dcut2 = dcut * dcut;

    for (int i = 0; i < natom - 1; ++i)
    {
        for (int l = 0; l < npt.ipot[i]; ++l)
        {
            int j = npt.npot[i][l];

            if (pdb.nat[i] > 1 && pdb.nat[j] > 1)
            {
                double rij1 = dbox(i, j, 0, xoc);
                double rij2 = dbox(i, j, 1, xoc);
                double rij3 = dbox(i, j, 2, xoc);

                double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                if (rmod2 < dcut2)
                    intr.inter[i][j] = 1;
            }
        }
    }

    int nres = pdb.ind2[natom - 1] + OFFSET;

    for (int i = 0; i < nres - 1; ++i)
    {
        int n1 = atpres.in[i];
        int n2 = atpres.in[i + 1];

        intr.inter[n1][n2] = 0;
    }

    for (int i = 1; i < nres; ++i)
    {
        int n1 = atpres.ico[i - 1];
        int n2 = atpres.ico[i];

        intr.inter[n1][n2] = 0;
    }

    return;

}

auto creapouhb(int n1, int n2, double rmin, double r0, double r1, double rmax, double ehb, Pous &pous, Intr &intr) -> void
{

    intr.inter[n1][n2] = 1;
    intr.istruct[n1][n2] = 1;
    intr.nstep[n1][n2] = 2;

    pous.rstep[n1][n2][0] = rmin;
    pous.rstep[n1][n2][1] = rmax;

    pous.estep[n1][n2][0] = -1.5 * ehb;
    pous.estep[n1][n2][1] = 1.5 * ehb;

    return;

}

auto chgmom(int mem1, int mem2, double rij1, double rij2, double rij3, Xoc &xoc) -> void
{

    double vdmod = 0.0;

    vdmod += (xoc.v[mem2][0] - xoc.v[mem1][0]) * rij1;
    vdmod += (xoc.v[mem2][1] - xoc.v[mem1][1]) * rij2;
    vdmod += (xoc.v[mem2][2] - xoc.v[mem1][2]) * rij3;

    double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

    vdmod /= rmod2;

    double xsum = 0.5 * (1.0 / xoc.xm[mem1] + 1.0 / xoc.xm[mem2]);

    /// modul del moment transferit en la colisio
    double dp = vdmod / xsum;

    xoc.v[mem1][0] += dp / xoc.xm[mem1] * rij1;
    xoc.v[mem2][0] -= dp / xoc.xm[mem2] * rij1;

    xoc.v[mem1][1] += dp / xoc.xm[mem1] * rij2;
    xoc.v[mem2][1] -= dp / xoc.xm[mem2] * rij2;

    xoc.v[mem1][2] += dp / xoc.xm[mem1] * rij3;
    xoc.v[mem2][2] -= dp / xoc.xm[mem2] * rij3;

    return;

}

auto chgmomene(int mem1, int mem2, double rij1, double rij2, double rij3, double dpot, int &ich, Xoc &xoc) -> void
{

    double a = 1E-10;

    double vdmod = 0.0;

    vdmod += (xoc.v[mem2][0] - xoc.v[mem1][0]) * rij1;
    vdmod += (xoc.v[mem2][1] - xoc.v[mem1][1]) * rij2;
    vdmod += (xoc.v[mem2][2] - xoc.v[mem1][2]) * rij3;

    double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

    /// projeccio del moment en l'eix que uneix les dues particules
    vdmod /= rmod2;
    double xsum = 0.5 * (1.0 / xoc.xm[mem1] + 1.0 / xoc.xm[mem2]);

    /// modul del moment transferit/distancia en un xoc elastic
    double dp = vdmod / xsum;
    double sto = dp * dp / 4.0 - dpot / (rmod2 * xsum * a * a);

    if (sto > 0.0)
    {
        /// sempre es la resta dels dos valors absoluts
        if (vdmod > 0.0)
        {
            dp = dp / 2.0 - std::sqrt(sto);
        }
        else
        {
            dp = dp / 2.0 + std::sqrt(sto);
        }

        ich = 1;
    }
    else
    {
        ich = 0;
    }

    xoc.v[mem1][0] += dp / xoc.xm[mem1] * rij1;
    xoc.v[mem2][0] -= dp / xoc.xm[mem2] * rij1;

    xoc.v[mem1][1] += dp / xoc.xm[mem1] * rij2;
    xoc.v[mem2][1] -= dp / xoc.xm[mem2] * rij2;

    xoc.v[mem1][2] += dp / xoc.xm[mem1] * rij3;
    xoc.v[mem2][2] -= dp / xoc.xm[mem2] * rij3;

    return;

}

auto dmdshake(int natom, int nbound, int mem1, int mem2, Xoc &xoc, Cov &cov, Shake &shake) -> void
{

    int ierr = 0;

    /// particules NO enllaçades

    for (int i = 0; i < natom - 1; ++i)
    {
        for (int l = 0; l < shake.ishk[i]; ++l)
        {
            int j = shake.nshk[i][l];

            double rij1 = dbox(j, i, 0, xoc);
            double rij2 = dbox(j, i, 1, xoc);
            double rij3 = dbox(j, i, 2, xoc);

            double rij = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

            double vij1 = xoc.v[i][0] - xoc.v[j][0];
            double vij2 = xoc.v[i][1] - xoc.v[j][1];
            double vij3 = xoc.v[i][2] - xoc.v[j][2];

            double prod = rij1 * vij1 + rij2 * vij2 + rij3 * vij3;

            double dmin = cov.rhc[i] + cov.rhc[j];

            /// xoc frontal entre particules no enllaçades
            if (rij < dmin && prod < 0.0)
            {
                chgmom(i, j, rij1, rij2, rij3, xoc);
                ++ierr;
            }
        }
    }

    /// particules enllaçades
    for (int k = 0; k < nbound; ++k)
    {
        int i = cov.ibound[k][0];
        int j = cov.ibound[k][1];

        double rbmin = cov.rbound[k] * (1.0 - cov.sigma);
        double rbmax = cov.rbound[k] * (1.0 + cov.sigma);

        double rbmin2 = rbmin * rbmin;
        double rbmax2 = rbmax * rbmax;

        double rij1 = dbox(j, i, 0, xoc);
        double rij2 = dbox(j, i, 1, xoc);
        double rij3 = dbox(j, i, 2, xoc);

        double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

        double vij1 = xoc.v[i][0] - xoc.v[j][0];
        double vij2 = xoc.v[i][1] - xoc.v[j][1];
        double vij3 = xoc.v[i][2] - xoc.v[j][2];

        double prod = rij1 * vij1 + rij2 * vij2 + rij3 * vij3;

        if (rmod2 > rbmax2 && prod > 0.0)
        {
            ++ierr;
            chgmom(i, j, rij1, rij2, rij3, xoc);
        }
        else
        {
            if (rmod2 < rbmin2 && prod < 0.0)
            {
                ++ierr;
                chgmom(i, j, rij1, rij2, rij3, xoc);
            }
        }
    }

    return;

}

auto rnd_gauss(double &fi, double xm, double T) -> void
{
    double R = 8.314;

    double std_dev = std::sqrt((T * R) / xm);
    double pi = 3.14159265359;

    double rnd1 = getUniformRandom();
    double rnd2 = getUniformRandom();

    fi = std_dev * std::sqrt(-2.0 * std::log(rnd1)) * std::cos(2.0 * pi * rnd2);

    return;
}

auto get_molecule_info(std::string file_path, Input &input) -> std::tuple<int, Pdb, Xoc, Other>
{

    struct Atom
    {
        std::string record_type;
        int serial_number;
        std::string name;
        std::string residue_name;
        std::string chain_identifier;
        int residue_sequence_number;
        double x_orthogonal_coord;
        double y_orthogonal_coord;
        double z_orthogonal_coord;
    };

    std::ifstream file(file_path, std::ios::in);

    std::vector<Atom> atoms;

    for (int n = 0; true; ++n)
    {
        /// Each line/record in the .pdb file we're reading
        std::string atom_record;

        //file_path >> atom_record;
        std::getline(file, atom_record);

        if (file.eof())
            break;

        /// Fields of each record
        std::string record_type;
        int serial_number;
        std::string name;
        std::string residue_name;
        std::string chain_identifier;
        int residue_sequence_number;
        double x_orthogonal_coord;
        double y_orthogonal_coord;
        double z_orthogonal_coord;

        std::stringstream ss(atom_record);
        ss >> record_type >> serial_number >> name >> residue_name >> chain_identifier >> residue_sequence_number >> x_orthogonal_coord >> y_orthogonal_coord >> z_orthogonal_coord;

        Atom curr_atom;
        curr_atom.record_type = record_type;
        curr_atom.serial_number = serial_number - OFFSET;
        curr_atom.name = name;
        curr_atom.residue_name = residue_name;
        curr_atom.chain_identifier = chain_identifier;
        curr_atom.residue_sequence_number = residue_sequence_number - OFFSET;
        curr_atom.x_orthogonal_coord = x_orthogonal_coord;
        curr_atom.y_orthogonal_coord = y_orthogonal_coord;
        curr_atom.z_orthogonal_coord = z_orthogonal_coord;

        atoms.push_back(curr_atom);
    }

    int natom = atoms.size();

    Pdb pdb = Pdb(natom);
    Xoc xoc = Xoc(natom);
    xoc.rbox = input.rbox;
    Other other = Other(natom, NATP);

    for (int i = 0; i < natom; ++i)
    {
        pdb.atom[i] = atoms[i].name;
        pdb.res[i] = atoms[i].residue_name;

        other.cad[i] = atoms[i].chain_identifier;
        other.ind1[i] = atoms[i].residue_sequence_number;

        xoc.r[i] = {atoms[i].x_orthogonal_coord, atoms[i].y_orthogonal_coord, atoms[i].z_orthogonal_coord};
    }

// Correcto
//    for (int i = 0; i < natom; ++i)
//    {
//        std::cout << "atom[" << i << "]=" << pdb.atom[i] << std::endl;
//        std::cout << "res[" << i << "]=" << pdb.res[i] << std::endl;
//        std::cout << "cad[" << i << "]=" << other.cad[i] << std::endl;
//        std::cout << "ind1[" << i << "]=" << other.ind1[i] << std::endl;
//    }

    return {natom, pdb, xoc, other};

}

int main()
{
	//std::vector<std::string> dmdcg = getInput(std::cin);
	std::vector<std::string> dmdcg;
	Input input{dmdcg};
	Distancies distancies{dmdcg};

    if (input.tene > input.tsnap)
        input.tene = input.tsnap;

    if (input.tact > input.tene)
        input.tact = input.tene;

    if (input.rbox < 1E-10)
    {
        input.icm = 1;
        input.rbox = 300.0;
    }

    double rbox2 = 0.5 * input.rbox;
    double rshake2 = input.rshake * input.rshake;
    double rpot2 = input.rpot * input.rpot;

    //call random_seed()

    /// llegeix les coordenades (fitxer pdb)
    auto [natom, pdb, xoc, other] = get_molecule_info(input.file9, input);

    /// We instantiate all the structs with their vectors sized according to "natom"
    /*Xoc xoc(natom);
    Pous pous(natom);
    Intr intr(natom);
    Cov cov(natom, NPAIRMAX);
    Pdb pdb(natom);
    Atpres atpres(natom);
    Shake shake(natom);
    Npt npt(natom);
    Fisic fisic(natom);
    Param param;
    Parmsolv parmsolv(natom);*/

    int kk = 0; /// Temp variable for storing the sum of previous molecule's residue sequence numbers
    int im = 0; /// Number of molecules (A, B, C, ...)

    Atpres atpres = Atpres(natom);

    for (int n = 0; n < natom; ++n)
    {
        /// Test if we start handling another molecule
        if (n > 0 && other.ind1[n] < other.ind1[n - 1])
            kk += other.ind1[n - 1] + OFFSET;

        if (n == 0 || other.cad[n] != other.cad[n - 1])
            ++im;

        pdb.ind2[n] = other.ind1[n] + kk;

        pdb.imol[n] = im;

        int k1 = pdb.ind2[n];

        if (pdb.atom[n] == "N")  atpres.in[k1] = n;
        if (pdb.atom[n] == "H")  atpres.ih[k1] = n;
        if (pdb.atom[n] == "CA") atpres.ica[k1] = n;
        if (pdb.atom[n] == "C")  atpres.ico[k1] = n;
        if (pdb.atom[n] == "O")  atpres.io[k1] = n;
    }

    int nres = pdb.ind2[natom - 1] + OFFSET;

    /// With dynamic resize of vectors we dont need a compile time limit
    /*if (natom > NATMAX)
    {
        std::cout << "The number of particles exceeds the limit of " << NATMAX << "\n";
        return 1;
    }*/


    /// assigna un tipus a cada atom (atomtypes.dat)
    /** TODO: This section reads the file multiple times, so we can
        make it faster by reading only once into a vector and then read
        from that vector multiple times
    */
    /**
        TODO: pdb.nat[] holds the indexes of the 2nd dimension of other.atp[][],
        the more correct approach would be to ask the size of each 1st dimension
        of other.atp[][] to erase data redundancy
        Also we should push_back() the elements to each 1st dimension of other.atp[][]
        so we can save the data in much less memory space (currently the size of the matrix
        is natom x natom)
    */
    std::ifstream file16(input.file16, std::ios_base::in);

    for (int i = 0; i < natom; ++i)
    {
        pdb.nat[i] = 0;

        if (pdb.atom[i] == "N" || pdb.atom[i] == "H" || pdb.atom[i] == "C" || pdb.atom[i] == "O" || pdb.atom[i] == "OXT")
        {
            pdb.nat[i] = 1;

            /** Experimental: pass the conditionals of lines 218-223 to this section */
            if      (pdb.atom[i] == "N") other.atp[i][0] = "nh";
            else if (pdb.atom[i] == "H") other.atp[i][0] = "h";
            else if (pdb.atom[i] == "C") other.atp[i][0] = "co";
            else if (pdb.atom[i] == "O" || pdb.atom[i] == "OXT") other.atp[i][0] = "oc";

            continue;
        }

        while (true)
        {
            /// Each line/record in the .dat file we're reading
            std::string atom_type_record;

            std::getline(file16, atom_type_record);

            if (file16.eof())
                break;

            /// Fields of each record
            std::string name;
            std::string residue_name;
            std::string atp_name;

            std::stringstream ss(atom_type_record);
            ss >> name >> residue_name >> atp_name;
            ss.ignore();

            if (pdb.atom[i] == name && pdb.res[i] == residue_name)
            {
                int j = pdb.nat[i];
                other.atp[i][j] = atp_name;
                ++pdb.nat[i];
            }
        }

        /// Reset file reading to beggining
        file16.clear();
        file16.seekg(0, std::ios::beg);
    }

    file16.close();


    /// carrega els parametres de cada tipus d'atom (potentials.dat)
    /// TODO: Same as above, we are reading the file multiple times
    std::ifstream file17(input.file17, std::ios_base::in);

    /*std::vector<std::vector<double>> qa(natom);
    std::vector<std::vector<double>> gfreea(natom);
    std::vector<std::vector<double>> va(natom);
    std::vector<std::vector<double>> evdwa(natom);
    std::vector<std::vector<double>> rvdwa(natom);
    std::vector<std::vector<double>> rhca(natom);
    std::vector<std::vector<double>> xma(natom);*/

    for (int i = 0; i < natom; ++i)
    {
        for (int j = 0; j < pdb.nat[i]; ++j)
        {
            while (true)
            {
                /// Each line/record in the .dat file we're reading
                std::string atom_params_record;

                //file17 >> atom_params_record;
                std::getline(file17, atom_params_record);

                if (file17.eof())
                    break;

                /// Fields of each record
                std::string atp_name;
                double xq;
                double xfree;
                double xvol;
                double xevdw;
                double xrvdw;
                double xrhc;
                double xmassa;

                std::stringstream ss(atom_params_record);
                ss >> atp_name >> xq >> xfree >> xvol >> xevdw >> xrvdw >> xrhc >> xmassa;
                ss.ignore();

                if (other.atp[i][j] == atp_name)
                {
                    other.qa[i][j] = xq;
                    other.gfreea[i][j] = xfree;
                    other.va[i][j] = xvol;
                    other.evdwa[i][j] = xevdw;
                    other.rvdwa[i][j] = xrvdw;
                    other.rhca[i][j] = 0.8 * xrvdw; //bug? xrhc?
                    other.xma[i][j] = xmassa;

                    /// Reset file reading to beggining
                    file17.clear();
                    file17.seekg(0, std::ios::beg);

                    break;
                }
            }
        }
    }

    file17.close();


    /// propietats de les boles
    double xmassa = 0.0;

    Cov cov = Cov(natom, NPAIRMAX);
    cov.sigma = input.sigma;
    Fisic fisic = Fisic(natom);

    for (int i = 0; i < natom; ++i)
    {
        xoc.xm[i] = 0.0;
        fisic.qq[i] = 0.0;
        fisic.vol[i] = 0.0;
        fisic.gfree[i] = 0.0;
        fisic.evdw[i] = 0.0;

        double sumrhc = 0.0;
        double sumrvdw = 0.0;

        if (pdb.nat[i] == 1)
        {
            xoc.xm[i] = other.xma[i][0];

            fisic.qq[i] = other.qa[i][0];
            fisic.vol[i] = other.va[i][0];
            fisic.gfree[i] = other.gfreea[i][0];
            fisic.evdw[i] = other.evdwa[i][0];
            fisic.rvdw[i] = other.rvdwa[i][0];

            cov.rhc[i] = 0.8 * fisic.rvdw[i];
        }
        else
        {
            for (int j = 0; j < pdb.nat[i]; ++j)
            {
                xoc.xm[i] += other.xma[i][j];

                fisic.qq[i] += other.qa[i][j];
                fisic.vol[i] += other.va[i][j];
                fisic.gfree[i] += other.gfreea[i][j];
                fisic.evdw[i] += other.evdwa[i][j];

                sumrhc += std::pow(other.rhca[i][j], 3.0);
                sumrvdw += std::pow(other.rvdwa[i][j], 3.0);
            }

            fisic.rvdw[i] = input.factr * std::pow(sumrvdw, 0.333333);

            cov.rhc[i] = 0.8 * fisic.rvdw[i];
        }

        xmassa += xoc.xm[i];

        xoc.v[i] = {0.0, 0.0, 0.0};
    }


    /// llegeix la matriu de topologia (topologia.dat)

    /*std::vector<std::vector<int>> icov(natom, std::vector<int>(natom));
    std::vector<std::vector<int>> nstep(natom, std::vector<int>(natom));
    std::vector<std::vector<int>> inter(natom, std::vector<int>(natom));
    std::vector<std::vector<int>> istruct(natom, std::vector<int>(natom));

    std::vector<std::array<int, 2>> ibound(NPAIRMAX);
    std::vector<int> rbound(NPAIRMAX);*/

    std::ifstream file7(input.file7, std::ios_base::in);

    int npair = 0;

    /// The topology file is not necessary so it's not an error if it doesn't exist
    while (file7)
    {
        /// Each line/record in the .dat file we're reading
        std::string topology_record;

        std::getline(file7, topology_record);

        if (file17.eof())
            break;

        ///Fields of each record
        int i;
        int j;
        int rij;

        std::stringstream ss(topology_record);
        ss >> i >> j >> rij;

        cov.icov[i][j] = 1;
        cov.ibound[npair][0] = i - 1;
        cov.ibound[npair][1] = j - 1;
        cov.rbound[npair] = rij - 1;

        ++npair;
    }


    /** TODO: Next we write file "dmd.out" with contains all input
        parameters, for now i'm going to skip it
    */

    /// reconeix estructura secundaria i estableix ponts d'hidrogen [324]

    //std::vector<int> ihb(natom);
    int nhb = 0;

    /// ponts hidrogen O --> H
    for (int i = 0; i < nres - 4; ++i)
    {
        int ii = atpres.io[i];

        for (int j = i + 4; j < nres; ++j)
        {
            if (pdb.res[atpres.ica[j]] == "PRO")
                continue;

            int jj = atpres.ih[j];

            if (atpres.ihb[ii] != 0 || atpres.ihb[jj] != 0)
                continue;

            int n1 = atpres.io[i];
            int n2 = atpres.ih[j];

            double rij1 = dbox(n2, n1, 0, xoc);
            double rij2 = dbox(n2, n1, 1, xoc);
            double rij3 = dbox(n2, n1, 2, xoc);

            double roh = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

            if (roh > distancies.rohmax || roh < distancies.rohmin)
                continue;

            n1 = atpres.io[i];
            n2 = atpres.in[j];

            rij1 = dbox(n2, n1, 0, xoc);
            rij2 = dbox(n2, n1, 1, xoc);
            rij3 = dbox(n2, n1, 2, xoc);

            double rno = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

            if (rno > distancies.rnomax || rno < distancies.rnomin)
                continue;

            n1 = atpres.ico[i];
            n2 = atpres.ih[j];

            rij1 = dbox(n2, n1, 0, xoc);
            rij2 = dbox(n2, n1, 1, xoc);
            rij3 = dbox(n2, n1, 2, xoc);

            double rch = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

            if (rch > distancies.rchmax || rch < distancies.rchmin)
                continue;

            ++nhb;

            n1 = atpres.io[i];
            n2 = atpres.ih[j];

            atpres.ihb[n1] = 1;
            atpres.ihb[n2] = 1;

            /// write to dmd.out
            std::cout << "HBOND " << pdb.atom[ii] << " " << pdb.res[ii] << " " << pdb.ind2[ii] + OFFSET << " " << pdb.atom[jj] << " " << pdb.res[jj] << " " << pdb.ind2[jj] + OFFSET << std::endl;
        }
    }

    /// ponts hidrogen H --> O
    for (int i = 0; i < nres - 4; ++i)
    {
        if (pdb.res[atpres.ica[i]] == "PRO")
                continue;

        int ii = atpres.ih[i];

        for (int j = i + 4; j < nres; ++j)
        {
            int jj = atpres.io[j];

            if (atpres.ihb[ii] != 0 || atpres.ihb[jj] != 0)
                continue;

            int n1 = atpres.ih[i];
            int n2 = atpres.io[j];

            double rij1 = dbox(n2, n1, 0, xoc);
            double rij2 = dbox(n2, n1, 1, xoc);
            double rij3 = dbox(n2, n1, 2, xoc);

            double roh = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

            if (roh > distancies.rohmax || roh < distancies.rohmin)
                continue;

            n1 = atpres.in[i];
            n2 = atpres.io[j];

            rij1 = dbox(n2, n1, 0, xoc);
            rij2 = dbox(n2, n1, 1, xoc);
            rij3 = dbox(n2, n1, 2, xoc);

            double rno = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

            if (rno > distancies.rnomax || rno < distancies.rnomin)
                continue;

            n1 = atpres.ih[i];
            n2 = atpres.ico[j];

            rij1 = dbox(n2, n1, 0, xoc);
            rij2 = dbox(n2, n1, 1, xoc);
            rij3 = dbox(n2, n1, 2, xoc);

            double rch = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

            if (rch > distancies.rchmax || rch < distancies.rchmin)
                continue;

            ++nhb;

            n1 = atpres.ih[i];
            n2 = atpres.io[j];

            atpres.ihb[n1] = 1;
            atpres.ihb[n2] = 1;

            /// write to dmd.out
            std::cout << "HBOND " << pdb.atom[ii] << " " << pdb.res[ii] << " " << pdb.ind2[ii] + OFFSET << " " << pdb.atom[jj] << " " << pdb.res[jj] << " " << pdb.ind2[jj] + OFFSET << std::endl;
        }
    }


    Pous pous = Pous(natom);
    Intr intr = Intr(natom);

    Param param = Param(); // TODO: assign parameters from input to 'param'. A better alternative is to pass input to the constructor
    param.fvdw = input.fvdw;
    param.fsolv = input.fsolv;
    param.eps = input.eps;
    param.xlamb = input.xlamb;

    Parmsolv parmsolv = Parmsolv(natom); //TODO: same as above, assign input to 'paramsolv'. A better alternative is to pass input to the constructor
    parmsolv.rsolv = input.rsolv;
    parmsolv.asolv = input.asolv;
    parmsolv.bsolv = input.bsolv;
    parmsolv.dwat = input.dwat;

    potencial(natom, xoc, pous, intr, cov, pdb, fisic, param, parmsolv);

    /// llista de solapaments plausibles

    Shake shake = Shake(natom);
    Npt npt = Npt(natom);

    for (int i = 0; i < natom - 1; ++i)
    {
        shake.ishk[i] = 0;
        npt.ipot[i] = 0;

        for (int j = i + 1; j < natom; ++j)
        {
            if (cov.icov[i][j] == 0)
            {
                double rij1 = dbox(i, j, 0, xoc);
                double rij2 = dbox(i, j, 1, xoc);
                double rij3 = dbox(i, j, 2, xoc);

                double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                if (rmod2 < rshake2)
                {
                    int k = shake.ishk[i];  ///NOTE: different order due to C's indexes starting at 0
                    shake.nshk[i][k] = j;
                    ++shake.ishk[i];
                }

                if (intr.istruct[i][j] != 1)
                {
                    if (rmod2 < rpot2)
                    {
                        int k = npt.ipot[i];
                        npt.npot[i][k] = j;
                        ++npt.ipot[i];
                    }
                }
            }
        }
    }

    if (input.isolv != 0)
        enchufa(natom, input.dcut, xoc, intr, cov, atpres, pdb, npt);


    /// assigna la regio on es troba la interaccio entre dues particules

    for (int i = 0; i < natom - 1; ++i)
    {
        for (int j = i + 1; j < natom; ++j)
        {
            if (cov.icov[i][j] == 0)
            {
                double rij1 = dbox(j, i, 0, xoc);
                double rij2 = dbox(j, i, 1, xoc);
                double rij3 = dbox(j, i, 2, xoc);

                double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                double rij = std::sqrt(rmod2);

                int k = 1 - OFFSET;
                while (rij > pous.rstep[i][j][k] && k <= intr.nstep[i][j])
                    ++k;

                other.ireg[i][j] = k;
            }
        }
    }


    /// suma l'energia potencial de la conformacio inicial

    double epot0 = 0.0;
    double epotmol0 = 0.0;
    double epothb0 = 0.0;
    double epothbmol0 = 0.0;

    for (int i = 0; i < natom - 1; ++i)
    {//std::cout << "i=" << i << " ";
        for (int j = i + 1; j < natom; ++j)
        {//std::cout << "j=" << j << " ";
            double rij1 = dbox(i, j, 0, xoc);
            double rij2 = dbox(i, j, 1, xoc);
            double rij3 = dbox(i, j, 2, xoc);

            double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

            double dist = std::sqrt(rmod2);

            if (intr.inter[i][j] == 1)
            {
                int k = intr.nstep[i][j];

                while (dist < pous.rstep[i][j][k] && k > 0 - OFFSET)
                {
                    epot0 -= pous.estep[i][j][k];

                    if (pdb.imol[i] != pdb.imol[j])
                        epotmol0 -= pous.estep[i][j][k];

                    if (intr.istruct[i][j] == 1)
                    {
                        epothb0 -= pous.estep[i][j][k];

                        if (pdb.imol[i] != pdb.imol[j])
                            epothbmol0 -= pous.estep[i][j][k];
                    }

                    --k;
                }
            }
        }
    }


    /// assigna velocitats aleatories

    for (int j = 0; j < 3; ++j)
    {
        other.vcm[j] = 0.0;

        for (int i = 0; i < natom; ++i)
        {
            /// call random_number(fi)
            xoc.v[i][j] = getUniformRandom();
            other.vcm[j] += xoc.xm[i] * xoc.v[i][j];
        }

        other.vcm[j] /= xmassa;
    }


    /// ajusta l'energia cinetica a la temperatura requerida

    double ekin = 0.0;

    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < natom; ++i)
        {
            xoc.v[i][j] -= other.vcm[j];
            ekin += 0.5 * xoc.xm[i] * std::pow(xoc.v[i][j] * input.a, 2.0);
        }
    }

    double sto = 1.5 * 8.314 * natom * input.temp / ekin;
    double ekin0 = 0.0;

    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < natom; ++i)
        {
            xoc.v[i][j] *= std::sqrt(sto);
            ekin0 += 0.5 * xoc.xm[i] * std::pow(xoc.v[i][j] * input.a, 2.0);
        }
    }

    ekin0 /= FACTE;

    double etot0 = epot0 + ekin0;


    /// ara busca el CM

    for (int j = 0; j < 3; ++j)
    {
        other.rcm[j] = 0.0;

        for (int i = 0; i < natom; ++i)
        {
            other.rcm[j] += xoc.xm[i] * xoc.r[i][j];
        }

        other.rcm[j] /= xmassa;
    }

    /// escriu les coordenades en el SRCM

    int ibloc = 0;

    std::ofstream file_input("input.pdb", std::ios_base::out); ///TODO: rename to file12 and add it to input (input.file12)
    std::ofstream file20(input.file20, std::ios_base::out);
    std::ofstream file21(input.file21, std::ios_base::out);

    file20 << boost::format("%5s       %5d") % "MODEL" % ibloc << std::endl;
	file21 << boost::format("%5s       %5d") % "MODEL" % ibloc << std::endl;

	for (int i = 0; i < natom; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            xoc.r[i][j] += -other.rcm[j] + rbox2;

            if (xoc.r[i][j] > input.rbox)
                xoc.r[i][j] -= input.rbox;

            if (xoc.r[i][j] < 0.0)
                xoc.r[i][j] += input.rbox;
        }

        file_input << boost::format("%4s  %5d  %-3s %3s %c %3d    %8.3f%8.3f%8.3f%8.3f %2d") % "ATOM" % (i + 1) % pdb.atom[i] % pdb.res[i] % other.cad[i] % (other.ind1[i] + 1) %
                                    xoc.r[i][0] % xoc.r[i][1] % xoc.r[i][2] % parmsolv.fcont[i] % parmsolv.icont[i] << std::endl;

        file20 << boost::format("%4s  %5d  %-3s %3s %c %3d    %8.3f%8.3f%8.3f") % "ATOM" % (i + 1) % pdb.atom[i] % pdb.res[i] % other.cad[i] % (other.ind1[i] + 1) %
                                    xoc.r[i][0] % xoc.r[i][1] % xoc.r[i][2] << std::endl;

        if (pdb.atom[i] == "CA")
        {
            file21 << boost::format("%4s  %5d  %-3s %3s %c %3d    %8.3f%8.3f%8.3f") % "ATOM" % (i + 1) % pdb.atom[i] % pdb.res[i] % other.cad[i] % (other.ind1[i] + 1) %
                                    xoc.r[i][0] % xoc.r[i][1] % xoc.r[i][2] << std::endl;
        }
    }

    file_input.close();

    file20 << boost::format("%s") % "ENDMDL" << std::endl;
	file21 << boost::format("%s") % "ENDMDL" << std::endl;


    for (int i = 0; i < natom; ++i)
        for (int j = 0; j < 3; ++j)
            other.rant[i][j] = xoc.r[i][j];

    for (int i = 0; i < natom - 1; ++i)
        for (int j = i + 1; j < natom; ++j)
            other.timp[i][j] = 1.0;

    std::cout << "natom=" << natom << " nres=" << nres << " tsnap=" << input.tsnap << " nbloc=" << input.nbloc << std::endl;
    std::cout << "epot=" << epot0 << " ekin=" << ekin0 << " etot=" << etot0 << " epothb=" << epothb0 << " nhb=" << nhb << std::endl;

    std::ofstream file_output("dmd.out", std::ios_base::out); ///TODO: rename to file8 if it doesn't exist
	//write(8,INPUT)
	//write(8,DISTANCIES)
	file_output << "natom=" << natom << " nres=" << nres << " tsnap=" << input.tsnap << " nbloc=" << input.nbloc << std::endl;
	file_output << "epot=" << epot0 << " ekin=" << ekin0 << " etot=" << etot0 << " epothb=" << epothb0 << " nhb=" << nhb << std::endl;


    ///  comença a iterar per buscar l'event mes proper
    /// la variable ixoc indica quin tipus d'event ocorre:
    ///  0 -> enllaç
    ///  1 -> xoc entre atoms no enllaçats

    std::ofstream file12(input.file12, std::ios_base::out);

    double temps = 0.0;
    double temps0 = 0.0;
    double tevent = 0.0;

    file12 << "#Energia inicial " << epot0 << " " << epotmol0 << " " << epothb0 << " " << epothbmol0 << " " << ekin0 << " " << etot0 << " " << nhb << std::endl;
/// Checkpoint
    for (int ibloc = 0; ibloc < input.nbloc; ++ibloc)
    {
        double tacum = 0.0;

        double epot_i;
        double etot_i;

        while (tacum < input.tsnap)
        {
            for (int i = 0; i < natom - 1; ++i)
                for (int j = i + 1; j < natom; ++j)
                    if (intr.istruct[i][j] == 0)
                        intr.inter[i][j] = 0;

            if (input.isec == 0)
            {
                for (int i = 0; i < natom - 1; ++i)
                    for (int j = i + 1; j < natom; ++j)
                        intr.inter[i][j] = 0;

                potencial(natom, xoc, pous, intr, cov, pdb, fisic, param, parmsolv);

                int nhb = 0;

                for (int i = 0; i < natom; ++i)
                    atpres.ihb[i] = 0;

                for (int i = 0; i < nres - 4; ++i)
                {
                    int ii = atpres.io[i];

                    for (int j = i + 4; j < nres; ++j)
                    {
                        if (pdb.res[atpres.ica[j]] != "PRO")
                        {
                            int jj = atpres.ih[j];

                            /// If non of the two atoms have a hydrogen bond
                            if (atpres.ihb[ii] + atpres.ihb[jj] == 0)
                            {
                                int ic = 0;

                                int n1 = atpres.io[i];
                                int n2 = atpres.ih[j];

                                double rij1 = dbox(n2, n1, 0, xoc);
                                double rij2 = dbox(n2, n1, 1, xoc);
                                double rij3 = dbox(n2, n1, 2, xoc);

                                double roh = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                                if (roh < distancies.rohmax && roh > distancies.rohmin)
                                {
                                    ++ic;
                                }
                                else
                                {
                                    if (intr.istruct[n1][n2] == 1)
                                        if (other.ireg[n1][n2] <= intr.nstep[n1][n2] && other.ireg[n1][n2] > 1)
                                            ++ic;
                                }

                                n1 = atpres.io[i];
                                n2 = atpres.in[j];

                                rij1 = dbox(n2, n1, 0, xoc);
                                rij2 = dbox(n2, n1, 1, xoc);
                                rij3 = dbox(n2, n1, 2, xoc);

                                double rno = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                                if (rno < distancies.rnomax && rno > distancies.rnomin)
                                {
                                    ++ic;
                                }
                                else
                                {
                                    if (intr.istruct[n1][n2] == 1)
                                        if (other.ireg[n1][n2] <= intr.nstep[n1][n2] && other.ireg[n1][n2] > 1)
                                            ++ic;
                                }

                                n1 = atpres.ico[i];
                                n2 = atpres.ih[j];

                                rij1 = dbox(n2, n1, 0, xoc);
                                rij2 = dbox(n2, n1, 1, xoc);
                                rij3 = dbox(n2, n1, 2, xoc);

                                double rch = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                                if (rch < distancies.rchmax && rch > distancies.rchmin)
                                {
                                    ++ic;
                                }
                                else
                                {
                                    if (intr.istruct[n1][n2] == 1)
                                        if (other.ireg[n1][n2] <= intr.nstep[n1][n2] && other.ireg[n1][n2] > 1)
                                            ++ic;
                                }

                                if (ic == 3)
                                {
                                    ++nhb;

                                    n1 = atpres.io[i];
                                    n2 = atpres.ih[j];

                                    atpres.ihb[n1] = 1;
                                    atpres.ihb[n2] = 1;

                                    double sto = parmsolv.fcont[n1] * parmsolv.fcont[n2];
                                    double ene = input.ehb * sto + input.ehbc * (1.0 - sto);

                                    creapouhb(n1, n2, distancies.rohmin, distancies.roha, distancies.rohb, distancies.rohmax, ene, pous, intr);

                                    other.ireg[n1][n2] = 2;

                                    n1 = atpres.in[i];
                                    n2 = atpres.io[j];

                                    sto = parmsolv.fcont[n1] * parmsolv.fcont[n2];
                                    ene = input.ehb * sto + input.ehbc * (1.0 - sto);

                                    creapouhb(n1, n2, distancies.rnomin, distancies.rnoa, distancies.rnob, distancies.rnomax, ene, pous, intr);

                                    other.ireg[n1][n2] = 2;

                                    n1 = atpres.ih[i];
                                    n2 = atpres.ico[j];

                                    sto = parmsolv.fcont[n1] * parmsolv.fcont[n2];
                                    ene = input.ehb * sto + input.ehbc * (1.0 - sto);

                                    creapouhb(n1, n2, distancies.rchmin, distancies.rcha, distancies.rchb, distancies.rchmax, ene, pous, intr);

                                    other.ireg[n1][n2] = 2;
                                }
                                else
                                {
                                    n1 = atpres.ih[i];
                                    n2 = atpres.io[j];

                                    intr.istruct[n1][n2] = 0;
                                    other.ireg[n1][n2] = 0;

                                    n1 = atpres.in[i];
                                    n2 = atpres.io[j];

                                    intr.istruct[n1][n2] = 0;
                                    other.ireg[n1][n2] = 0;

                                    n1 = atpres.ih[i];
                                    n2 = atpres.ico[j];

                                    intr.istruct[n1][n2] = 0;
                                    other.ireg[n1][n2] = 0;
                                }
                            }
                        }
                    }
                }

                for (int i = 0; i < nres - 4; ++i)
                {
                    if (pdb.res[atpres.ica[i]] != "PRO")
                    {
                        int ii = atpres.ih[i];

                        for (int j = i + 4; j < nres; ++j)
                        {
                            int jj = atpres.io[j];

                            if (atpres.ihb[ii] + atpres.ihb[jj] == 0)
                            {
                                int ic = 0;

                                int n1 = atpres.ih[i];
                                int n2 = atpres.io[j];

                                double rij1 = dbox(n2, n1, 0, xoc);
                                double rij2 = dbox(n2, n1, 1, xoc);
                                double rij3 = dbox(n2, n1, 2, xoc);

                                double roh = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                                if (roh < distancies.rohmax && roh > distancies.rohmin)
                                {
                                    ++ic;
                                }
                                else
                                {
                                    if (intr.istruct[n1][n2] == 1)
                                        if (other.ireg[n1][n2] <= intr.nstep[n1][n2] && other.ireg[n1][n2] > 1)
                                            ++ic;
                                }

                                n1 = atpres.in[i];
                                n2 = atpres.io[j];

                                rij1 = dbox(n2, n1, 0, xoc);
                                rij2 = dbox(n2, n1, 1, xoc);
                                rij3 = dbox(n2, n1, 2, xoc);

                                double rno = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                                if (rno < distancies.rnomax && rno > distancies.rnomin)
                                {
                                    ++ic;
                                }
                                else
                                {
                                    if (intr.istruct[n1][n2] == 1)
                                        if (other.ireg[n1][n2] <= intr.nstep[n1][n2] && other.ireg[n1][n2] > 1)
                                            ++ic;
                                }

                                n1 = atpres.ih[i];
                                n2 = atpres.ico[j];

                                rij1 = dbox(n2, n1, 0, xoc);
                                rij2 = dbox(n2, n1, 1, xoc);
                                rij3 = dbox(n2, n1, 2, xoc);

                                double rch = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                                if (rch < distancies.rchmax && rch > distancies.rchmin)
                                {
                                    ++ic;
                                }
                                else
                                {
                                    if (intr.istruct[n1][n2] == 1)
                                        if (other.ireg[n1][n2] <= intr.nstep[n1][n2] && other.ireg[n1][n2] > 1)
                                            ++ic;
                                }

                                if (ic == 3)
                                {
                                    ++nhb;

                                    n1 = atpres.ih[i];
                                    n2 = atpres.io[j];

                                    atpres.ihb[n1] = 1;
                                    atpres.ihb[n2] = 1;

                                    double sto = parmsolv.fcont[n1] * parmsolv.fcont[n2];
                                    double ene = input.ehb * sto + input.ehbc * (1.0 - sto);

                                    creapouhb(n1, n2, distancies.rohmin, distancies.roha, distancies.rohb, distancies.rohmax, ene, pous, intr);

                                    other.ireg[n1][n2] = 2;

                                    n1 = atpres.in[i];
                                    n2 = atpres.io[j];

                                    sto = parmsolv.fcont[n1] * parmsolv.fcont[n2];
                                    ene = input.ehb * sto + input.ehbc * (1.0 - sto);

                                    creapouhb(n1, n2, distancies.rnomin, distancies.rnoa, distancies.rnob, distancies.rnomax, ene, pous, intr);

                                    other.ireg[n1][n2] = 2;

                                    n1 = atpres.ih[i];
                                    n2 = atpres.ico[j];

                                    sto = parmsolv.fcont[n1] * parmsolv.fcont[n2];
                                    ene = input.ehb * sto + input.ehbc * (1.0 - sto);

                                    creapouhb(n1, n2, distancies.rchmin, distancies.rcha, distancies.rchb, distancies.rchmax, ene, pous, intr);

                                    other.ireg[n1][n2] = 2;
                                }
                                else
                                {
                                    n1 = atpres.ih[i];
                                    n2 = atpres.io[j];

                                    intr.istruct[n1][n2] = 0;
                                    other.ireg[n1][n2] = 0;

                                    n1 = atpres.in[i];
                                    n2 = atpres.io[j];

                                    intr.istruct[n1][n2] = 0;
                                    other.ireg[n1][n2] = 0;

                                    n1 = atpres.ih[i];
                                    n2 = atpres.ico[j];

                                    intr.istruct[n1][n2] = 0;
                                    other.ireg[n1][n2] = 0;
                                }
                            }
                        }
                    }
                }
            }

            potencial(natom, xoc, pous, intr, cov, pdb, fisic, param, parmsolv);

            for (int i = 0; i < natom - 1; ++i)
            {
                for (int j = i + 1; j < natom; ++j)
                {
                    if (other.ireg[i][j] == 0 && cov.icov[i][j] == 0)
                    {
                        double rmin2 = pous.rstep[i][j][0] * pous.rstep[i][j][0];
                        double rmax2 = pous.rstep[i][j][1] * pous.rstep[i][j][1];

                        double rij1 = dbox(i, j, 0, xoc);
                        double rij2 = dbox(i, j, 1, xoc);
                        double rij3 = dbox(i, j, 2, xoc);

                        double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                        if (rmod2 < rmin2)
                        {
                            other.ireg[i][j] = 1;
                        }
                        else if (rmod2 > rmax2)
                        {
                            other.ireg[i][j] = 3;
                        }
                        else
                        {
                            other.ireg[i][j] = 2;
                        }
                    }
                }
            }

            /// llista de solapaments plausibles

            for (int i = 0; i < natom -1; ++i)
            {
                shake.ishk[i] = 0;
                npt.ipot[i] = 0;

                for (int j = i + 1; j < natom; ++j)
                {
                    if (cov.icov[i][j] == 0)
                    {
                        double rij1 = dbox(i, j, 0, xoc);
                        double rij2 = dbox(i, j, 1, xoc);
                        double rij3 = dbox(i, j, 2, xoc);

                        double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                        if (rmod2 < rshake2)
                        {
                            int k = shake.ishk[i];
                            shake.nshk[i][k] = j;
                            ++shake.ishk[i];
                        }

                        if (intr.istruct[i][j] != 1)
                        {
                            if (rmod2 < rpot2)
                            {
                                int k = npt.ipot[i];
                                npt.npot[i][k] = j;
                                ++npt.ipot[i];
                            }
                        }
                    }
                }
            }

            double mem1 = 0;
            double mem2 = 0;

            int iev = 0;
            int ierr = 0;
            int ierr2 = 0;

            temps0 = temps;

            for (int i = 0; i < natom; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    xoc.r[i][j] = other.rant[i][j];
                }
            }

            double tacene = 0.0;

            while (tacene < input.tene)
            {
                dmdshake(natom, npair, mem1, mem2, xoc, cov, shake);

                for (int i = 0; i < natom - 1; ++i)
                {
                    for (int j = i + 1; j < natom; ++j)
                    {
                        if (intr.istruct[i][j] == 0)
                            intr.inter[i][j] = 0;
                    }
                }

                if (input.isolv != 0)
                    enchufa(natom, input.dcut, xoc, intr, cov, atpres, pdb, npt);

                /// fa tornar a la seva regio els parells que han canviat de regio sense que el programa se n'adoni
                int icont = 0;

                for (int i = 0; i < natom - 1; ++i)
                {
                    for (int j = i + 1; j < natom; ++j)
                    {
                        if (intr.inter[i][j] == 1)
                        {
                            double rij1 = dbox(j, i, 0, xoc);
                            double rij2 = dbox(j, i, 1, xoc);
                            double rij3 = dbox(j, i, 2, xoc);

                            double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                            double vij1 = xoc.v[i][0] - xoc.v[j][0];
                            double vij2 = xoc.v[i][1] - xoc.v[j][1];
                            double vij3 = xoc.v[i][2] - xoc.v[j][2];

                            double vmod2 = vij1 * vij1 + vij2 * vij2 + vij3 * vij3;

                            double prod = rij1 * vij1 + rij2 * vij2 + rij3 * vij3;

                            double rij = std::sqrt(rmod2);

                            int k = 0;

                            while (rij > pous.rstep[i][j][k] && k <= intr.nstep[i][j])
                            {
                                ++k;
                            }

                            int k0 = other.ireg[i][j];

                            int ich = 0;

                            if (k0 > k && prod < 0.0)
                            {
                                ++icont;
                                double dpot = -pous.estep[i][j][k0 - 1] * FACTE;

                                chgmomene(i, j, rij1, rij2, rij3, dpot, ich, xoc);

                                if (ich == 1)
                                    --other.ireg[i][j];

                                ++ierr2;
                            }
                            else if (k0 < k && prod > 0.0)
                            {
                                ++icont;
                                double dpot = pous.estep[i][j][k0] * FACTE;

                                chgmomene(i, j, rij1, rij2, rij3, dpot, ich, xoc);

                                if (ich == 1)
                                    ++other.ireg[i][j];
                            }
                        }
                    }
                }

                /// translacio i variacio dels temps
                for (int j = 0; j < 3; ++j)
                {
                    for (int i = 0; i < natom; ++i)
                    {
                        xoc.r[i][j] += input.tact * xoc.v[i][j];
                    }
                }

                if (input.icm == 0)
                {
                    for (int i = 0; i < natom; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            if (xoc.r[i][j] > input.rbox)
                                xoc.r[i][j] -= input.rbox;
                            if (xoc.r[i][j] < 0.0)
                                xoc.r[i][j] += input.rbox;
                        }
                    }
                }

                tacum += input.tact;
                tacene += input.tact;
                //tacterm += input.tact;
                temps += input.tact;

                if (input.iterm == 1)
                {
                    /// termostat Andersen
                    /// selecciona una particula que termalitzar
                    double fi = getUniformRandom();

                    int i = int(natom * fi);
                    if (i == natom)
                        --natom;

                    for (int j = 0; j < 3; ++j)
                    {
                        rnd_gauss(fi, xoc.xm[i], input.temp);
                        xoc.v[i][j] = fi / input.a;
                    }
                }
            }

            ekin = 0.0;

            for (int j = 0; j < 3; ++j)
            {
                for (int i = 0; i < natom; ++i)
                {
                    ekin += 0.5 * xoc.xm[i] * std::pow(xoc.v[i][j] * input.a, 2.0);
                }
            }

            double ekin2 = ekin / FACTE;

            if (input.icm == 1)
            {
                for (int j = 0; j < 3; ++j)
                {
                    other.rcm[j] = 0.0;

                    for (int i = 0; i < natom; ++i)
                    {
                        other.rcm[j] += xoc.xm[i] * xoc.r[i][j];
                    }

                    other.rcm[j] /= xmassa;
                }

                for (int i = 0; i < natom; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        other.rant[i][j] = xoc.r[i][j] - other.rcm[j] + rbox2;
                    }
                }
            }
            else
            {
                for (int i = 0; i < natom; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        other.rant[i][j] = xoc.r[i][j];
                    }
                }
            }

            /// suma l'energia potencial de la conformacio inicial
            double epothb = 0.0;
            double epot = 0.0;
            static double epotmol = 0.0; ///TODO: in fortran the variable appears to mimic static behaviour the way it's declared
            static double epothbmol = 0.0; ///TODO: same here

            for (int i = 0; i < natom - 1; ++i)
            {
                for (int j = i + 1; j < natom; ++j)
                {
                    double rij1 = dbox(j, i, 0, xoc);
                    double rij2 = dbox(j, i, 1, xoc);
                    double rij3 = dbox(j, i, 2, xoc);

                    double rmod2 = rij1 * rij1 + rij2 * rij2 + rij3 * rij3;

                    double dist = std::sqrt(rmod2);

                    if (intr.inter[i][j] == 1)
                    {
                        int k = intr.nstep[i][j];

                        while (dist < pous.rstep[i][j][k] && k > 0 - OFFSET)
                        {
                            epot -= pous.estep[i][j][k];

                            if (pdb.imol[i] != pdb.imol[j])
                                epotmol -= pous.estep[i][j][k];

                            if (intr.istruct[i][j] == 1)
                            {
                                epothb -= pous.estep[i][j][k];

                                if (pdb.imol[i] != pdb.imol[j])
                                    epothbmol -= pous.estep[i][j][k];
                            }

                            --k;
                        }
                    }
                }
            }

            double etot = epot + ekin2;

            ///NEW: we added epot_i and etot_i to print epot and etot outside their current scope
            epot_i = epot;
            etot_i = etot;

            if (input.iprint == 1)
                file12 << temps << " " << epot << " " << epotmol << " " << epothb << " " << epothbmol << " " << ekin2 << " " << etot << " " << nhb << std::endl;
        }

        file20 << boost::format("%5s       %5d") % "MODEL" % ibloc << std::endl;
        file21 << boost::format("%5s       %5d") % "MODEL" % ibloc << std::endl;

        for (int i = 0; i < natom; ++i)
        {
            file20 << boost::format("%4s  %5d  %-3s %3s %c %3d    %8.3f%8.3f%8.3f") % "ATOM" % (i + 1) % pdb.atom[i] % pdb.res[i] % other.cad[i] % (other.ind1[i] + 1) %
                                    xoc.r[i][0] % xoc.r[i][1] % xoc.r[i][2] << std::endl;

            if (pdb.atom[i] == "CA")
            {
                file21 << boost::format("%4s  %5d  %-3s %3s %c %3d    %8.3f%8.3f%8.3f") % "ATOM" % (i + 1) % pdb.atom[i] % pdb.res[i] % other.cad[i] % (other.ind1[i] + 1) %
                                        xoc.r[i][0] % xoc.r[i][1] % xoc.r[i][2] << std::endl;
            }
        }

        file20 << boost::format("%s") % "ENDMDL" << std::endl;
        file21 << boost::format("%s") % "ENDMDL" << std::endl;

        std::cout << "Temps " << temps << " hbonds " << nhb << " epot " << epot_i << " " << etot_i << std::endl;
        file_output << "Temps " << temps << " hbonds " << nhb << " epot " << epot_i << " " << etot_i << std::endl;
    }

    file12.close();
    file20.close();
    file_output.close();

    ekin /= FACTE;

    std::ofstream file19(input.file19, std::ios_base::out);

    for (int i = 0; i < natom; ++i)
    {
        file19 << boost::format("%4s  %5d  %-3s %3s %c %3d    %8.3f%8.3f%8.3f") % "ATOM" % (i + 1) % pdb.atom[i] % pdb.res[i] % other.cad[i] % (other.ind1[i] + 1) %
                                    xoc.r[i][0] % xoc.r[i][1] % xoc.r[i][2] << std::endl;
    }

    file19.close();

}
