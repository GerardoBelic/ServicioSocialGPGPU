#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <tuple>

#define NATMAX 5000
#define NPAIRMAX 40000
#define NATP 20

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
    int isolv=1;
    //fpot;
    int ehb=3;
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
    int isec=0;
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
    int iprint=1;
    double rsolv=3.5;
    int asolv=10;
    double bsolv=0.5;
    double dwat=6.0;

    int icm = 0;

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

struct Xoc
{
    Xoc(int natom) : r(natom), v(natom), xm(natom) { }

    std::vector<std::array<double, 3>> r;
    std::vector<std::array<double, 3>> v;
    std::vector<double> xm;
    double rbox;
    int ierr; /// Apparently unused
};

struct Pous
{
    Pous(int natom) : rstep(natom, std::vector<std::array<double, 2>>(natom)), estep(natom, std::vector<std::array<double, 2>>(natom)) { }

    std::vector<std::vector<std::array<double, 2>>> rstep;
    std::vector<std::vector<std::array<double, 2>>> estep;
};

struct Intr
{
    Intr(int natom) : nstep(natom, std::vector<int>(natom)), istruct(natom, std::vector<int>(natom)), inter(natom, std::vector<int>(natom)) { }

    std::vector<std::vector<int>> nstep;
    std::vector<std::vector<int>> istruct;
    std::vector<std::vector<int>> inter;
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

    std::vector<std::string> atom;
    std::vector<std::string> res;
    std::vector<int> ind2;
    std::vector<int> nat;
    std::vector<int> imol;
};

struct Atpres
{
    Atpres(int natom) : ihb(natom), ica(natom), io(natom), ih(natom), ico(natom), in(natom), icb(natom) { }

    std::vector<int> ihb;
    std::vector<int> ica;
    std::vector<int> io;
    std::vector<int> ih;
    std::vector<int> ico;
    std::vector<int> in;
    std::vector<int> icb;
};

struct Shake
{
    Shake(int natom) : ishk(natom), nshk(natom, std::vector<int>(natom)) { }

    std::vector<int> ishk;
    std::vector<std::vector<int>> nshk;
};

struct Npt
{
    Npt(int natom) : ipot(natom), npot(natom, std::vector<int>(natom)) { }

    std::vector<int> ipot;
    std::vector<std::vector<int>> npot;
};

struct Fisic
{
    Fisic(int natom) : evdw(natom), rvdw(natom), qq(natom), gfree(natom), vol(natom) { }

    std::vector<double> evdw;
    std::vector<double> rvdw;
    std::vector<double> qq;
    std::vector<double> gfree;
    std::vector<double> vol;
};

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

    double rsolv;
    double asolv;
    double bsolv;
    double dwat;
    std::vector<int> icont;
    std::vector<int> fcont;
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
    std::vector<std::vector<std::string>> atp;

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

auto potential(int natom, Xoc &xoc, Pous &pous, Intr &intr, Cov &cov, Pdb &pdb, Fisic &fisic, Param &param, Parmsolv &parmsolv) -> void
{

    double rsolv2 = parmsolv.rsolv * parmsolv.rsolv;
    double rbox2 = 0.5 * xoc.rbox;
    int nres = pdb.ind2[natom - 1];
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
                else if (rij1 < -rbox2) //bug?
                    rij1 += xoc.rbox;

                double rij2 = r[r_ind][1] - xoc.r[j][1];

                if (rij2 > rbox2)
                    rij2 -= xoc.rbox;
                else if (rij2 < -rbox2) //bug?
                    rij2 += xoc.rbox;

                double rij3 = r[r_ind][2] - xoc.r[j][2];

                if (rij3 > rbox2)
                    rij3 -= xoc.rbox;
                else if (rij3 < -rbox2) //bug?
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
                else if (rij1 < -rbox2) //bug?
                    rij1 += xoc.rbox;

                double rij2 = rv[rv_ind][1] - xoc.r[j][1];

                if (rij2 > rbox2)
                    rij2 -= xoc.rbox;
                else if (rij2 < -rbox2) //bug?
                    rij2 += xoc.rbox;

                double rij3 = rv[rv_ind][2] - xoc.r[j][2];

                if (rij3 > rbox2)
                    rij3 -= xoc.rbox;
                else if (rij3 < -rbox2) //bug?
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

                    intr.nstep[i][j] = 1; //Originaly was 2, but it is an index so I decremented it
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

auto get_molecule_info(std::string file_path) -> std::tuple<int, Pdb, Xoc, Other>
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
        curr_atom.serial_number = serial_number - 1;
        curr_atom.name = name;
        curr_atom.residue_name = residue_name;
        curr_atom.chain_identifier = chain_identifier;
        curr_atom.residue_sequence_number = residue_sequence_number - 1;
        curr_atom.x_orthogonal_coord = x_orthogonal_coord;
        curr_atom.y_orthogonal_coord = y_orthogonal_coord;
        curr_atom.z_orthogonal_coord = z_orthogonal_coord;

        atoms.push_back(curr_atom);
    }

    int natom = atoms.size();

    Pdb pdb = Pdb(natom);
    Xoc xoc = Xoc(natom);
    Other other = Other(natom, NATP);

    for (int i = 0; i < natom; ++i)
    {
        pdb.atom[i] = atoms[i].name;
        pdb.res[i] = atoms[i].residue_name;

        other.cad[i] = atoms[i].chain_identifier;
        other.ind1[i] = atoms[i].residue_sequence_number;

        xoc.r[i] = {atoms[i].x_orthogonal_coord, atoms[i].y_orthogonal_coord, atoms[i].z_orthogonal_coord};
    }

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
std::cout << "a\n";
    /// llegeix les coordenades (fitxer pdb)
    auto [natom, pdb, xoc, other] = get_molecule_info(input.file9);

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

    int kk = 0;
    int im = 0;

    for (int n = 0; n < natom; ++n)
    {
        /// Test if we start handling another molecule
        if (n > 0 && atoms[n].residue_sequence_number < atoms[n-1].residue_sequence_number)
            kk += atoms[n-1].residue_sequence_number;

        if (n > 0 && atoms[n].chain_identifier != atoms[n-1].chain_identifier)
            ++im;

        ind2[n] = atoms[n].residue_sequence_number + kk;

        imol[n] = im;

        int k1 = ind2[n];

        r[n] = {atoms[n].x_orthogonal_coord, atoms[n].y_orthogonal_coord, atoms[n].z_orthogonal_coord});

        if (atoms[n].name == "N")  in[k1] = n;
        if (atoms[n].name == "H")  ih[k1] = n;
        if (atoms[n].name == "CA") ica[k1] = n;
        if (atoms[n].name == "C")  ico[k1] = n;
        if (atoms[n].name == "O")  io[k1] = n;
    }



    int nres = ind2[natom - 1];

    if (natom > NATMAX)
    {
        std::cout << "The number of particles exceeds the limit of " << NATMAX << "\n";
        return 1;
    }

std::cout << "b\n";
    /// assigna un tipus a cada atom (atomtypes.dat)
    /** TODO: This section reads the file multiple times, so we can
        make it faster by reading only once into a vector and then read
        from that fector multiple times
    */
    std::ifstream file16(input.file16, std::ios_base::in);

    std::vector<int> nat(natom);
    std::vector<std::vector<std::string>> atp(natom);

    for (int i = 0; i < natom; ++i)
    {
        nat[i] = 0;

        if (atoms[i].name == "N" || atoms[i].name == "H" || atoms[i].name == "C" || atoms[i].name == "O" || atoms[i].name == "OXT")
        {
            /** Experimental: pass the conditionals of lines 218-223 to this section */
            nat[i] = 1;

            if      (atoms[i].name == "N") atp[i].push_back("nh");
            else if (atoms[i].name == "H") atp[i].push_back("h");
            else if (atoms[i].name == "C") atp[i].push_back("co");
            else if (atoms[i].name == "O" || atoms[i].name == "OXT") atp[i].push_back("oc");

            continue;
        }

        while (true)
        {
            /// Each line/record in the .dat file we're reading
            std::string atom_type_record;

            //file16 >> atom_type_record;
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

            if (atoms[i].name == name && atoms[i].residue_name == residue_name)
            {
                atp[i].push_back(atp_name);
                ++nat[i];
            }
        }

        /// Reset file reading to beggining
        file16.clear();
        file16.seekg(0, std::ios::beg);
    }

    file16.close();

std::cout << "c\n";
    /// carrega els parametres de cada tipus d'atom (potentials.dat)
    /// Same as above, we are reading the file multiple times
    std::ifstream file17(input.file17, std::ios_base::in);

    std::vector<std::vector<double>> qa(natom);
    std::vector<std::vector<double>> gfreea(natom);
    std::vector<std::vector<double>> va(natom);
    std::vector<std::vector<double>> evdwa(natom);
    std::vector<std::vector<double>> rvdwa(natom);
    std::vector<std::vector<double>> rhca(natom);
    std::vector<std::vector<double>> xma(natom);

    for (int i = 0; i < natom; ++i)
    {
        for (int j = 0; j < nat[i]; ++j)
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

                if (atp[i][j] == atp_name)
                {
                    qa[i].push_back(xq);
                    gfreea[i].push_back(xfree);
                    va[i].push_back(xvol);
                    evdwa[i].push_back(xevdw);
                    rvdwa[i].push_back(xrvdw);
                    rhca[i].push_back(0.8*xrvdw); //bug? xrhc?
                    xma[i].push_back(xmassa);

                    /// Reset file reading to beggining
                    file17.clear();
                    file17.seekg(0, std::ios::beg);

                    break;
                }
            }
        }
    }

    file17.close();

std::cout << "d\n";
    /// propietats de les boles
    double xmassa = 0.0;

    std::vector<double> xm(natom);
    std::vector<double> qq(natom);
    std::vector<double> vol(natom);
    std::vector<double> gfree(natom);
    std::vector<double> evdw(natom);
    //std::vector<double> sumrhc(natom);
    //std::vector<double> sumrvdw(natom);

    std::vector<double> rvdw(natom);
    std::vector<double> rhc(natom);

    std::vector<std::array<double, 3>> v(natom);

    for (int i = 0; i < natom; ++i)
    {
        xm[i] = 0.0;
        qq[i] = 0.0;
        vol[i] = 0.0;
        gfree[i] = 0.0;
        evdw[i] = 0.0;
        double sumrhc = 0.0;
        double sumrvdw = 0.0;

        if (nat[i] == 1)
        {
            xm[i] = xma[i][0];
            qq[i] = qa[i][0];
            vol[i] = va[i][0];
            gfree[i] = gfreea[i][0];
            evdw[i] = evdwa[i][0];
            rvdw[i] = rvdwa[i][0];
            rhc[i] = 0.8*rvdw[i];
        }
        else
        {
            for (int j = 0; j < nat[i]; ++j)
            {
                xm[i] += xma[i][j];
                qq[i] += qa[i][j];
                vol[i] += va[i][j];
                gfree[i] += gfreea[i][j];
                evdw[i] += evdwa[i][j];
                sumrhc += std::pow(rhca[i][j], 3.0);
                sumrvdw += std::pow(rvdwa[i][j], 3.0);
            }

            rvdw[i] = input.factr * std::pow(sumrvdw, 0.333333);
            rhc[i] = 0.8 * rvdw[i];
        }

        xmassa += xm[i];

        v[i] = {0.0, 0.0, 0.0};
    }

std::cout << "e\n";
    /// llegeix la matriu de topologia (topologia.dat)

    std::vector<std::vector<int>> icov(natom, std::vector<int>(natom));
    std::vector<std::vector<int>> nstep(natom, std::vector<int>(natom));
    std::vector<std::vector<int>> inter(natom, std::vector<int>(natom));
    std::vector<std::vector<int>> istruct(natom, std::vector<int>(natom));

    std::vector<std::array<int, 2>> ibound(NPAIRMAX);
    std::vector<int> rbound(NPAIRMAX);

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

        icov[i][j] = 1;
        ibound[npair][0] = i;
        ibound[npair][1] = j;
        rbound[npair] = rij;

        ++npair;
    }


    /** TODO: Next we write file "dmd.out" with contains all input
        parameters, for now i'm going to skip it
    */
std::cout << "f\n";
    /// reconeix estructura secundaria i estableix ponts d'hidrogen [324]

    std::vector<int> ihb(natom);
    int nhb = 0;

    Xoc xoc;
    xoc.r = r;
    xoc.v = v;
    xoc.xm = xm;
    xoc.rbox = input.rbox;

    for (int i = 0; i < nres - 4; ++i)
    {

        int ii = io[i];

        for (int j = i + 4; j < nres; ++j)
        {
            if (atoms[ica[j]].residue_name != "PRO")
            {
                int jj = ih[j];

                if (ihb[ii] == 0 && ihb[jj] == 0)
                {
                    int n1 = io[i];
                    int n2 = ih[j];

                    double rij1 = dbox(n2, n1, 0, xoc);
                    double rij2 = dbox(n2, n1, 1, xoc);
                    double rij3 = dbox(n2, n1, 2, xoc);

                    double roh = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                    if (roh < distancies.rohmax && roh > distancies.rohmin)
                    {
                        n1 = io[i];
                        n2 = in[j];

                        rij1 = dbox(n2, n1, 0, xoc);
                        rij2 = dbox(n2, n1, 1, xoc);
                        rij3 = dbox(n2, n1, 2, xoc);

                        double rno = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                        if (rno < distancies.rnomax && rno > distancies.rnomin)
                        {
                            n1 = ico[i];
                            n2 = ih[j];

                            rij1 = dbox(n2, n1, 0, xoc);
                            rij2 = dbox(n2, n1, 1, xoc);
                            rij3 = dbox(n2, n1, 2, xoc);

                            double rch = std::sqrt(rij1 * rij1 + rij2 * rij2 + rij3 * rij3);

                            if (rch < distancies.rchmax && rch > distancies.rchmin)
                            {
                                ++nhb;

                                n1 = io[i];
                                n2 = ih[j];

                                ihb[n1] = 1;
                                ihb[n2] = 1;

                                /// write to dmd.out
                                /// write to stdout
                            }
                        }
                    }
                }
            }
        }
    }
std::cout << "g\n";

    Pous pous;
    pous.rstep = std::vector<std::vector<std::array<double, 2>>>(natom);
    pous.estep = std::vector<std::vector<std::array<double, 2>>>(natom);

    Intr intr;

    potential(natom, xoc, pous, intr, cov, pdb, fisic, param, parmsolv);


}
