using Parameters
using Formatting
using Bits

@with_kw mutable struct Distancies{T <: Real}

    rohmin::T = 1.75
    rohmax::T = 2.50
    rnomin::T = 2.75
    rnomax::T = 3.50
    rchmin::T = 2.90
    rchmax::T = 3.75

    roha::T = 2.15
    rohb::T = 2.34
    rnoa::T = 3.10
    rnob::T = 3.20
    rcha::T = 3.27
    rchb::T = 3.40

end

#TODO: parse input instead of commenting the original value and manually hardcodding it
@with_kw mutable struct Input{T <: Real}

    #tsnap::T = 1e-11
    tsnap::T = 1e-12
    tene::T = 1e-12
    tact::T = 2e-14
    sigma::T = 0.05
    temp::T = 300.00
    #kkk::Int32 = 2381
    kkk::Int32 = 2857

    #file7::String = "topologia.dat"
    file7::String = "topcg.dat"
    file8::String = "dmd.out"
    #file9::String = "nativain.pdb"
    file9::String = "structurecg.pdb"
    file10::String = "res"
    #file12::String = "energia.dat"
    file12::String = "energy.dat"
    file15::String = "res"
    #file16::String = "atomtypes.dat"
    file16::String = "beads.dat"
    file17::String = "potentials.dat"
    file19::String = "output.pdb"
    #file20::String = "snapca.pdb"
    file20::String = "snapshots.pdb"
    file21::String = "snapca.pdb"

    nbloc::Int32 = 10000
    rcutgo::T = 8.0
    tmin::T = 1e-30
    dcut::T = 10.0
    ehb::T = 3.0
    ehbc::T = 4.0
    rshake::T = 50.0
    rpot::T = 50.0
    dstep::T = 1e-4
    dijmin::T = 1e-4
    tpush::T = 5e-4
    factm::T = 1.0
    fvdw::T = 8.0
    #fsolv::T = 15.0
    fsolv::T = 12.0
    eps::T = 16.5
    rbox::T = 0.0
    facthc::T = 0.8
    factr::T = 0.9
    rsolv::T = 3.5
    asolv::T = 10.0
    bsolv::T = 0.5
    dwat::T = 6.0

    isolv::Int = 1
    isec::Int = 0
    iterm::Int = 1
    iwr::Int = 0
    icons::Int = 1
    iprint::Int = 1

end

@with_kw mutable struct Constants{T <: Real}

    xlamb::T = 3.5
    a::T = 1e-10
    facte::T = 4186.0

    icm::Int = 0

end

#TODO: parse input info
#function getInput()

let st_rnd = 0.0
    global function getUniformRandom()

        st_rnd += 0.0625

        if (st_rnd > 0.99)
            st_rnd = 0.0625
        end

        return st_rnd

    end
end

function getAtomInRange(natom::Int)::Int

    return rand(1:natom)

end

@with_kw mutable struct Xoc{T <: Real}

    Xoc{T}(natom::Int) where {T<:Real} = 
    begin
        
        #r = Array{T, 2}(undef, natom, 3)
        r = zeros(T, natom, 3)
        #v = Array{T, 2}(undef, natom, 3)
        v = zeros(T, natom, 3)
        #xm = Array{T, 1}(undef, natom)
        xm = zeros(T, natom)

        new(r, v, xm)

    end

    r::Array{T, 2}
    v::Array{T, 2}
    xm::Array{T, 1}
    rbox::T
    ierr::Int

end

@with_kw mutable struct Pous{T <: Real}

    Pous{T}(natom::Int) where {T<:Real} = 
    begin
        
        #rstep = Array{T, 3}(undef, natom, natom, 2)
        rstep = zeros(T, natom, natom, 2)
        #estep = Array{T, 3}(undef, natom, natom, 2)
        estep = zeros(T, natom, natom, 2)

        new(rstep, estep)

    end

    rstep::Array{T, 3}
    estep::Array{T, 3}

end

@with_kw mutable struct Intr

    Intr(natom::Int) =
    begin
        
        #nstep = Array{Int, 2}(undef, natom, natom)
        nstep = zeros(Int, natom, natom)

        #istruct = Array{Int, 2}(undef, natom, natom)
        istruct = zeros(Int, natom, natom)

        #inter = Array{Int, 2}(undef, natom, natom)
        inter = zeros(Int, natom, natom)

        new(nstep, istruct, inter)

    end

    nstep::Array{Int, 2}
    istruct::Array{Int, 2}
    inter::Array{Int, 2}

end

@with_kw mutable struct Cov{T <: Real}

    Cov{T}(natom::Int)  where {T<:Real} = 
    begin
        
        #icov = Array{Int, 2}(undef, natom, natom)
        icov = zeros(Int, natom, natom)
        rbound = Array{T, 1}(undef, 0)
        ibound = Array{Int, 2}(undef, 0, 2)
        #rhc = Array{T, 1}(undef, natom)
        rhc = zeros(T, natom)

        new(icov, rbound, ibound, rhc)

    end

    icov::Array{Int, 2}
    rbound::Array{T, 1}
    ibound::Array{Int, 2}
    rhc::Array{T, 1}
    sigma::T

end

@with_kw mutable struct Pdb

    Pdb(natom::Int) = 
    begin
        
        atom = Array{String, 1}(undef, natom)
        res = Array{String, 1}(undef, natom)
        ind2 = Array{Int, 1}(undef, natom)
        nat = Array{Int, 1}(undef, natom)
        imol = Array{Int, 1}(undef, natom)

        new(atom, res, ind2, nat, imol)
        
    end

    atom::Array{String, 1}
    res::Array{String, 1}
    ind2::Array{Int, 1}
    nat::Array{Int, 1}
    imol::Array{Int, 1}

end

@with_kw mutable struct Atpres

    Atpres(natom::Int) = 
    begin
        
        #ihb = Array{Int, 1}(undef, natom)
        ihb = zeros(Int, natom)
        #ica = Array{Int, 1}(undef, natom)
        ica = zeros(Int, natom)
        #io = Array{Int, 1}(undef, natom)
        io = zeros(Int, natom)
        #ih = Array{Int, 1}(undef, natom)
        ih = zeros(Int, natom)
        #ico = Array{Int, 1}(undef, natom)
        ico = zeros(Int, natom)
        #in = Array{Int, 1}(undef, natom)
        in = zeros(Int, natom)

        new(ihb, ica, io, ih, ico, in)

    end

    ihb::Array{Int, 1}
    ica::Array{Int, 1}
    io::Array{Int, 1}
    ih::Array{Int, 1}
    ico::Array{Int, 1}
    in::Array{Int, 1}

end

@with_kw mutable struct Shake

    Shake(natom::Int) = 
    begin
        
        #ishk = Array{Int, 1}(undef, natom)
        ishk = zeros(Int, natom)

        #nshk = Array{Int, 2}(undef, natom, natom)
        nshk = zeros(Int, natom, natom)

        new(ishk, nshk)

    end

    ishk::Array{Int, 1}
    nshk::Array{Int, 2}

end

mutable struct Npt

    Npt(natom::Int) = 
    begin

        #ipot = Array{Int, 1}(undef, natom)
        ipot = zeros(Int, natom)

        #npot = Array{Int, 2}(undef, natom, natom)
        npot = zeros(Int, natom, natom)

        new(ipot, npot)

    end

    ipot::Array{Int, 1}
    npot::Array{Int, 2}

end

mutable struct Fisic{T <: Real}

    Fisic{T}(natom::Int)  where {T<:Real} = 
    begin

        #evdw = Array{T, 1}(undef, natom)
        evdw = zeros(T, natom)
        #rvdw = Array{T, 1}(undef, natom)
        rvdw = zeros(T, natom)
        #qq = Array{T, 1}(undef, natom)
        qq = zeros(T, natom)
        #gfree = Array{T, 1}(undef, natom)
        gfree = zeros(T, natom)
        #vol = Array{T, 1}(undef, natom)
        vol = zeros(T, natom)

        new(evdw, rvdw, qq, gfree, vol)

    end

    evdw::Array{T, 1}
    rvdw::Array{T, 1}
    qq::Array{T, 1}
    gfree::Array{T, 1}
    vol::Array{T, 1}

end

mutable struct Param{T <: Real}

    Param{T}() where {T<:Real} = 
    begin
        
        new(0.0, 0.0, 0.0, 0.0)

    end

    fvdw::T
    fsolv::T
    eps::T
    xlamb::T

end

mutable struct Parmsolv{T <: Real}

    Parmsolv{T}(natom::Int) where {T<:Real} = 
    begin

        #icont = Array{Int, 1}(undef, natom)
        icont = zeros(Int, natom)
        #fcont = Array{T, 1}(undef, natom)
        fcont = zeros(T, natom)

        new(icont, fcont)

    end

    icont::Array{Int, 1}
    fcont::Array{T, 1}

    rsolv::T
    asolv::T
    bsolv::T
    dwat::T

end

mutable struct Other{T <: Real}

    Other{T}(natom::Int)  where {T<:Real} = 
    begin

        #qa = fill(Array{T}(undef, 0), natom)
        qa = [[] for _ in 1 : natom]
        #gfreea = fill(Array{T}(undef, 0), natom)
        gfreea = [[] for _ in 1 : natom]
        #va = fill(Array{T}(undef, 0), natom)
        va = [[] for _ in 1 : natom]
        #evdwa = fill(Array{T}(undef, 0), natom)
        evdwa = [[] for _ in 1 : natom]
        #rvdwa = fill(Array{T}(undef, 0), natom)
        rvdwa = [[] for _ in 1 : natom]
        #rhca = fill(Array{T}(undef, 0), natom)
        rhca = [[] for _ in 1 : natom]
        #xma = fill(Array{T}(undef, 0), natom)
        xma = [[] for _ in 1 : natom]

        #rant = Array{T, 2}(undef, natom, 3)
        rant = zeros(T, natom, 3)

        ind1 = Array{Int, 1}(undef, natom)

        #rcm = Array{T, 1}(undef, 3)
        rcm = zeros(T, 3)
        vcm = Array{T, 1}(undef, 3)

        #ireg = Array{Int, 2}(undef, natom, natom)
        ireg = zeros(Int, natom, natom)

        cad = Array{String, 1}(undef, natom)

        #atp = fill(Array{String}(undef, 0), natom)
        atp = [[] for _ in 1 : natom]

        new(qa, gfreea, va, evdwa, rvdwa, rhca, xma, rant, ind1, rcm, vcm, ireg, cad, atp)

    end

    qa::Array{Array{T, 1}, 1}
    gfreea::Array{Array{T, 1}, 1}
    va::Array{Array{T, 1}, 1}
    evdwa::Array{Array{T, 1}, 1}
    rvdwa::Array{Array{T, 1}, 1}
    rhca::Array{Array{T, 1}, 1}
    xma::Array{Array{T, 1}, 1}

    rant::Array{T, 2}

    ind1::Array{Int, 1}

    rcm::Array{T, 1}
    vcm::Array{T, 1}

    ireg::Array{Int, 2}

    cad::Array{String, 1}

    atp::Array{Array{String}, 1}   #TODO: this should be array of arrays instead of matrix, check if other variables are the same

end

function dbox(n1, n2, k, xoc::Xoc)

    rbox2 = 0.5 * xoc.rbox
    r12 = xoc.r[n2, k] - xoc.r[n1, k]

    if (r12 > rbox2)
        r12 += -xoc.rbox
    elseif (r12 < -rbox2)
        r12 += xoc.rbox

    end

    return r12

end

function potencial(natom, xoc::Xoc, pous::Pous, intr::Intr, cov::Cov, pdb::Pdb, fisic::Fisic, param::Param, parmsolv::Parmsolv)

    rsolv2 = parmsolv.rsolv ^ 2.0
    rbox2 = 0.5 * xoc.rbox
    nres = pdb.ind2[natom]
    dwatd = parmsolv.dwat / sqrt(3.0)

    for i in 1 : natom
        r =
        [
            [xoc.r[i, 1] + parmsolv.dwat, xoc.r[i, 2]                , xoc.r[i, 3]                ],
            [xoc.r[i, 1] - parmsolv.dwat, xoc.r[i, 2]                , xoc.r[i, 3]                ],
            [xoc.r[i, 1]                , xoc.r[i, 2] + parmsolv.dwat, xoc.r[i, 3]                ],
            [xoc.r[i, 1]                , xoc.r[i, 2] - parmsolv.dwat, xoc.r[i, 3]                ],
            [xoc.r[i, 1]                , xoc.r[i, 2]                , xoc.r[i, 3] + parmsolv.dwat],
            [xoc.r[i, 1]                , xoc.r[i, 2]                , xoc.r[i, 3] - parmsolv.dwat]
        ]

        rv =
        [
            [xoc.r[i, 1] + dwatd, xoc.r[i, 2] + dwatd, xoc.r[i, 3] + dwatd],
            [xoc.r[i, 1] + dwatd, xoc.r[i, 2] + dwatd, xoc.r[i, 3] - dwatd],
            [xoc.r[i, 1] + dwatd, xoc.r[i, 2] - dwatd, xoc.r[i, 3] + dwatd],
            [xoc.r[i, 1] + dwatd, xoc.r[i, 2] - dwatd, xoc.r[i, 3] - dwatd],
            [xoc.r[i, 1] - dwatd, xoc.r[i, 2] + dwatd, xoc.r[i, 3] + dwatd],
            [xoc.r[i, 1] - dwatd, xoc.r[i, 2] + dwatd, xoc.r[i, 3] - dwatd],
            [xoc.r[i, 1] - dwatd, xoc.r[i, 2] - dwatd, xoc.r[i, 3] + dwatd],
            [xoc.r[i, 1] - dwatd, xoc.r[i, 2] - dwatd, xoc.r[i, 3] - dwatd]
        ]

        parmsolv.icont[i] = 0

        for r_ind in 1 : length(r)
            for j in 1 : natom
                if (pdb.imol[i] != pdb.imol[j])
                    continue
                end

                if (i == j)
                    continue
                end

                rij1 = r[r_ind][1] - xoc.r[j, 1]

                if (rij1 > rbox2)
                    rij1 += -xoc.rbox
                elseif (rij1 < -rbox2)
                    rij1 += xoc.rbox
                end

                rij2 = r[r_ind][2] - xoc.r[j, 2]

                if (rij2 > rbox2)
                    rij2 += -xoc.rbox
                elseif (rij2 < -rbox2)
                    rij2 += xoc.rbox
                end

                rij3 = r[r_ind][3] - xoc.r[j, 3]

                if (rij3 > rbox2)
                    rij3 += -xoc.rbox
                elseif (rij3 < -rbox2)
                    rij3 += xoc.rbox
                end

                rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

                if (rmod2 < rsolv2)
                    parmsolv.icont[i] += 1
                    break
                end
            end
        end

        for rv_ind in 1 : length(rv)
            for j in 1 : natom
                if (pdb.imol[i] != pdb.imol[j])
                    continue
                end

                if (i == j)
                    continue
                end

                rij1 = rv[rv_ind][1] - xoc.r[j, 1]

                if (rij1 > rbox2)
                    rij1 += -xoc.rbox
                elseif (rij1 < -rbox2)
                    rij1 += xoc.rbox
                end

                rij2 = rv[rv_ind][2] - xoc.r[j, 2]

                if (rij2 > rbox2)
                    rij2 += -xoc.rbox
                elseif (rij2 < -rbox2)
                    rij2 += xoc.rbox
                end

                rij3 = rv[rv_ind][3] - xoc.r[j, 3]

                if (rij3 > rbox2)
                    rij3 += -xoc.rbox
                elseif (rij3 < -rbox2)
                    rij3 += xoc.rbox
                end

                rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

                if (rmod2 < rsolv2)
                    parmsolv.icont[i] += 1
                    break
                end
            end
        end
    end

    for i in 1 : natom

        parmsolv.fcont[i] = 1.0 / (1.0 + exp((parmsolv.icont[i] - parmsolv.asolv) / parmsolv.bsolv))

    end

    for i in 1 : (natom - 1)

        for j in (i + 1) : natom

            if (cov.icov[i, j] != 0)
                continue
            end

            if (intr.istruct[i, j] != 0)
                continue
            end

            rvdwij = fisic.rvdw[i] + fisic.rvdw[j]
            sto = (2.0 / (pdb.nat[i] ^ 0.33 + pdb.nat[j] ^ 0.33)) ^ 6.0
            potvdw = sqrt(fisic.evdw[i] * fisic.evdw[j]) * sto * (sto - 2.0)
            potlk = -0.09 / param.xlamb * (fisic.gfree[i] * fisic.vol[j] + fisic.gfree[j] * fisic.vol[i]) / (rvdwij ^ 2.0 * exp((rvdwij / param.xlamb) ^ 2.0))
            eij = param.fvdw * potvdw + param.fsolv * potlk * parmsolv.fcont[i] * parmsolv.fcont[j] + param.eps * fisic.qq[i] * fisic.qq[j] / rvdwij

            intr.nstep[i, j] = 2
            pous.rstep[i, j, 1] = 0.9 * rvdwij
            pous.rstep[i, j, 2] = 1.1 * rvdwij

            if (eij < 0.0)
                pous.estep[i, j, 1] = 3.0 * eij
                pous.estep[i, j, 2] = -eij
            else
                pous.estep[i, j, 1] = -eij
                pous.estep[i, j, 2] = -eij
            end

        end

    end

    return

end

function enchufa(natom, dcut, xoc::Xoc, intr::Intr, cov::Cov, atpres::Atpres, pdb::Pdb, npt::Npt)

    dcut2 = dcut ^ 2.0

    for i in 1 : (natom - 1)
        for l in 1 : npt.ipot[i]

            j = npt.npot[i, l]

            if (pdb.nat[i] > 1 && pdb.nat[j] > 1)

                rij1 = dbox(i, j, 1, xoc)
                rij2 = dbox(i, j, 2, xoc)
                rij3 = dbox(i, j, 3, xoc)

                rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0
                
                if (rmod2 < dcut2)
                    intr.inter[i, j] = 1
                end
            end
        end
    end

    nres = pdb.ind2[natom]

    for i in 1 : (nres - 1)

        n1 = atpres.in[i]
        n2 = atpres.in[i + 1]

        intr.inter[n1, n2] = 0

    end

    for i in 2 : nres

        n1 = atpres.ico[i - 1]
        n2 = atpres.ico[i]

        intr.inter[n1, n2] = 0

    end

    return

end

function creapouhb(n1, n2, rmin, r0, r1, rmax, ehb, pous::Pous, intr::Intr)

    intr.inter[n1, n2] = 1
    intr.istruct[n1, n2] = 1
    intr.nstep[n1, n2] = 2

    pous.rstep[n1, n2, 1] = rmin
    pous.rstep[n1, n2, 2] = rmax

    pous.estep[n1, n2, 1] = -1.5 * ehb
    pous.estep[n1, n2, 2] = 1.5 * ehb

    #println(n1, " ", n2, " ", rmin, " ", rmax, " ", ehb)
    #readline()

    return

end

function chgmom(mem1, mem2, rij1, rij2, rij3, xoc::Xoc)

    vdmod = 0.0

    vdmod += (xoc.v[mem2, 1] - xoc.v[mem1, 1]) * rij1
    vdmod += (xoc.v[mem2, 2] - xoc.v[mem1, 2]) * rij2
    vdmod += (xoc.v[mem2, 3] - xoc.v[mem1, 3]) * rij3

    rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

    vdmod = vdmod / rmod2

    xsum = 0.5 * (1.0 / xoc.xm[mem1] + 1.0 / xoc.xm[mem2])

    dp = vdmod / xsum

    xoc.v[mem1, 1] = xoc.v[mem1, 1] + dp / xoc.xm[mem1] * rij1
    xoc.v[mem2, 1] = xoc.v[mem2, 1] - dp / xoc.xm[mem2] * rij1

    xoc.v[mem1, 2] = xoc.v[mem1, 2] + dp / xoc.xm[mem1] * rij2
    xoc.v[mem2, 2] = xoc.v[mem2, 2] - dp / xoc.xm[mem2] * rij2

    xoc.v[mem1, 3] = xoc.v[mem1, 3] + dp / xoc.xm[mem1] * rij3
    xoc.v[mem2, 3] = xoc.v[mem2, 3] - dp / xoc.xm[mem2] * rij3

    return

end

function chgmomene(mem1, mem2, rij1, rij2, rij3, dpot, ich, xoc::Xoc)

    a = 1e-10

    vdmod = 0.0

    vdmod = vdmod + (xoc.v[mem2, 1] - xoc.v[mem1, 1]) * rij1
    vdmod = vdmod + (xoc.v[mem2, 2] - xoc.v[mem1, 2]) * rij2
    vdmod = vdmod + (xoc.v[mem2, 3] - xoc.v[mem1, 3]) * rij3

    rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

    vdmod /= rmod2

    xsum = 0.5 * (1.0 / xoc.xm[mem1] + 1.0 / xoc.xm[mem2])

    dp = vdmod / xsum
    sto = dp ^ 2.0 / 4.0 - dpot / (rmod2 * xsum * a ^ 2.0)

    if (sto > 0.0)
        if (vdmod > 0.0)
            dp = dp / 2.0 - sqrt(sto)
            
        else
            dp = dp / 2.0 + sqrt(sto)

        end

        ich[] = 1

    else
        ich[] = 0

    end

    xoc.v[mem1, 1] = xoc.v[mem1, 1] + (dp / xoc.xm[mem1]) * rij1
    xoc.v[mem2, 1] = xoc.v[mem2, 1] - (dp / xoc.xm[mem2]) * rij1

    xoc.v[mem1, 2] = xoc.v[mem1, 2] + (dp / xoc.xm[mem1]) * rij2
    xoc.v[mem2, 2] = xoc.v[mem2, 2] - (dp / xoc.xm[mem2]) * rij2

    xoc.v[mem1, 3] = xoc.v[mem1, 3] + (dp / xoc.xm[mem1]) * rij3
    xoc.v[mem2, 3] = xoc.v[mem2, 3] - (dp / xoc.xm[mem2]) * rij3

    #println(dpot, " ", mem1, " ", mem2, " ", rmod2)
    #readline()

    return

end

#TODO: delete unused vars
function dmdshake(natom, nbound, mem1, mem2, xoc::Xoc, cov::Cov, shake::Shake)

    ierr = 0

    for i in 1 : (natom - 1)
        for l in 1 : shake.ishk[i]
            j = shake.nshk[i, l]

            rij1 = dbox(j, i, 1, xoc)
            rij2 = dbox(j, i, 2, xoc)
            rij3 = dbox(j, i, 3, xoc)

            rij = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

            vij1 = xoc.v[i, 1] - xoc.v[j, 1]
            vij2 = xoc.v[i, 2] - xoc.v[j, 2]
            vij3 = xoc.v[i, 3] - xoc.v[j, 3]

            prod = rij1 * vij1 + rij2 * vij2 + rij3 * vij3

            dmin = cov.rhc[i] + cov.rhc[j]

            if (rij < dmin && prod < 0.0)
                chgmom(i, j, rij1, rij2, rij3, xoc)

                ierr += 1

            end

        end
    end

    for k in 1 : nbound
        i = cov.ibound[k, 1]
        j = cov.ibound[k, 2]

        rbmin = cov.rbound[k] * (1.0 - cov.sigma)
        rbmax = cov.rbound[k] * (1.0 + cov.sigma)

        rbmin2 = rbmin ^ 2.0
        rbmax2 = rbmax ^ 2.0

        rij1 = dbox(j, i, 1, xoc)
        rij2 = dbox(j, i, 2, xoc)
        rij3 = dbox(j, i, 3, xoc)

        rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

        vij1 = xoc.v[i, 1] - xoc.v[j, 1]
        vij2 = xoc.v[i, 2] - xoc.v[j, 2]
        vij3 = xoc.v[i, 3] - xoc.v[j, 3]

        prod = rij1 * vij1 + rij2 * vij2 + rij3 * vij3

        if (rmod2 > rbmax2 && prod > 0.0)
            ierr += 1
            chgmom(i, j, rij1, rij2, rij3, xoc)

        elseif (rmod2 < rbmin2 && prod < 0.0)
            ierr += 1
            chgmom(i, j, rij1, rij2, rij3, xoc)

        end
    end
end

function rnd_gauss(fi, xm, T)

    R = 8.314

    std_dev = sqrt((T * R) / xm)
    pi = 4.0 * atan(1.0)

    rnd1 = getUniformRandom()
    rnd2 = getUniformRandom()

    fi = std_dev * sqrt(-2.0 * log(rnd1) * cos(2.0 * pi * rnd2))

end

struct Atom
    record_type::String
    serial_number::Int
    name::String
    residue_name::String
    chain_identifier::String
    residue_sequence_number::Int
    x_orthogonal_coord::Float64
    y_orthogonal_coord::Float64
    z_orthogonal_coord::Float64
end

function get_molecule_info(molecule_file_path::String, input::Input)

    atoms = Array{Atom, 1}(undef, 0)

    file = open(molecule_file_path, "r")

    for atom_record in eachline(file)

        atom_record = split(atom_record, r" +")

        record_type = atom_record[1]
        serial_number = parse(Int, atom_record[2])
        name = atom_record[3]
        residue_name = atom_record[4]
        chain_identifier = atom_record[5]
        residue_sequence_number = parse(Int, atom_record[6])
        x_orthogonal_coord = parse(Float64, atom_record[7])
        y_orthogonal_coord = parse(Float64, atom_record[8])
        z_orthogonal_coord = parse(Float64, atom_record[9])
        
        curr_atom = Atom(record_type, serial_number, name, residue_name, chain_identifier, residue_sequence_number, x_orthogonal_coord, y_orthogonal_coord, z_orthogonal_coord)

        push!(atoms, curr_atom)
        #atoms = [atoms, curr_atom]

    end

    close(file)

    natom = length(atoms)

    pdb = Pdb(natom)
    xoc = Xoc{Float64}(natom)
    xoc.rbox = input.rbox
    other = Other{Float64}(natom)

    for i in 1 : natom

        pdb.atom[i] = atoms[i].name
        pdb.res[i] = atoms[i].residue_name

        other.cad[i] = atoms[i].chain_identifier
        other.ind1[i] = atoms[i].residue_sequence_number

        #xoc.r[i] = [atoms[i].x_orthogonal_coord, atoms[i].y_orthogonal_coord, atoms[i].z_orthogonal_coord]
        xoc.r[i, 1] = atoms[i].x_orthogonal_coord
        xoc.r[i, 2] = atoms[i].y_orthogonal_coord
        xoc.r[i, 3] = atoms[i].z_orthogonal_coord

    end

    return natom, pdb, xoc, other

end

function set_residue_indexes(natom::Int, pdb::Pdb, other::Other)

    kk = 0
    im = 0

    atpres = Atpres(natom)

    for n in 1 : natom

        if (n > 1 && other.ind1[n] < other.ind1[n - 1])
            kk += other.ind1[n - 1]
        end

        if (n == 1 || other.cad[n] != other.cad[n - 1])
            im += 1
        end

        pdb.ind2[n] = other.ind1[n] + kk

        pdb.imol[n] = im

        k1 = pdb.ind2[n]

        if (pdb.atom[n] == "N")
            atpres.in[k1] = n
        elseif (pdb.atom[n] == "H")
            atpres.ih[k1] = n
        elseif (pdb.atom[n] == "CA")
            atpres.ica[k1] = n
        elseif (pdb.atom[n] == "C")
            atpres.ico[k1] = n
        elseif (pdb.atom[n] == "O")
            atpres.io[k1] = n
        end

    end

    nres = pdb.ind2[natom]

    return nres, atpres

end

function get_atom_types(atom_types_file_path::String, natom::Int, pdb::Pdb, other::Other)

    for i in 1 : natom

        pdb.nat[i] = 0

        name_atoms = ["N", "H", "C", "O", "OXT"]

        if (pdb.atom[i] in name_atoms)

            pdb.nat[i] = 1

            if     (pdb.atom[i] == "N") push!(other.atp[i], "nh")
            elseif (pdb.atom[i] == "H") push!(other.atp[i], "h")
            elseif (pdb.atom[i] == "C") push!(other.atp[i], "co")
            elseif (pdb.atom[i] == "O" || pdb.atom[i] == "OXT") push!(other.atp[i], "oc")
            end

            continue

        end

        file16 = open(atom_types_file_path, "r")

        for atom_type_record in eachline(file16)

            atom_type_record = split(atom_type_record, r" +")

            name = atom_type_record[1]
            residue_name = atom_type_record[2]
            atp_name = atom_type_record[3]

            if (pdb.atom[i] == name && pdb.res[i] == residue_name)

                push!(other.atp[i], atp_name)
                pdb.nat[i] += 1

            end

        end

        close(file16)

    end
    

end

function get_atom_parameters(atom_parameters_file_path::String, natom::Int, pdb::Pdb, other::Other)

    for i in 1 : natom

        for j in 1 : pdb.nat[i]

            file17 = open(atom_parameters_file_path, "r")

            for atom_params_record in eachline(file17)

                atom_params_record = split(atom_params_record, r" +")

                atp_name = atom_params_record[1]
                xq = parse(Float64, atom_params_record[2])
                xfree = parse(Float64, atom_params_record[3])
                xvol = parse(Float64, atom_params_record[4])
                xevdw = parse(Float64, atom_params_record[5])
                xrvdw = parse(Float64, atom_params_record[6])
                xrhc = parse(Float64, atom_params_record[7])
                xmassa = parse(Float64, atom_params_record[8])

                if (other.atp[i][j] == atp_name)

                    push!(other.qa[i], xq)
                    push!(other.gfreea[i], xfree)
                    push!(other.va[i], xvol)
                    push!(other.evdwa[i], xevdw)
                    push!(other.rvdwa[i], xrvdw)
                    push!(other.rhca[i], 0.8 * xrvdw)
                    push!(other.xma[i], xmassa)

                    break

                end

            end

            close(file17)

        end

    end

end

function set_atom_parameters(natom::Int, xoc::Xoc, pdb::Pdb, other::Other, input::Input)

    xmassa = 0.0

    cov = Cov{Float64}(natom)
    cov.sigma = input.sigma
    fisic = Fisic{Float64}(natom)

    for i in 1 : natom

        xoc.xm[i] = 0.0
        fisic.qq[i] = 0.0
        fisic.vol[i] = 0.0
        fisic.gfree[i] = 0.0
        fisic.evdw[i] = 0.0

        sumrhc = 0.0
        sumrvdw = 0.0

        if (pdb.nat[i] == 1)

            xoc.xm[i] = other.xma[i][1]

            fisic.qq[i] = other.qa[i][1]
            fisic.vol[i] = other.va[i][1]
            fisic.gfree[i] = other.gfreea[i][1]
            fisic.evdw[i] = other.evdwa[i][1]
            fisic.rvdw[i] = other.rvdwa[i][1]

            cov.rhc[i] = 0.8 * fisic.rvdw[i]

        else

            for j in 1 : pdb.nat[i]

                xoc.xm[i] += other.xma[i][j]

                fisic.qq[i] += other.qa[i][j]
                fisic.vol[i] += other.va[i][j]
                fisic.gfree[i] += other.gfreea[i][j]
                fisic.evdw[i] += other.evdwa[i][j]

                sumrhc += other.rhca[i][j] ^ 3.0
                sumrvdw += other.rvdwa[i][j] ^ 3.0

            end

            fisic.rvdw[i] = input.factr * sumrvdw ^ 0.3333

            cov.rhc[i] = 0.8 * fisic.rvdw[i]

        end

        xmassa += xoc.xm[i]

        xoc.v[i, 1] = 0.0
        xoc.v[i, 2] = 0.0
        xoc.v[i, 3] = 0.0

    end

    return xmassa, cov, fisic

end

function get_topology_matrix(topology_file_path::String, cov::Cov)

    file7 = open(topology_file_path, "r")

    npair = 0

    for topology_record in eachline(file7)

        topology_record = split(topology_record, r" +")
        filter!(!isempty, topology_record)

        i = parse(Int, topology_record[1])
        j = parse(Int, topology_record[2])
        rij = parse(Float64, topology_record[3])

        cov.icov[i, j] = 1

        ibound_ij = [i, j]'
        cov.ibound = vcat(cov.ibound, ibound_ij)

        append!(cov.rbound, rij)

        npair += 1

    end

    close(file7)

    return npair

end

#TODO: check why it doesn't work
function stablish_hydrogen_bonds(nres::Int, xoc::Xoc, pdb::Pdb, atpres::Atpres, distancies::Distancies)

    nhb = 0

    for i in 1 : (nres - 4)

        ii = atpres.io[i]

        for j in (i + 4) : nres

            if (pdb.res[atpres.ica[j]] == "PRO")
                continue
            end

            jj = atpres.ih[j]

            if (atpres.ihb[ii] != 0 || atpres.ihb[jj] != 0)
                continue
            end

            n1 = atpres.io[i]
            n2 = atpres.ih[j]

            rij1 = dbox(n2, n1, 1, xoc)
            rij2 = dbox(n2, n1, 2, xoc)
            rij3 = dbox(n2, n1, 3, xoc)

            roh = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

            if (roh > distancies.rohmax || roh < distancies.rohmin)
                continue
            end

            n1 = atpres.io[i]
            n2 = atpres.in[j]

            rij1 = dbox(n2, n1, 1, xoc)
            rij2 = dbox(n2, n1, 2, xoc)
            rij3 = dbox(n2, n1, 3, xoc)

            rno = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

            if (rno > distancies.rnomax || rno < distancies.rnomin)
                continue
            end

            n1 = atpres.ico[i]
            n2 = atpres.ih[j]

            rij1 = dbox(n2, n1, 1, xoc)
            rij2 = dbox(n2, n1, 2, xoc)
            rij3 = dbox(n2, n1, 3, xoc)

            rch = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

            if (rch > distancies.rchmax || rch < distancies.rchmin)
                continue
            end

            nhb += 1

            n1 = atpres.io[i]
            n2 = atpres.ih[j]

            atpres.ihb[n1] = 1
            atpres.ihb[n2] = 1

            #TODO: write to dmd.out
            println("HBOND ", pdb.atom[ii], " ", pdb.res[ii], " ", pdb.ind2[ii], " ", pdb.atom[jj], " ", pdb.res[jj], " ", pdb.ind2[jj])

        end

    end

    for i in 1 : (nres - 4)

        if (pdb.res[atpres.ica[i]] == "PRO")
            continue
        end

        ii = atpres.ih[i]

        for j in (i + 4) : nres

            jj = atpres.io[j]

            if (atpres.ihb[ii] != 0 || atpres.ihb[jj] != 0)
                continue
            end

            n1 = atpres.ih[i]
            n2 = atpres.io[j]

            rij1 = dbox(n2, n1, 1, xoc)
            rij2 = dbox(n2, n1, 2, xoc)
            rij3 = dbox(n2, n1, 3, xoc)

            roh = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

            if (roh > distancies.rohmax || roh < distancies.rohmin)
                continue
            end

            n1 = atpres.in[i]
            n2 = atpres.io[j]

            rij1 = dbox(n2, n1, 1, xoc)
            rij2 = dbox(n2, n1, 2, xoc)
            rij3 = dbox(n2, n1, 3, xoc)

            rno = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

            if (rno > distancies.rnomax || rno < distancies.rnomin)
                continue
            end

            n1 = atpres.ih[i]
            n2 = atpres.ico[j]

            rij1 = dbox(n2, n1, 1, xoc)
            rij2 = dbox(n2, n1, 2, xoc)
            rij3 = dbox(n2, n1, 3, xoc)

            rch = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

            if (rch > distancies.rchmax || rch < distancies.rchmin)
                continue
            end

            nhb += 1

            n1 = atpres.ih[i]
            n2 = atpres.io[j]

            atpres.ihb[n1] = 1
            atpres.ihb[n2] = 1

            #TODO: write to dmd.out
            println("HBOND ", pdb.atom[ii], " ", pdb.res[ii], " ", pdb.ind2[ii], " ", pdb.atom[jj], " ", pdb.res[jj], " ", pdb.ind2[jj])

        end

    end

    return nhb

end

function assign_initial_overlaps(natom::Int, rshake2, rpot2, xoc::Xoc, intr::Intr, cov::Cov, shake::Shake, npt::Npt)

    for i in 1 : (natom - 1)
        
        shake.ishk[i] = 0
        npt.ipot[i] = 0

        for j in (i + 1) : natom

            if (cov.icov[i, j] != 0)
                continue
            end

            rij1 = dbox(i, j, 1, xoc)
            rij2 = dbox(i, j, 2, xoc)
            rij3 = dbox(i, j, 3, xoc)

            rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

            if (rmod2 < rshake2)
            
                shake.ishk[i] += 1
                k = shake.ishk[i]
                shake.nshk[i, k] = j
                
            end

            if (intr.istruct[i, j] != 1 && rmod2 < rpot2)

                npt.ipot[i] += 1
                k = npt.ipot[i]
                npt.npot[i, k] = j

            end

        end

    end

end

function assign_initial_interaction_regions(natom::Int, xoc::Xoc, intr::Intr, pous::Pous, cov::Cov, other::Other)

    for i in 1 : (natom - 1)
        for j in (i + 1) : natom

            if (cov.icov[i, j] != 0)
                continue
            end
            
            rij1 = dbox(j, i, 1, xoc)
            rij2 = dbox(j, i, 2, xoc)
            rij3 = dbox(j, i, 3, xoc)

            rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

            rij = sqrt(rmod2)

            k = 1
            while (k <= intr.nstep[i, j] && rij > pous.rstep[i, j, k])
                k += 1
            end

            other.ireg[i, j] = k

        end
    end

end

function obtain_initial_conformation_potential_energy(natom::Int, xoc::Xoc, intr::Intr, pous::Pous, pdb::Pdb)

    epot0 = 0.0
    epotmol0 = 0.0
    epothb0 = 0.0
    epothbmol0 = 0.0

    for i in 1 : (natom - 1)
        for j in (i + 1) : natom

            rij1 = dbox(i, j, 1, xoc)
            rij2 = dbox(i, j, 2, xoc)
            rij3 = dbox(i, j, 3, xoc)

            rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

            dist = sqrt(rmod2)

            if (intr.inter[i, j] != 1)
                continue
            end

            k = intr.nstep[i, j]

            while (k > 0 && dist < pous.rstep[i, j, k])
                epot0 -= pous.estep[i, j, k]

                if (pdb.imol[i] != pdb.imol[j])
                    epotmol0 -= pous.estep[i, j, k]
                end

                if (intr.istruct[i, j] == 1)
                    epothb0 -= pous.estep[i, j, k]

                    if (pdb.imol[i] != pdb.imol[j])
                        epothbmol0 -= pous.estep[i, j, k]
                    end
                end

                k -= 1
            end

        end
    end

    return epot0, epotmol0, epothb0, epothbmol0

end

function adjust_initial_velocities_and_energies(natom::Int, xmassa, epot0, xoc::Xoc, other::Other, input::Input, constants::Constants)

    for j in 1 : 3

        other.vcm[j] = 0.0

        for i in 1 : natom

            xoc.v[i, j] = getUniformRandom()
            other.vcm[j] = other.vcm[j] + xoc.xm[i] * xoc.v[i, j]

        end

        other.vcm[j] /= xmassa

    end

    ekin = 0.0

    for j in 1 : 3
        for i in 1 : natom

            xoc.v[i, j] -= other.vcm[j]
            ekin = ekin + 0.5 * xoc.xm[i] * (xoc.v[i, j] * constants.a) ^ 2.0

        end
    end

    sto = 1.5 * 8.314 * natom * input.temp / ekin

    ekin0 = 0.0

    for j in 1 : 3
        for i in 1 : natom

            xoc.v[i, j] *= sqrt(sto)
            ekin0 = ekin0 + 0.5 * xoc.xm[i] * ((xoc.v[i, j] * constants.a) ^ 2.0)

        end
    end

    ekin0 /= constants.facte

    etot0 = epot0 + ekin0

    for j in 1 : 3

        other.rcm[j] = 0.0

        for i in 1 : natom

            other.rcm[j] = other.rcm[j] + xoc.xm[i] * xoc.r[i, j]

        end

        other.rcm[j] /= xmassa

    end

    return ekin, sto, ekin0, etot0

end

function main()

    #TODO: parse input
    input = Input()
    distancies = Distancies()
    constants = Constants()

    if (input.tene > input.tsnap)
        input.tene = input.tsnap
    end

    if (input.tact > input.tene)
        input.tact = input.tene
    end

    if (input.rbox < 1E-10)
        constants.icm = 1
        input.rbox = 300.0
    end

    rbox2 = 0.5 * input.rbox
    rshake2 = input.rshake * input.rshake
    rpot2 = input.rpot * input.rpot

    natom, pdb, xoc, other = get_molecule_info(input.file9, input)

    nres, atpres = set_residue_indexes(natom, pdb, other)

    get_atom_types(input.file16, natom, pdb, other)

    get_atom_parameters(input.file17, natom, pdb, other)

    xmassa, cov, fisic = set_atom_parameters(natom, xoc, pdb, other, input)

    npair = get_topology_matrix(input.file7, cov)

    nhb = stablish_hydrogen_bonds(nres, xoc, pdb, atpres, distancies)

    pous = Pous{Float64}(natom)
    intr = Intr(natom)

    param = Param{Float64}()
    param.fvdw = input.fvdw
    param.fsolv = input.fsolv
    param.eps = input.eps
    param.xlamb = constants.xlamb

    parmsolv = Parmsolv{Float64}(natom)
    parmsolv.rsolv = input.rsolv
    parmsolv.asolv = input.asolv
    parmsolv.bsolv = input.bsolv
    parmsolv.dwat = input.dwat

    potencial(natom, xoc, pous, intr, cov, pdb, fisic, param, parmsolv)

    shake = Shake(natom)
    npt = Npt(natom)

    assign_initial_overlaps(natom, rshake2, rpot2, xoc, intr, cov, shake, npt)

    if (input.isolv != 0)
        enchufa(natom, input.dcut, xoc, intr, cov, atpres, pdb, npt)
    end

    assign_initial_interaction_regions(natom, xoc, intr, pous, cov, other)

    epot0, epotmol0, epothb0, epothbmol0 = obtain_initial_conformation_potential_energy(natom, xoc, intr, pous, pdb)

    ekin, sto, ekin0, etot0 = adjust_initial_velocities_and_energies(natom, xmassa, epot0, xoc, other, input, constants)

    #ibloc = 0

    #TODO: DELETE SECTION BEGIN
    #printfmtln("rbox2: {:.10f} rshake2: {:.10f} rpot2: {:.10f}", rbox2, rshake2, rpot2)
    #println("natom: ", natom, " nres: ", nres, " npair: ", npair, " nhb: ", nhb)

    #printfmtln("xmassa: {:.10f}", xmassa)
    #println("xmassa: ", xmassa)

    #printfmtln("epot0: {:.15f} epotmol0: {:.15f} epothb0: {:.15f} epothbmol0: {:.15f}", epot0, epotmol0, epothb0, epothbmol0)
    #println("epot0: ", epot0, " epotmol0: ", epotmol0, " epothb0: ", epothb0, " epothbmol0: ", epothbmol0)

    #printfmtln("ekin: {:.15f} sto: {:.15f} ekin0: {:.15f} etot0: {:.15f}", ekin, sto*1e-25, ekin0, etot0)
    #println("ekin: ", ekin, " sto: ", sto, " ekin0: ", ekin0, " etot0: ", etot0)

    #println(bits(epot0))
    #TODO: DELETE SECTION END

    file_input = open("input.pdb", "w")
    file20 = open(input.file20, "w")
    file21 = open(input.file21, "w")

    printfmtln(file20, "{:5s}       {:5d}", "MODEL", 0)
    printfmtln(file21, "{:5s}       {:5d}", "MODEL", 0)

    for i in 1 : natom

        for j in 1 : 3

            xoc.r[i, j] = xoc.r[i, j] - other.rcm[j] + rbox2

            if (xoc.r[i, j] > input.rbox)
                xoc.r[i, j] += -input.rbox
            end

            if (xoc.r[i, j] < 0.0)
                xoc.r[i, j] += input.rbox
            end

        end

        printfmtln(file_input, "{:4s}  {:5d}  {:<3s} {:3s} {:1s} {:3d}    {:8.3f}{:8.3f}{:8.3f}{:8.3f} {:2d}",
                   "ATOM", i, pdb.atom[i], pdb.res[i], other.cad[i], other.ind1[i], xoc.r[i, 1], xoc.r[i, 2],
                   xoc.r[i, 3], parmsolv.fcont[i], parmsolv.icont[i])

        printfmtln(file20, "{:4s}  {:5d}  {:<3s} {:3s} {:1s} {:3d}    {:8.3f}{:8.3f}{:8.3f}",
                   "ATOM", i, pdb.atom[i], pdb.res[i], other.cad[i], other.ind1[i], xoc.r[i, 1], xoc.r[i, 2],
                   xoc.r[i, 3])

        if (pdb.atom[i] == "CA")
            printfmtln(file21, "{:4s}  {:5d}  {:<3s} {:3s} {:1s} {:3d}    {:8.3f}{:8.3f}{:8.3f}",
                       "ATOM", i, pdb.atom[i], pdb.res[i], other.cad[i], other.ind1[i], xoc.r[i, 1], xoc.r[i, 2],
                       xoc.r[i, 3])
        end
    end

    close(file_input)

    printfmtln(file20, "{:s}", "ENDMDL")
    printfmtln(file21, "{:s}", "ENDMDL")

    for i in 1 : natom, j in 1 : 3
        other.rant[i, j] = xoc.r[i, j]
    end

    println("natom=", natom, " nres=", nres, " tsnap=", input.tsnap, " nbloc=", input.nbloc)
    println("epot=", epot0, " ekin=", ekin0, " etot=", etot0, " epothb=", epothb0, " nhb=", nhb)

    file8 = open(input.file8, "w")
    #TODO: write INPUT and DISTANCIES to file8
    println(file8, "natom=", natom, " nres=", nres, " tsnap=", input.tsnap, " nbloc=", input.nbloc)
    println(file8, "epot=", epot0, " ekin=", ekin0, " etot=", etot0, " epothb=", epothb0, " nhb=", nhb)

    file12 = open(input.file12, "w")

    temps = 0.0
    temps0 = 0.0

    println(file12, "#Energia inicial ", epot0, " ", epotmol0, " ", epothb0, " ", epothbmol0, " ", ekin0, " ", etot0, " ", nhb)

    epotmol = 0.0
    epothbmol = 0.0
    
    # Main loop were iterations occur
    for ibloc in 1 : input.nbloc

        tacum = 0.0

        epot_i = 0.0
        etot_i = 0.0

        while (tacum < input.tsnap)
            
            for i in 1 : (natom - 1), j in (i + 1) : natom
                if (intr.istruct[i, j] == 0)
                    intr.inter[i, j] = 0
                end
            end

            if (input.isec == 0)

                for i in 1 : (natom - 1), j in (i + 1) : natom
                    intr.inter[i, j] = 0
                end

                potencial(natom, xoc, pous, intr, cov, pdb, fisic, param, parmsolv)

                nhb = 0

                for i in 1 : natom
                    atpres.ihb[i] = 0
                end

                for i in 1 : (nres - 4)

                    ii = atpres.io[i]

                    for j in (i + 4) : nres

                        if (pdb.res[atpres.ica[j]] == "PRO")
                            continue
                        end

                        jj = atpres.ih[j]

                        if (atpres.ihb[ii] + atpres.ihb[jj] != 0)
                            continue
                        end

                        ic = 0

                        n1 = atpres.io[i]
                        n2 = atpres.ih[j]

                        rij1 = dbox(n2, n1, 1, xoc)
                        rij2 = dbox(n2, n1, 2, xoc)
                        rij3 = dbox(n2, n1, 3, xoc)

                        roh = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

                        if (roh < distancies.rohmax && roh > distancies.rohmin)
                            ic += 1
                        elseif (intr.istruct[n1, n2] == 1)
                            if (other.ireg[n1, n2] <= intr.nstep[n1, n2] && other.ireg[n1, n2] > 1)
                                ic += 1
                            end
                        end

                        n1 = atpres.io[i]
                        n2 = atpres.in[j]

                        rij1 = dbox(n2, n1, 1, xoc)
                        rij2 = dbox(n2, n1, 2, xoc)
                        rij3 = dbox(n2, n1, 3, xoc)

                        rno = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

                        if (rno < distancies.rnomax && rno > distancies.rnomin)
                            ic += 1
                        elseif (intr.istruct[n1, n2] == 1)
                            if (other.ireg[n1, n2] <= intr.nstep[n1, n2] && other.ireg[n1, n2] > 1)
                                ic += 1
                            end
                        end

                        n1 = atpres.ico[i]
                        n2 = atpres.ih[j]

                        rij1 = dbox(n2, n1, 1, xoc)
                        rij2 = dbox(n2, n1, 2, xoc)
                        rij3 = dbox(n2, n1, 3, xoc)

                        rch = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

                        if (rch < distancies.rchmax && rch > distancies.rchmin)
                            ic += 1
                        elseif (intr.istruct[n1, n2] == 1)
                            if (other.ireg[n1, n2] <= intr.nstep[n1, n2] && other.ireg[n1, n2] > 1)
                                ic += 1
                            end
                        end

                        if (ic == 3)

                            nhb += 1

                            n1 = atpres.io[i]
                            n2 = atpres.ih[j]

                            atpres.ihb[n1] = 1
                            atpres.ihb[n2] = 1

                            sto = parmsolv.fcont[n1] * parmsolv.fcont[n2]
                            ene = input.ehb * sto + input.ehbc * (1.0 - sto)

                            #println(nhb, " ", sto, " ", ene)
                            #println(input.ehb, " ", input.ehbc)
                            #readline()

                            creapouhb(n1, n2, distancies.rohmin, distancies.roha, distancies.rohb, distancies.rohmax, ene, pous, intr)

                            other.ireg[n1, n2] = 2

                            n1 = atpres.io[i]
                            n2 = atpres.in[j]

                            sto = parmsolv.fcont[n1] * parmsolv.fcont[n2]
                            ene = input.ehb * sto + input.ehbc * (1.0 - sto)

                            creapouhb(n1, n2, distancies.rnomin, distancies.rnoa, distancies.rnob, distancies.rnomax, ene, pous, intr)

                            other.ireg[n1, n2] = 2

                            n1 = atpres.ico[i]
                            n2 = atpres.ih[j]

                            sto = parmsolv.fcont[n1] * parmsolv.fcont[n2]
                            ene = input.ehb * sto + input.ehbc * (1.0 - sto)

                            creapouhb(n1, n2, distancies.rchmin, distancies.rcha, distancies.rchb, distancies.rchmax, ene, pous, intr)

                            other.ireg[n1, n2] = 2

                        else

                            n1 = atpres.io[i]
                            n2 = atpres.ih[j]

                            intr.istruct[n1, n2] = 0
                            other.ireg[n1, n2] = 0

                            n1 = atpres.io[i]
                            n2 = atpres.in[j]

                            intr.istruct[n1, n2] = 0
                            other.ireg[n1, n2] = 0

                            n1 = atpres.ico[i]
                            n2 = atpres.ih[j]

                            intr.istruct[n1, n2] = 0
                            other.ireg[n1, n2] = 0

                        end

                    end

                end

                for i in 1 : (nres - 4)

                    if (pdb.res[atpres.ica[i]] == "PRO")
                        continue
                    end

                    ii = atpres.ih[i]

                    for j in (i + 4) : nres

                        jj = atpres.io[j]

                        if (atpres.ihb[ii] + atpres.ihb[jj] != 0)
                            continue
                        end

                        ic = 0

                        n1 = atpres.ih[i]
                        n2 = atpres.io[j]

                        rij1 = dbox(n2, n1, 1, xoc)
                        rij2 = dbox(n2, n1, 2, xoc)
                        rij3 = dbox(n2, n1, 3, xoc)

                        roh = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

                        if (roh < distancies.rohmax && roh > distancies.rohmin)
                            ic += 1
                        elseif (intr.istruct[n1, n2] == 1)
                            if (other.ireg[n1, n2] <= intr.nstep[n1, n2] && other.ireg[n1, n2] > 1)
                                ic += 1
                            end
                        end

                        n1 = atpres.in[i]
                        n2 = atpres.io[j]

                        rij1 = dbox(n2, n1, 1, xoc)
                        rij2 = dbox(n2, n1, 2, xoc)
                        rij3 = dbox(n2, n1, 3, xoc)

                        rno = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

                        if (rno < distancies.rnomax && rno > distancies.rnomin)
                            ic += 1
                        elseif (intr.istruct[n1, n2] == 1)
                            if (other.ireg[n1, n2] <= intr.nstep[n1, n2] && other.ireg[n1, n2] > 1)
                                ic += 1
                            end
                        end

                        n1 = atpres.ih[i]
                        n2 = atpres.ico[j]

                        rij1 = dbox(n2, n1, 1, xoc)
                        rij2 = dbox(n2, n1, 2, xoc)
                        rij3 = dbox(n2, n1, 3, xoc)

                        rch = sqrt(rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0)

                        if (rch < distancies.rchmax && rch > distancies.rchmin)
                            ic += 1
                        elseif (intr.istruct[n1, n2] == 1)
                            if (other.ireg[n1, n2] <= intr.nstep[n1, n2] && other.ireg[n1, n2] > 1)
                                ic += 1
                            end
                        end

                        if (ic == 3)

                            nhb += 1

                            n1 = atpres.ih[i]
                            n2 = atpres.io[j]

                            atpres.ihb[n1] = 1
                            atpres.ihb[n2] = 1

                            sto = parmsolv.fcont[n1] * parmsolv.fcont[n2]
                            ene = input.ehb * sto + input.ehbc * (1.0 - sto)

                            creapouhb(n1, n2, distancies.rohmin, distancies.roha, distancies.rohb, distancies.rohmax, ene, pous, intr)

                            other.ireg[n1, n2] = 2

                            n1 = atpres.in[i]
                            n2 = atpres.io[j]

                            sto = parmsolv.fcont[n1] * parmsolv.fcont[n2]
                            ene = input.ehb * sto + input.ehbc * (1.0 - sto)

                            creapouhb(n1, n2, distancies.rnomin, distancies.rnoa, distancies.rnob, distancies.rnomax, ene, pous, intr)

                            other.ireg[n1, n2] = 2

                            n1 = atpres.ih[i]
                            n2 = atpres.ico[j]

                            sto = parmsolv.fcont[n1] * parmsolv.fcont[n2]
                            ene = input.ehb * sto + input.ehbc * (1.0 - sto)

                            creapouhb(n1, n2, distancies.rchmin, distancies.rcha, distancies.rchb, distancies.rchmax, ene, pous, intr)

                            other.ireg[n1, n2] = 2

                        else

                            n1 = atpres.ih[i]
                            n2 = atpres.io[j]

                            intr.istruct[n1, n2] = 0
                            other.ireg[n1, n2] = 0

                            n1 = atpres.in[i]
                            n2 = atpres.io[j]

                            intr.istruct[n1, n2] = 0
                            other.ireg[n1, n2] = 0

                            n1 = atpres.ih[i]
                            n2 = atpres.ico[j]

                            intr.istruct[n1, n2] = 0
                            other.ireg[n1, n2] = 0

                        end

                    end

                end

            end

            potencial(natom, xoc, pous, intr, cov, pdb, fisic, param, parmsolv)

            for i in 1 : (natom - 1), j in (i + 1) : natom

                if (other.ireg[i, j] != 0 || cov.icov[i, j] != 0)
                    continue
                end

                rmin2 = pous.rstep[i, j, 1] * pous.rstep[i, j, 1]
                rmax2 = pous.rstep[i, j, 2] * pous.rstep[i, j, 2]

                rij1 = dbox(i, j, 1, xoc)
                rij2 = dbox(i, j, 2, xoc)
                rij3 = dbox(i, j, 3, xoc)

                rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

                if (rmod2 < rmin2)
                    other.ireg[i, j] = 1
                elseif (rmod2 > rmax2)
                    other.ireg[i, j] = 3
                else
                    other.ireg[i, j] = 2
                end

            end

            # llista de solapaments plausibles
            for i in 1 : (natom - 1)

                shake.ishk[i] = 0
                npt.ipot[i] = 0

                for j in (i + 1) : natom

                    if (cov.icov[i, j] != 0)
                        continue
                    end

                    rij1 = dbox(i, j, 1, xoc)
                    rij2 = dbox(i, j, 2, xoc)
                    rij3 = dbox(i, j, 3, xoc)

                    rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

                    #TODO: check if this is correct
                    if (rmod2 < rshake2)
                        shake.ishk[i] += 1
                        k = shake.ishk[i]
                        shake.nshk[i, k] = j
                    end

                    if (intr.istruct[i, j] != 1 && rmod2 < rpot2)
                        npt.ipot[i] += 1
                        k = npt.ipot[i]
                        npt.npot[i, k] = j
                    end

                end

            end

            mem1 = 0
            mem2 = 0

            ierr2 = 0

            temps0 = temps

            for i in 1 : natom
                for j in 1 : 3
                    xoc.r[i, j] = other.rant[i, j]
                end
                #println("rant: ", other.rant[i, 1], " ", other.rant[i, 2], " ", other.rant[i, 3])
            end

            #=if (ibloc == 2)
                readline()
            end=#

            tacene = 0.0

            while (tacene < input.tene)

                dmdshake(natom, npair, mem1, mem2, xoc, cov, shake)

                for i in 1 : (natom - 1), j in (i + 1) : natom
                    if (intr.istruct[i, j] == 0)
                        intr.inter[i, j] = 0
                    end
                end

                if (input.isolv != 0)
                    enchufa(natom, input.dcut, xoc, intr, cov, atpres, pdb, npt)
                end

                icont = 0

                for i in 1 : (natom - 1)
                    for j in (i + 1) : natom
                        if (intr.inter[i, j] != 1)
                            continue
                        end

                        rij1 = dbox(j, i, 1, xoc)
                        rij2 = dbox(j, i, 2, xoc)
                        rij3 = dbox(j, i, 3, xoc)

                        rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

                        vij1 = xoc.v[i, 1] - xoc.v[j, 1]
                        vij2 = xoc.v[i, 2] - xoc.v[j, 2]
                        vij3 = xoc.v[i, 3] - xoc.v[j, 3]

                        prod = rij1 * vij1 + rij2 * vij2 + rij3 * vij3

                        rij = sqrt(rmod2)

                        k = 1

                        while (k <= intr.nstep[i, j] && rij > pous.rstep[i, j, k])
                            k += 1
                        end

                        k0 = other.ireg[i, j]

                        #print(k,"->",k0," ")
                        #printfmt("{:2.8f} ", prod*1e-13)

                        # Dirty and lazy trick to pass by reference
                        ich = Ref{Int}(0)

                        if (k0 > k && prod < 0.0)

                            icont += 1
                            dpot = -pous.estep[i, j, k0 - 1] * constants.facte

                            #=println(rij)
                            println(k0, " ", pous.estep[i, j, k0 - 1])
                            println(i, " ", j, " ", dpot)
                            readline()=#

                            chgmomene(i, j, rij1, rij2, rij3, dpot, ich, xoc)
                            ich = ich[]

                            #print(ich, " ")

                            if (ich == 1)
                                other.ireg[i, j] -= 1
                            end

                            ierr2 += 1

                        elseif (k0 < k && prod > 0.0)

                            icont += 1
                            dpot = pous.estep[i, j, k0] * constants.facte

                            chgmomene(i, j, rij1, rij2, rij3, dpot, ich, xoc)
                            ich = ich[]

                            #print(ich, " ")

                            if (ich == 1)
                                other.ireg[i, j] += 1
                            end

                            ierr2 += 1

                        end
                    end
                    #println()
                end
                #readline()

                #=if (input.tact * 1 - 1e-15 < tacene)
                    for i in 1 : natom
                        for j in 1 : 3
                            printfmt("{:3.10f} ", xoc.v[i, j] / 1e12)
                        end
                        println()
                    end
                    readline()
                end=#

                #println(icont)

                for j in 1 : 3, i in 1 : natom
                    xoc.r[i, j] = xoc.r[i, j] + input.tact * xoc.v[i, j]
                end

                if (constants.icm == 0)
                    for i in 1 : natom, j in 1 : 3
                        if (xoc.r[i, j] > input.rbox)
                            xoc.r[i, j] -= input.rbox
                        end

                        if (xoc.r[i, j] < 0.0)
                            xoc.r[i, j] += input.rbox
                        end
                    end
                end

                tacum += input.tact
                tacene += input.tact
                temps += input.tact

                #TODO: delete this to return randomness to the program
                input.iterm = 0

                if (input.iterm == 1)

                    fi = getUniformRandom()

                    i = trunc(Int, fi) + 1
                    if (i == natom)
                        i -= 1
                    end

                    #TODO: delete the non-random function above and call this function
                    #i = getAtomInRange(natom)

                    for j in 1 : 3
                        rnd_gauss(fi, xoc.xm[i], input.temp)
                        xoc.v[i, j] = fi / constants.a
                    end

                end

                #=sum = 0.0
                for i in 1 : natom, j in 1 : 3
                    sum += xoc.v[i, j]
                end
                printfmtln("{:2.8f}", sum * 1e-13)
                #readline()=#

            end

            ekin = 0.0
            temp_delete = 0.0

            for j in 1 : 3, i in 1 : natom
                ekin = ekin + 0.5 * xoc.xm[i] * (xoc.v[i, j] * constants.a) ^ 2.0
            end

            ekin2 = ekin / constants.facte

            println("ekin: ", ekin)

            #=for i in 1 : natom
                for j in 1 : 3
                    printfmt("{:3.10f} ", xoc.v[i, j] / 1e12)
                end
                println()
            end
            readline()=#

            if (constants.icm == 1)

                for j in 1 : 3
                    other.rcm[j] = 0.0

                    # TODO: problem somehow lies here, with the xoc.r part
                    for i in 1 : natom
                        #println(other.rcm[j] + xoc.xm[i], " ", xoc.r[i, j])
                        other.rcm[j] = other.rcm[j] + xoc.xm[i] * xoc.r[i, j]
                    end

                    #=if (ibloc == 2)
                        readline()
                    end=#

                    other.rcm[j] /= xmassa
                end

                for i in 1 : natom
                    #println("rant before: ", other.rant[i, 1], " ", other.rant[i, 2], " ", other.rant[i, 3])
                    #println("xoc.r: ", xoc.r[i, 1], " ", xoc.r[i, 2], " ", xoc.r[i, 3])
                    #println("rbox2: ", rbox2)
                    #println("other.rcm: ", other.rcm[1], " ", other.rcm[2], " ", other.rcm[3])
                    for j in 1 : 3
                        other.rant[i, j] = xoc.r[i, j] - other.rcm[j] + rbox2
                    end
                    #println("rant after: ", other.rant[i, 1], " ", other.rant[i, 2], " ", other.rant[i, 3])
                    #println()
                end
                #readline()

                #println(other.rcm[1], " ", other.rcm[2], " ", other.rcm[3])

            else

                for i in 1 : natom, j in 1 : 3

                    other.rant[i, j] = xoc.r[i, j]

                end

            end

            epothb = 0.0
            epot = 0.0
            dist_sum = 0.0

            for i in 1 : (natom - 1)
                for j in (i + 1) : natom

                    if (intr.inter[i, j] != 1)
                        continue
                    end

                    rij1 = dbox(j, i, 1, xoc)
                    rij2 = dbox(j, i, 2, xoc)
                    rij3 = dbox(j, i, 3, xoc)

                    rmod2 = rij1 ^ 2.0 + rij2 ^ 2.0 + rij3 ^ 2.0

                    dist = sqrt(rmod2)
                    dist_sum += dist

                    k = intr.nstep[i, j]

                    while (k > 0 && dist < pous.rstep[i, j, k])

                        epot += -pous.estep[i, j, k]

                        if (pdb.imol[i] != pdb.imol[j])
                            epotmol += -pous.estep[i, j, k]
                        end

                        if (intr.istruct[i, j] == 1)
                            epothb += -pous.estep[i, j, k]

                            if (pdb.imol[i] != pdb.imol[j])
                                epothbmol += -pous.estep[i, j, k]
                            end
                        end

                        k -= 1

                    end
                end
            end

            sum = 0.0
            for i in 1 : natom
                for j in 1 : natom
                    sum += pous.estep[i, j, 1]
                end
            end
            #println("estep: ", sum)
            #readline()

            sum2 = 0
            for i in 1 : natom
                for j in 1 : natom
                    sum2 += pous.rstep[i, j, 1]
                end
            end
            #println("rstep: ", sum2)

            sum3 = 0
            for i in 1 : natom
                for j in 1 : 3
                    sum3 += xoc.r[i, j] ^ 2.0
                end
            end
            #println("xoc.r: ", sum3)

            sum4 = 0
            for i in 1 : natom
                for j in 1 : 3
                    sum4 += xoc.v[i, j] ^ 2.0
                end
            end
            #println("xoc.v: ", sum4)

            sum5 = 0
            for i in 1 : natom
                sum5 += xoc.xm[i]
            end
            #println("xoc.xm: ", sum5)

            #println("dist: ", dist_sum)

            #println(epot, " ", epotmol, " ", epothb, " ", epothbmol)

            etot = epot + ekin2

            epot_i = epot
            etot_i = etot

            if (input.iprint == 1)
                println(file12, temps, " ", epot, " ", epotmol, " ", epothb, " ", epothbmol, " ", ekin2, " ", etot, " ", nhb)
            end

        end

        printfmtln(file20, "{:5s}       {:5d}", "MODEL", ibloc)
        printfmtln(file21, "{:5s}       {:5d}", "MODEL", ibloc)

        for i in 1 : natom

            printfmtln(file20, "{:4s}  {:5d}  {:<3s} {:3s} {:1s} {:3d}    {:8.3f}{:8.3f}{:8.3f}",
                       "ATOM", i, pdb.atom[i], pdb.res[i], other.cad[i], other.ind1[i], xoc.r[i, 1], xoc.r[i, 2],
                       xoc.r[i, 3])

            if (pdb.atom[i] == "CA")
                printfmtln(file21, "{:4s}  {:5d}  {:<3s} {:3s} {:1s} {:3d}    {:8.3f}{:8.3f}{:8.3f}",
                           "ATOM", i, pdb.atom[i], pdb.res[i], other.cad[i], other.ind1[i], xoc.r[i, 1], xoc.r[i, 2],
                           xoc.r[i, 3])
            end

        end

        printfmtln(file20, "{:s}", "ENDMDL")
        printfmtln(file21, "{:s}", "ENDMDL")

        println("Temps ", temps, " hbonds ", nhb, " epot ", epot_i, " ", etot_i)
        println(file8, "Temps ", temps, " hbonds ", nhb, " epot ", epot_i, " ", etot_i)

    end

    close(file8)
    close(file12)
    close(file20)
    close(file21)

    ekin = ekin / constants.facte

    file19 = open(input.file19, "w")

    for i in 1 : natom
        printfmtln(file19, "{:4s}  {:5d}  {:<3s} {:3s} {:1s} {:3d}    {:8.3f}{:8.3f}{:8.3f}",
                   "ATOM", i, pdb.atom[i], pdb.res[i], other.cad[i], other.ind1[i], xoc.r[i, 1], xoc.r[i, 2],
                   xoc.r[i, 3])
    end

    close(file19)

end

main()