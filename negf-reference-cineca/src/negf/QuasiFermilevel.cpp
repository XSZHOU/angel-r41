#include "QuasiFermilevel.h"

QuasiFermilevel::QuasiFermilevel(const Geometry * xspace_,
                    const Kspace * kspace_,
                    const Energies * energies_,
                    const GreenFunctions * gf_) throw (Exception *):
     xspace(xspace_),
     kspace(kspace_),
     energies(energies_),
     gf(gf_)
{
    NEGF_ASSERT(xspace!=0 && kspace!=0 && energies!=0 && gf!=0, "null pointer encountered.");
    this->Nx   = xspace->get_num_internal_vertices();
    this->NxNn = Nx*Nn;
    this->Nk   = kspace->get_number_of_points();
    this->NE   = energies->get_number_of_points();
    this->myNE = energies->get_my_number_of_points();

    this->options = gf->get_options();

    if (mpi->get_rank()==constants::mpi_master_rank) {
        electron_qfl.resize(xspace->get_num_vertices(), 0.0);
        hole_qfl    .resize(xspace->get_num_vertices(), 0.0);
    }
}


QuasiFermilevel::~QuasiFermilevel() {}


void QuasiFermilevel::set_densities(const vector<double> & edens_, const vector<double> & hdens_)
{STACK_TRACE(
    NEGF_FASSERT(edens_.size()==xspace->get_num_internal_vertices()
             && hdens_.size()==xspace->get_num_internal_vertices(),
             "inconsistent array sizes: expected %d, got %d/%d.", xspace->get_num_vertices(), edens_.size(), hdens_.size());

    this->edens = edens_;
    this->hdens = hdens_;
);}


/** for each point we seek the quantities EFn satisfying
 *  \f$
 *      \sum_k \sum_E \sum_{n \in CB} f(EF_n,E,kT) * Im GR(xi,xi)   =   n(xi)
 *      \sum_k \sum_E \sum_{n \in VB} f(EF_p,E,kT) * Im GR(xi,xi)   =   p(xi)
 *  \f$
 */
void QuasiFermilevel::calculate() throw (Exception *)
{STACK_TRACE(
    logmsg->emit_small_header("Calculating Quasi-Fermilevels");
    logmsg->emit(LOG_INFO,"Calculating DOS...");
    this->calculate_dos();
    mpi->synchronize_processes();

    // only master process does calculation
    if (!mpi->get_rank()==constants::mpi_master_rank) return;

    // --------------------------------------------------------------
    // do a little Newton loop for each x:
    // \sum_E w_E     f(EF_n,E,kT) * electron_dos[ee][xx] = n[xx]
    // \sum_E w_E (1-f(EF_p,E,kT)) * hole_dos    [ee][xx] = p[xx]
    // --------------------------------------------------------------
    logmsg->emit(LOG_INFO,"Calculating EFn, EFp...");

    vector<uint> cb_dofs; options->get_conduction_degrees_of_freedom(cb_dofs);
    vector<uint> vb_dofs; options->get_valence_degrees_of_freedom(vb_dofs);

    uint max_iterations = 200;
    uint iter;

    for (uint xx=1; xx<Nx; xx++)
    {
        uint ii = xspace->get_global_vertex_index(xx-1); // xx is 1-based, ii is 0-based

        // electrons
        if (cb_dofs.size()>0)
        {
            double EF_n = electron_qfl[ii];
            for (iter = 1; iter <= max_iterations; iter++)
            {
                // find F(mu)
                double F = this->calculate_density(EF_n, xx-1, true, false) - edens[xx-1];

                // find dF_dmu(mu)
                double dF_dmu = this->calculate_density(EF_n, xx-1, true, true);
                logmsg->emit(LOG_INFO_L3," Iteration %d electrons: F=%g, dF_dmu = %g", iter, F, dF_dmu);

                // find update
                double update = -1.0/dF_dmu * F;

                // limit update!
                double max_update = constants::convert_from_SI(units::energy, constants::max_kpquasi_change_V * constants::SIec);
                if (fabs(update) > max_update) {
                    update = negf_math::sign(update) * max_update;
                    logmsg->emit(LOG_INFO_L3," Iteration %d electrons: Limiting to %e", iter, update);
                }

                // peform update
                EF_n += update;

                // convergence check
                if (fabs(update) < constants::convert_from_SI(units::energy, 1e-8*constants::SIec)) {
                    logmsg->emit(LOG_INFO_L2,"xx=%d: Quasi-Fermilevel convergence reached after %d iterations: EF_n=%.6g",xx,iter,EF_n);
                    break;
                }
            }
            electron_qfl[ii] = EF_n;
        } else {
            electron_qfl[ii] = -100.0;
        }

        // holes
        if (vb_dofs.size()>0)
        {
            double EF_p = hole_qfl[ii];
            for (iter = 1; iter <= max_iterations; iter++)
            {
                // find F(mu)
                double F = this->calculate_density(EF_p, xx-1, false, false) - hdens[xx-1];

                // find dF_dmu(mu)
                double dF_dmu = this->calculate_density(EF_p, xx-1, false, true);
                logmsg->emit(LOG_INFO_L3," Iteration %d holes: F=%g, dF_dmu = %g", iter, F, dF_dmu);

                // find update
                double update = -1.0/dF_dmu * F;

                // limit update!
                double max_update = constants::convert_from_SI(units::energy, constants::max_kpquasi_change_V * constants::SIec);
                if (fabs(update) > max_update) {
                    update = negf_math::sign(update) * max_update;
                    logmsg->emit(LOG_INFO_L3," Iteration %d holes: Limiting to %e", iter, update);
                }

                // peform update
                EF_p += update;

                // convergence check
                if (fabs(update) < constants::convert_from_SI(units::energy, 1e-8*constants::SIec)) {
                    logmsg->emit(LOG_INFO_L2,"xx=%d: Quasi-Fermilevel convergence reached after %d iterations: EF_p=%.6g",xx,iter,EF_p);
                    break;
                }
            }
            hole_qfl[ii] = EF_p;
        } else {
            hole_qfl[ii] = +100.0;
        }
    }
);}

// xx should be 0-based!
double QuasiFermilevel::calculate_density(double EF, uint xx, bool e_or_h, bool derivative)
{STACK_TRACE(
    NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank, "only master may call this function.");

    double kpmethod = options->get("kp_method");
    double kT   = constants::convert_from_SI(units::energy, constants::SIkb * options->get("temperature"));
    double spin = (e_or_h)
            ? get_spin_degeneracy(kpmethod, quantities::electron_density)
            : get_spin_degeneracy(kpmethod, quantities::hole_density);

    double result = 0.0;
    for (uint ee=0; ee<NE; ee++)
    {
        double E  = energies->get_energy_from_global_idx(ee);
        double dE = energies->get_weight_from_global_idx(ee);
        double nu = (E - EF) / kT;

        if (e_or_h==true)
        {
            // electrons
            if (derivative) {
                double f = 1.0 / (2.0 + negf_math::exp(nu) + negf_math::exp(-nu)) * 1.0/kT;
                result += spin*f * electron_dos[ee][xx] * dE;
            } else {
                double f = 1.0 / (1.0 + negf_math::exp(nu));
                result += spin*f * electron_dos[ee][xx] * dE;
            }
        } else
        {
            // holes
            if (derivative) {
                double f = 1.0 / (2.0 + negf_math::exp(nu) + negf_math::exp(-nu)) * 1.0/kT;
                result += -spin*f * hole_dos[ee][xx] * dE;
            } else {
                double f = 1.0 / (1.0 + negf_math::exp(nu));
                result +=  spin*(1.0-f) * hole_dos[ee][xx] * dE;
            }
        }
    }

    return result;
);}


/** calculate the LDOS from GR.
 *  result: NE*Nx matrices electron_dos and hole_dos in mater thread.
 */
void QuasiFermilevel::calculate_dos()
{STACK_TRACE(
    double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
    vector<uint> cb_dofs; options->get_conduction_degrees_of_freedom(cb_dofs);

    // -----------------------------------------------
    // calculate own LDOS
    // -----------------------------------------------
    logmsg->emit(LOG_INFO,"Calculating own LDOS...");
    vector<double> tmp(Nx, 0.0);
    vector< vector<double> > my_dos_n; my_dos_n.resize(myNE, tmp);
    vector< vector<double> > my_dos_p; my_dos_p.resize(myNE, tmp);

    for (uint ee2 = 0; ee2 < myNE; ee2++)
    {
        uint ee = energies->get_global_index(ee2);

        for (uint kk = 0; kk < Nk; kk++)
        {
            // get retarded GF
            Matc & GR = gf->get_retarded(kk,ee);

            // -----------------------------------------------
            // add diagonal to the density of states (GR-GA)
            // copy-paste from PostProcessing.cpp
            // -----------------------------------------------
            double wk = kspace->get_point(kk).get_weight() * kspace_factor;
            for (uint nn=1; nn<=Nn; nn++)
            {
                bool electron_band = false;
                for (uint jj=0; jj<cb_dofs.size(); jj++) {
                    if(cb_dofs[jj]+1==nn) {
                        electron_band = true;
                        break;
                    }
                }
                double spin = (electron_band)
                                ? get_spin_degeneracy(options->get("kp_method"), quantities::electron_density)
                                : get_spin_degeneracy(options->get("kp_method"), quantities::hole_density);

                for (uint xx=1; xx<=Nx; xx++) {
                    uint idx = get_mat_idx(xx,nn,Nx); // (xx-1)*Nn+nn;
                    // LDOS = i/2pi * (GR-GA)
                    if (electron_band) {
                        my_dos_n[ee2][xx-1] += wk * (-1.0/constants::pi) * spin * GR(idx,idx).imag();
                    } else {
                        my_dos_p[ee2][xx-1] += wk * (-1.0/constants::pi) * spin * GR(idx,idx).imag();
                    }
                }
            }
        }

        // correct with geometrical factor in the case of dumb orthogonal-basis-hamiltonian
        if (constants::old_orthogonal
            && (fabs(options->get("kp_method") - 0.0) < 1e-14 || fabs(options->get("kp_method") - 3.0) < 1e-14))
        {
            for (uint xx=1; xx<=Nx; xx++) {
                double x_cell_length = 0.0;
                Vertex * v = xspace->get_vertex(xspace->get_global_vertex_index(xx-1));
                for (uint ii=0; ii < xspace->get_edges_near(v).size(); ii++) {
                    x_cell_length += 0.5 * xspace->get_edges_near(v)[ii]->get_length();
                }
                my_dos_n[ee2][xx-1] = my_dos_n[ee2][xx-1] / x_cell_length;
                my_dos_p[ee2][xx-1] = my_dos_p[ee2][xx-1] / x_cell_length;
            }
        }
    }
    mpi->synchronize_processes();

    // -----------------------------------------------------------------------
    // communicate own DOS to master thread, which puts everything together
    // -----------------------------------------------------------------------
    logmsg->emit(LOG_INFO,"Aggregating total DOS...");
    if (mpi->get_rank()==constants::mpi_master_rank)
    {
        tmp.assign(Nx, -1.0);
        this->electron_dos.resize(NE, tmp);
        this->hole_dos.resize(NE, tmp);

        for (int pp=0; pp<mpi->get_num_procs(); pp++)
        {
            if (pp==mpi->get_rank()) {
                // own contribution
                for (uint ee2=0; ee2 < myNE; ee2++) {
                    uint ee = energies->get_global_index(ee2);
                    for (uint xx=1; xx<=Nx; xx++) {
                        this->electron_dos[ee][xx-1] = my_dos_n[ee2][xx-1];
                        this->hole_dos    [ee][xx-1] = my_dos_p[ee2][xx-1];
                    }
                }
            } else {
                uint num_energy_points_pp = energies->get_number_of_points(pp);
                vector<double> ndos_pp; ndos_pp.resize(Nx, 0.0);
                vector<double> pdos_pp; pdos_pp.resize(Nx, 0.0);
                for (uint ee=energies->get_start_global_idx(pp); ee<=energies->get_stop_global_idx(pp); ee++)
                {
                    NEGF_ASSERT(num_energy_points_pp==energies->get_stop_global_idx(pp)-energies->get_start_global_idx(pp)+1, "inconsistent number of points.");

                    int tag = 987;
                    int source = pp;
                    uint Nx2 = Nx;
                    mpi->recv(ndos_pp, Nx2, source, tag);
                    mpi->recv(pdos_pp, Nx2, source, tag);
                    NEGF_ASSERT(ndos_pp.size()==Nx && pdos_pp.size()==Nx, "received DOS matrices did not have right size.");

                    NEGF_ASSERT(electron_dos[ee].size()==Nx && hole_dos[ee].size()==Nx, "total DOS matrices did not have right size.");
                    for (uint xx=1; xx<=Nx; xx++) {
                        this->electron_dos[ee][xx-1] = ndos_pp[xx-1];
                        this->hole_dos    [ee][xx-1] = pdos_pp[xx-1];
                    }
                }
            }
        }
        // check if every point was received
        NEGF_ASSERT(electron_dos.size()==hole_dos.size(), "something went wrong.");
        for (uint ee=0; ee < electron_dos.size(); ee++) {
            NEGF_ASSERT(electron_dos[ee].size()==hole_dos[ee].size(), "something went wrong.");
            for (uint ii=0; ii<electron_dos[ee].size(); ii++) {
                NEGF_FASSERT(electron_dos[ee][ii] != -1.0, "did not assign electron_dos[%d][%d]",ee,ii);
                NEGF_FASSERT(hole_dos    [ee][ii] != -1.0, "did not assign hole_dos[%d][%d]",ee,ii);
            }
        }
    } else {
        int tag = 987;
        int dest = constants::mpi_master_rank;
        for (uint ee2=0; ee2<myNE; ee2++) {
            NEGF_ASSERT(my_dos_n[ee2].size()==Nx && my_dos_p[ee2].size()==Nx, "own DOS matrices did not have right size.");
            mpi->send(my_dos_n[ee2], dest, tag);
            mpi->send(my_dos_p[ee2], dest, tag);
        }
    }
    logmsg->emit_all(LOG_INFO_L3,"p%d is finished calculating the DOS.",mpi->get_rank());
    mpi->synchronize_processes();
);}


vector<double> QuasiFermilevel::get_electron_qfl() const throw (Exception *)
{STACK_TRACE(
    NEGF_ASSERT(this->electron_qfl.size()>0, "electron QFL was not yet computed.");
    return this->electron_qfl;
);}


vector<double> QuasiFermilevel::get_hole_qfl() const throw (Exception *)
{STACK_TRACE(
    NEGF_ASSERT(this->hole_qfl.size()>0, "electron QFL was not yet computed.");
    return this->hole_qfl;
);}


void QuasiFermilevel::write_to_file(char * filename) const throw (Exception *)
{STACK_TRACE(
    NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank, "transmission spectrum is stored in master process only!");
    Matd output(2, xspace->get_num_vertices());
    for (uint xx=0; xx<xspace->get_num_vertices(); xx++) {
        output(1, xx+1) = electron_qfl[xx];
        output(2, xx+1) = hole_qfl[xx];
    }
    negf::write_matrix(filename, output);
);}
