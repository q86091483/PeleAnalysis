#include <string>
#include <iostream>
#include <algorithm>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_BLFort.H>
#include <mechanism.H>
#include <PelePhysics.H>

using namespace amrex;

extern "C" {
    void pushvtog(const int* lo,  const int* hi,
                const int* dlo, const int* dhi,
                Real* U, const int* Ulo, const int* Uhi,
                const int* nc);
    void gradient(const int* lo,  const int* hi,
                const Real* U, const int* Ulo, const int* Uhi,
                Real* G, const int* Glo, const int* Ghi,
                const Real* dx);
}

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<plotfilename> \n\tOptions:\n\tis_per=<L M N> gradVar=<name>\n";
    exit(1);
}

std::string basename(const std::string& filename)
{
    std::string tmp(filename);
    if (tmp[tmp.length() - 1] == '/') {
        tmp = tmp.substr(0, tmp.size()-1);
    }

    if (const char *slash = strrchr(tmp.c_str(), '/')) {
        // Got at least one slash -- return the following tail.
        return std::string(slash + 1);
    } else {
        // No leading directory portion to name.
        return tmp;
    }
}

int main (int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  if (argc < 2)
    print_usage(argc,argv);
  ParmParse pp;
  const bool isioproc = ParallelDescriptor::IOProcessor();
  const int ioproc = ParallelDescriptor::IOProcessorNumber();
  if (pp.contains("help"))
    print_usage(argc,argv);

  // Read data
  std::string infile; pp.get("infile",infile);
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(infile, fileType);
  if( ! dataServices.AmrDataOk()) {
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  }
  AmrData& amrData = dataServices.AmrDataRef();

  // Initialize transport data
  pele::physics::transport::TransportParams<
      pele::physics::PhysicsType::transport_type>
      trans_parms;
  trans_parms.allocate();

  int finestLevel = amrData.FinestLevel();
  pp.query("finestLevel",finestLevel);
  int Nlev = finestLevel + 1;

  const amrex::Vector<std::string>& plotVarNames = amrData.PlotVarNames();
  for (int i=0; i<plotVarNames.size(); ++i)
  {
    amrex::Print() << plotVarNames[i] << std::endl;
  }

  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);

  int idYin = -1; 
  int idTin = -1;
  int idRin = -1;
  int idZin = -1;
  int idUin = -1;
  int idVin = -1;
  int idWin = -1;
  const int nCompIn  = NUM_SPECIES + 6;

  const std::string spName= "Y(" + spec_names[0] + ")";
  const std::string TName = "temp";
  const std::string RName = "density";
  const std::string ZName = "mixture_fraction";
  const std::string UName = "x_velocity";
  const std::string VName = "y_velocity";
  const std::string WName = "z_velocity";
  for (int i=0; i<plotVarNames.size(); ++i)
  {
    if (plotVarNames[i] == spName) idYin = i;
    if (plotVarNames[i] == TName)  idTin = i;
    if (plotVarNames[i] == RName)  idRin = i;
    if (plotVarNames[i] == ZName)  idZin = i;
    if (plotVarNames[i] == UName)  idUin = i;
    if (plotVarNames[i] == VName)  idVin = i;
    if (plotVarNames[i] == WName)  idWin = i;
  }

  const int idDout   = 0;
  const int idUYout  = idDout   + NUM_SPECIES;
  const int idMuout  = idUYout  + NUM_SPECIES;
  const int idXiout  = idMuout  + 1;
  const int idLamout = idXiout  + 1;
  const int idZout   = idLamout + 1;
  const int nCompOut = idZout   + 1;

  Vector<std::string> outNames(nCompOut);
  Vector<std::string> inNames(nCompIn);
  Vector<int> destFillComps(nCompIn);
  const int idYlocal = 0; // Xs start here
  const int idTlocal = NUM_SPECIES;   // T tempereature starts here
  const int idRlocal = NUM_SPECIES+1; // R density starts here
  const int idZlocal = NUM_SPECIES+2; // Z mixture fraction starts here
  const int idUlocal = NUM_SPECIES+3; // Z mixture fraction starts here
  const int idVlocal = NUM_SPECIES+4; // Z mixture fraction starts here
  const int idWlocal = NUM_SPECIES+5; // Z mixture fraction starts here
  for (int i=0; i<NUM_SPECIES; ++i)
  {
    destFillComps[i] = idYlocal + i;
    inNames[i]          =  "Y(" + spec_names[i] + ")";
    outNames[i+idDout]  = "rhoD(" + spec_names[i] + ")";
    outNames[i+idUYout] = "uY(" + spec_names[i] + ")";
  }
  destFillComps[idTlocal] = idTlocal;
  destFillComps[idRlocal] = idRlocal;
  destFillComps[idZlocal] = idZlocal;
  destFillComps[idUlocal] = idUlocal;
  destFillComps[idVlocal] = idVlocal;
  destFillComps[idWlocal] = idWlocal;
  inNames[idTlocal] = TName;
  inNames[idRlocal] = RName;
  inNames[idZlocal] = ZName;
  inNames[idUlocal] = UName;
  inNames[idVlocal] = VName;
  inNames[idWlocal] = WName;
  outNames[idMuout]   = "mu";
  outNames[idXiout]   = "xi";
  outNames[idLamout]  = "lambda";
  outNames[idZout]    = "mixture_fraction";

  Vector<MultiFab*> outdata(Nlev);
  Vector<Geometry> geoms(Nlev);
  amrex::RealBox real_box({AMREX_D_DECL(amrData.ProbLo()[0],
                                        amrData.ProbLo()[1],
                                        amrData.ProbLo()[2])},
                          {AMREX_D_DECL(amrData.ProbHi()[0],
                                        amrData.ProbHi()[1],
                                        amrData.ProbHi()[2])});
  amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
  geoms[0] = amrex::Geometry((amrData.ProbDomain())[0],
                              real_box,
                              amrData.CoordSys(),
                              is_periodic);
  const int nGrow = 1;

  for (int lev=0; lev<Nlev; ++lev)
  {
    const BoxArray ba = amrData.boxArray(lev);
    const DistributionMapping dm(ba);
    outdata[lev] = new MultiFab(ba,dm,nCompOut,nGrow);
    const Vector<Real>& delta = amrData.DxLevel()[lev];
    if ( lev > 0 ) {
       geoms[lev] = amrex::refine(geoms[lev - 1], 2);
    }

    MultiFab indata(ba, dm, nCompIn, nGrow);
    MultiFab chi(ba, dm, NUM_SPECIES, nGrow);
    MultiFab grads(ba, dm, BL_SPACEDIM+1, nGrow);

    Print() << "Reading data (FillVar) for level " << lev << std::endl;
    amrData.FillVar(indata,lev,inNames,destFillComps);
    Print() << "Data has been read for level " << lev << std::endl;

    // Fix up grow cells.  Use extrap for guess
    const Box& dbox = amrData.ProbDomain()[lev];
    for (amrex::MFIter mfi(indata); mfi.isValid(); ++mfi)
    {
      FArrayBox& fab = indata[mfi];
      const Box& box = mfi.validbox();
      pushvtog(BL_TO_FORTRAN_BOX(box),
               BL_TO_FORTRAN_BOX(dbox),
               BL_TO_FORTRAN_ANYD(fab),
               &nCompIn);
    }

    // Fix up fine-fine and periodic
    indata.FillBoundary(geoms[lev].periodicity());

    // Get the transport data pointer
    auto const* ltransparm = trans_parms.device_trans_parm();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(indata, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const Box& bx = mfi.tilebox();
      Array4<Real> const& Y_a = indata.array(mfi,idYlocal);
      Array4<Real> const& T_a = indata.array(mfi,idTlocal);
      Array4<Real> const& Z_a = indata.array(mfi,idZlocal);
      Array4<Real> const& U_a = indata.array(mfi,idUlocal);
      Array4<Real> const& V_a = indata.array(mfi,idVlocal);
      Array4<Real> const& W_a = indata.array(mfi,idWlocal);
      Array4<Real> const& rho_a = indata.array(mfi,idRlocal);

      Array4<Real> const& D_a = outdata[lev]->array(mfi,idDout);
      Array4<Real> const& UY_a = outdata[lev]->array(mfi,idUYout);
      Array4<Real> const& mu_a = outdata[lev]->array(mfi,idMuout);
      Array4<Real> const& xi_a = outdata[lev]->array(mfi,idXiout);
      Array4<Real> const& lam_a = outdata[lev]->array(mfi,idLamout);
      Array4<Real> const& Zout_a = outdata[lev]->array(mfi,idZout);
      Array4<Real> const& chi_a = chi.array(mfi); // Zisen - not sure what chi fields means, so only use a temporaray variable to store it

      amrex::launch(bx, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
        auto trans = pele::physics::PhysicsType::transport();
        trans.get_transport_coeffs(
          tbx, Y_a, T_a, rho_a, D_a, chi_a, mu_a, xi_a, lam_a, ltransparm);
      });

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          Real Cpmix = 0.0;
          Real Yloc[NUM_SPECIES] = {0.0};
          for (int n=0; n<NUM_SPECIES; ++n) {
              Yloc[n] = Y_a(i,j,k,n);
          }
          Real Tlocal = T_a(i,j,k);
          CKCPBS(T_a(i,j,k), Yloc, Cpmix);

          // convection - <UY|Z>
          for (int n=0; n<NUM_SPECIES; ++n) {
              UY_a(i,j,k,n) = U_a(i,j,k) * Y_a(i,j,k,n);
              UY_a(i,j,k,n) = UY_a(i,j,k,n) + V_a(i,j,k) * Y_a(i,j,k,n);
              UY_a(i,j,k,n) = UY_a(i,j,k,n) + W_a(i,j,k) * Y_a(i,j,k,n);
          }

          Zout_a(i,j,k) = Z_a(i,j,k); // Copy mixture fraction.
          amrex::Print(ioproc) << Z_a(i,j,k) << std::endl;
      });




    } // mfi
    Print() << "Derived finished for level " << lev << std::endl;
  } // lev

  std::string outDir = "cond_ISRN"; pp.query("outputDir", outDir);
  std::string outfilename(outDir + "/" + basename(infile) + "_derived"); 
  Print() << "Writing ISRN derived data to " << outfilename << std::endl;
  bool verb = false;
  WritePlotFile(outdata,amrData,outfilename,verb,outNames);

  trans_parms.deallocate();
  amrex::Finalize();
  return 0;
}
