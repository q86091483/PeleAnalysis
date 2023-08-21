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
  int idYin = -1;
  int idTin = -1;
  int idRin = -1;
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
  const std::string spName= "Y(" + spec_names[0] + ")";
  const std::string TName = "temp";
  const std::string RName = "density";
  for (int i=0; i<plotVarNames.size(); ++i)
  {
    if (plotVarNames[i] == spName) idYin = i;
    if (plotVarNames[i] == TName)  idTin = i;
    if (plotVarNames[i] == RName)  idRin = i;
  }

  const int nCompIn  = NUM_SPECIES + 2;
  const int idDout   = 0;
  const int idMuout  = idDout + NUM_SPECIES;
  const int idXiout  = idMuout + 1;
  const int idLamout = idXiout + 1;
  const int nCompOut = idLamout + 1;

  Vector<std::string> outNames(nCompOut);
  Vector<std::string> inNames(nCompIn);
  Vector<int> destFillComps(nCompIn);
  const int idYlocal = 0; // Xs start here
  const int idTlocal = NUM_SPECIES;   // T starts here
  const int idRlocal = NUM_SPECIES+1; // R starts here
  for (int i=0; i<NUM_SPECIES; ++i)
  {
    destFillComps[i] = idYlocal + i;
    inNames[i] =  "Y(" + spec_names[i] + ")";
    outNames[i] = "rhoD(" + spec_names[i] + ")";
  }
  destFillComps[idTlocal] = idTlocal;
  destFillComps[idRlocal] = idRlocal;
  inNames[idTlocal] = TName;
  inNames[idRlocal] = RName;
  outNames[NUM_SPECIES  ] = "mu";
  outNames[NUM_SPECIES+1] = "xi";
  outNames[NUM_SPECIES+2] = "lambda";

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
  const int nGrow = 0;

  for (int lev=0; lev<Nlev; ++lev)
  {
    const BoxArray ba = amrData.boxArray(lev);
    const DistributionMapping dm(ba);
    outdata[lev] = new MultiFab(ba,dm,nCompOut,nGrow);
    if ( lev > 0 ) {
       geoms[lev] = amrex::refine(geoms[lev - 1], 2);
    }
    MultiFab indata(ba,dm,nCompIn,nGrow);

    Print() << "Reading data for level " << lev << std::endl;
    amrData.FillVar(indata,lev,inNames,destFillComps);
    Print() << "Data has been read for level " << lev << std::endl;

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
      Array4<Real> const& rho_a = indata.array(mfi,idRlocal);
      Array4<Real> const& D_a = outdata[lev]->array(mfi,idDout);
      Array4<Real> const& mu_a = outdata[lev]->array(mfi,idMuout);
      Array4<Real> const& xi_a = outdata[lev]->array(mfi,idXiout);
      Array4<Real> const& lam_a = outdata[lev]->array(mfi,idLamout);
      Array4<Real> const& chi_a = outdata[lev]->array(mfi,idLamout+1);

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
          amrex::Print(ioproc) << D_a(i,j,k,5) << std::endl;
      });


    } // mfi
    Print() << "Derived finished for level " << lev << std::endl;
  } // lev

  std::string outfilename(basename(infile) + "_ISRNderived"); 
  Print() << "Writing ISRN derived data to " << outfilename << std::endl;
  bool verb = false;
  WritePlotFile(outdata,amrData,outfilename,verb,outNames);

  trans_parms.deallocate();
  amrex::Finalize();
  return 0;
}
