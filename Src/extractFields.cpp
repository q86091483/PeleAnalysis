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

  // Nlev to process
  int finestLevel = amrData.FinestLevel();
  pp.query("extract_finest_level", finestLevel);
  int Nlev = finestLevel + 1;

  // Fields to extract
  const int nVars(pp.countval("extract_vars"));
  amrex::Vector<std::string> extract_names(nVars);
  for (int i = 0; i < nVars; ++i) {
    pp.get("extract_vars", extract_names[i], i);
    amrex::Print() << "To extract - " << extract_names[i] << std::endl;
  }
  std::string extract_output_name;
  pp.get("extract_output_name", extract_output_name);

  const amrex::Vector<std::string>& plotVarNames = amrData.PlotVarNames();
  amrex::Print() << "The input file includes field:" << std::endl;
  for (int i=0; i<plotVarNames.size(); ++i)
  {
    amrex::Print() << plotVarNames[i] << std::endl;
  }

  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
  auto eos = pele::physics::PhysicsType::eos();

  const int nCompIn  = nVars;
  const int nCompOut = nVars;

  Vector<int> destFillComps(nCompIn);
  Vector<std::string> inNames(nCompIn);
  Vector<std::string> outNames(nCompOut);
  for (int i=0; i<nVars; ++i)
  {
    destFillComps[i]  = i;
    inNames[i]        = extract_names[i];
    outNames[i]       = extract_names[i];
  }

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
    const Vector<Real>& delta = amrData.DxLevel()[lev];
    if ( lev > 0 ) {
       geoms[lev] = amrex::refine(geoms[lev - 1], 2);
    }

    MultiFab indata(ba, dm, nCompIn, nGrow);

    // Loading data at lev.
    Print() << "Reading data (FillVar) for level " << lev << std::endl;
    amrData.FillVar(indata, lev, inNames, destFillComps);
    Print() << "Data has been read for level " << lev << std::endl;

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
      Array4<Real> const& in_a = indata.array(mfi, 0);
      Array4<Real> const& out_a = outdata[lev]->array(mfi, 0);

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          // Copy mixture fraction.
          out_a(i,j,k,0) = in_a(i,j,k,0); 
          out_a(i,j,k,1) = in_a(i,j,k,1); 
      });

    } // mfi
    Print() << "Extraction finished for level " << lev << std::endl;
  } // lev

  std::string outDir = "data_extracted"; pp.query("outputDir", outDir);
  std::string outfilename(outDir + "/" + basename(infile) + "_" + extract_output_name); 
  Print() << "Writing extracted data to " << outfilename << std::endl;
  bool verb = false;
  WritePlotFile(outdata, amrData, outfilename, verb, outNames);

  trans_parms.deallocate();
  amrex::Finalize();
  return 0;
}
