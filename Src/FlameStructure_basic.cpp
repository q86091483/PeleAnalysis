#include <string>
#include <iostream>
#include <algorithm>
#include <map>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_BLFort.H>
#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>
#include <mechanism.H>
#include <PelePhysics.H>

#include <H5Cpp.h>

using namespace amrex;

// Set write precision
#ifdef BL_USE_DOUBLE
#    define H5T_REAL H5::PredType::NATIVE_DOUBLE
#else
#    define H5T_REAL H5::PredType::NATIVE_FLOAT
#endif
// Integer and string type definition
#define H5T_INT      H5::PredType::NATIVE_INT
#define H5T_UINT     H5::PredType::NATIVE_UINT
#define H5T_STR      H5::StrType(H5::PredType::C_S1, H5T_VARIABLE)

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

enum weight_t
{
    w_none = 0,
    w_volume = 1,
    w_mass = 2,
    w_density = 3,
    w_notdefined = -1
};
const static Vector<std::string> weight_map = {"none", "volume", "mass", "density"};

enum bin_t
{
    b_lin = 0,
    b_log = 1
};
const static Vector<std::string> bin_map = {"lin", "log"};


void writeIntData
(
    H5::Group& h5g,
    const std::string& name,
    const Vector<int>& nBins,
    const Vector<int>& values,
    const int total,
    const weight_t weight = w_notdefined,
    const int cmpr_lvl = 0
) 
{
    // Modify dataset creation property to enable chunking and compression
    H5::DSetCreatPropList plist;
    if (cmpr_lvl > 0) {
        hsize_t chunk_dims[nBins.size()];
        for (int i=0; i<nBins.size(); i++) {
            int c = nBins[i];
            while (c > 22) {
                c = c/2 + 1;
            }
            chunk_dims[i] = c;
        }
        plist.setChunk(nBins.size(), chunk_dims);
        plist.setDeflate(cmpr_lvl); // ZLIB compression
    }

    // Create and write dataspace and dataset
    hsize_t dims[nBins.size()];
    for (int i=0; i<nBins.size(); i++) {
        dims[i] = nBins[i];
    }
    H5::DataSpace dataspace(nBins.size(), dims);

    H5::DataSet dataset
    (
        h5g.createDataSet(name, H5T_UINT, dataspace, plist)
    );
    dataset.write(values.data(), H5T_UINT);

    // Create total attribute
    H5::Attribute attributeTot =
        dataset.createAttribute("total", H5T_UINT, H5S_SCALAR);
    attributeTot.write(H5T_UINT, &total);

    // Create weight attribute
    if (weight != w_notdefined) {
        H5::Attribute attribute =
            dataset.createAttribute("weight", H5T_STR, H5S_SCALAR);
        attribute.write(H5T_STR, weight_map[weight]);
    }
} // void writeIntData


void writeRealData
(
    H5::Group& h5g,
    const std::string& name,
    const Vector<int>& nBins,
    const Vector<Real>& values,
    const Real total,
    const weight_t weight = w_notdefined,
    const int cmpr_lvl = 0
)
{
    // Modify dataset creation property to enable chunking and compression
    H5::DSetCreatPropList plist;
    if (cmpr_lvl > 0) {
        hsize_t chunk_dims[nBins.size()];
        for (int i=0; i<nBins.size(); i++) {
            int c = nBins[i];
            while (c > 22) {
                c = c/2 + 1;
            }
            chunk_dims[i] = c;
        }
        plist.setChunk(nBins.size(), chunk_dims);
        plist.setDeflate(cmpr_lvl); // ZLIB compression
    }

    // Create and write dataspace and dataset
    hsize_t dims[nBins.size()];
    for (int i=0; i<nBins.size(); i++) {
        dims[i] = nBins[i];
    }
    H5::DataSpace dataspace(nBins.size(), dims);

    H5::DataSet dataset
    (
        h5g.createDataSet(name, H5T_REAL, dataspace, plist)
    );
    dataset.write(values.data(), H5T_REAL);

    // Create total attribute
    H5::Attribute attributeTot =
        dataset.createAttribute("total", H5T_REAL, H5S_SCALAR);
    attributeTot.write(H5T_REAL, &total);

    // Create weight attribute
    if (weight != w_notdefined) {
        H5::Attribute attribute =
            dataset.createAttribute("weight", H5T_STR, H5S_SCALAR);
        attribute.write(H5T_STR, weight_map[weight]);
    }
} // void writeRealData


inline int computeBin
(
    const Real value,
    const Real lowerBound,
    const Real upperBound,
    const int N,
    const bin_t type
)
{
    int idx;
    switch (type)
    {
        case b_lin:
        {
            idx = floor((value - lowerBound)/(upperBound - lowerBound)*N);
            break;
        }
        case b_log:
        {
            idx = floor(std::log10(value/lowerBound)/std::log10(upperBound/lowerBound)*N);
            break;
        }
    }

    return std::max(std::min(idx, N - 1), 0);
} // computeBin


static void 
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
} // basename


int main (int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  ParmParse pp;
  const bool isioproc = ParallelDescriptor::IOProcessor();
  const int ioproc = ParallelDescriptor::IOProcessorNumber();

  if (argc < 2) {
    print_usage(argc,argv);
  }

  if (pp.contains("help")) {
    print_usage(argc,argv);
  }

  int verbose = 0;
  pp.query("verbose", verbose);
  if (verbose > 2) {
      AmrData::SetVerbose(true);
  }
  if (! isioproc) {
      verbose = 0;
  }

  // Number of files to read
  const int nPlotFiles = pp.countval("infile");
  amrex::Print() << nPlotFiles << std::endl;
  if (nPlotFiles < 1) {
    Print(ioproc) << "Bad nPlotFiles, exiting ..." << std::endl; 
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  } 
  const int nPlotFilesDerived = pp.countval("infile_derived");
  if (nPlotFilesDerived < 1) {
    Print(ioproc) << "No derived field is needed." << std::endl; 
  }
  if (verbose > 0) {
    Print(ioproc) << "Processing " << nPlotFiles << " + " << nPlotFilesDerived << " plotfiles ..." << std::endl;
  }

  // Plot file names
  amrex::Vector<std::string> plotFileNames(nPlotFiles);
  for (int i = 0; i < nPlotFiles; ++i) {
    pp.get("infile", plotFileNames[i], i);
    if (verbose > 1) {
      amrex::Print(ioproc) << "   " << basename(plotFileNames[i]) << std::endl;
    }
  }

  // Plot file names for derived fields
  amrex::Vector<std::string> plotFileNamesDerived;
  for (int i = 0; i < nPlotFilesDerived; ++i) {
    std::string strt;
    pp.get("infile_derived", strt, i);
    plotFileNamesDerived.push_back(strt);
    if (verbose > 1) {
      amrex::Print(ioproc) << "    " << plotFileNamesDerived[i] << std::endl;
    }
  }

  // Settings of HDF IO
  std::string outDir = "cond_ISRN";
  pp.query("outputDir", outDir);
  UtilCreateDirectory(outDir, 0755);
  bool writeGrid = false;
  pp.query("writeGrid", writeGrid);
  bool writeCoordinates = false;
  pp.query("writeCoordinates", writeCoordinates);
  std::string outputLabel;
  pp.query("outputLabel", outputLabel);

  // Finest level
  int finestLevel_in(-1);
  pp.query("finestLevel", finestLevel_in);
  int Nlev_in = finestLevel_in + 1;

  // Get species names
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);

  // Initialize transport data
  auto eos = pele::physics::PhysicsType::eos();
  pele::physics::transport::TransportParams<
      pele::physics::PhysicsType::transport_type>
      trans_parms;
  trans_parms.allocate();

  // Get the transport data pointer
  auto const* ltransparm = trans_parms.device_trans_parm();

  // Initialize dataServicePtrVec & amrDataPtrVec
  amrex::Print() << "Initializing dataServicePtrVec: " << std::endl;
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  amrex::Vector<DataServices*>  dataServicePtrVec(nPlotFiles);                                         
  amrex::Vector<AmrData*>           amrDataPtrVec(nPlotFiles);
  amrex::Vector<Real>                        time(nPlotFiles);
  for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {
    dataServicePtrVec[iPlot] = new DataServices(plotFileNames[iPlot], fileType);
    if( ! dataServicePtrVec[iPlot]->AmrDataOk()) {                        
      amrex::Print() << "   " << "Initialize dataServicePtrVec failed for " 
                     << plotFileNames[iPlot] << std::endl; 
	    DataServices::Dispatch(DataServices::ExitRequest, NULL); 
    }
    amrDataPtrVec[iPlot] = &(dataServicePtrVec[iPlot]->AmrDataRef());
    time[iPlot] = amrDataPtrVec[iPlot]->Time();
  }
  amrex::Print() << "   Done." << std::endl;

  std::string fn;

  // Input index and names
  amrex::Vector<std::string> inNames;
  std::map<std::string,int> mi;
  int nCompIn = 0;
  amrex::Vector<int> destFillComps;

  const int IYSP = 0;
  for (const auto &spn : spec_names) {
    fn = "Y(" + spn + ")";
    inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1;
  }
  fn = "temp";  
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int ITEMP = mi[fn];
  fn = "HeatRelease";  
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IHRR = mi[fn];
  fn = "density";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IRHO = mi[fn];
  fn = "x_velocity";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IU = mi[fn];
  fn = "y_velocity";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IV = mi[fn];
  fn = "z_velocity"; 
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IW = mi[fn];

  for (int i = 0; i < nCompIn; i++) {
    destFillComps.emplace_back(i);
  }
  amrex::Print() << "Read " << nCompIn 
    << " fields (inNames[i] -> destFillComps[i]):" << std::endl;
  for (int i=0; i<inNames.size(); ++i) {
    amrex::Print() << "   " << inNames[i] << " -> " << destFillComps[i] << std::endl;
  }

  // Intermediate fields (middle) index and names
  amrex::Vector<std::string> midNames;
  std::map<std::string, int> mm;
  int nCompMid = 0;
  fn = "x";
  midNames.emplace_back(fn); nCompMid = midNames.size(); mi[fn] = nCompMid + nCompIn - 1; 
  fn = "y";
  midNames.emplace_back(fn); nCompMid = midNames.size(); mi[fn] = nCompMid + nCompIn - 1;
  fn = "z";
  midNames.emplace_back(fn); nCompMid = midNames.size(); mi[fn] = nCompMid + nCompIn - 1; 

  // Output index and names (Y of <Y|X>)
  amrex::Vector<std::string> outNames;
  int nCompOut;
  std::map<std::string, int> mo;
  fn = "rho";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "mixture_fraction";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "temp";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "HeatRelease";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "ts11";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "ts22";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "ts33";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "ts12";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "ts13";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "ts23";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "ts_sum";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "mu";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 

  // Conditional variables (X of <Y|X>)
  const int nVars(pp.countval("vars"));
  if (nVars < 1) {
    Error("Needs to specify at least one conditioning variable.");
  }
  amrex::Vector<std::string>          varNames(nVars);
  amrex::Vector<int>                  bins(nVars);
  amrex::Vector<int>                  nBins(nVars);
  amrex::Vector<bin_t>                binType(nVars, b_lin);
  amrex::Vector<int>                  incOutOfBounds(nVars, 0);
  amrex::Vector<std::array<Real,2>>   varBounds(nVars);
  for (int i = 0; i < nVars; ++i) {
    pp.get("vars",      varNames[i],    i);
    pp.get("nBins",     nBins [i],      i);
    Real lo, hi;
    pp.get("varBounds", lo,             i*2);
    pp.get("varBounds", hi,             (i*2)+1);
    varBounds[i] = {lo, hi};
  }
  if (pp.contains("binType")) {
    for (int i = 0; i < nVars; i++) {
      std::string tmp;
      pp.get("binType", tmp, i);
      if (tmp == "lin") {
        binType[i] = b_lin;
      } else if (tmp == "log") {
        binType[i] = b_log;
      } else {
        std::string s = "Unknown binType: " + tmp;
        Error(s.c_str());
      }
    }
  }
  if (pp.contains("incOutOfBounds")) {
    for (int i = 0; i < nVars; i++) {
      pp.get("incOutOfBounds", incOutOfBounds[i], i);
    }
  }
  if (verbose > 0) {
    std::cout << "Bounds for bins:" << std::endl;
    for (int i=0; i<nVars; i++) {
      std::cout << "   " << varNames[i]
                << ": [" << varBounds[i][0]
                << ", " << varBounds[i][1]
                << "] (type = " << binType[i]
                << ", incOutOfBounds = " << incOutOfBounds[i] << ")"
                << std::endl;
    }
  }
  int nBinsTot = 1;
  for (int i = 0; i < nVars; i++) {
    nBinsTot *= nBins[i];
  }

  // Variables to be averaged (Y of <Y|X>)
  amrex::Vector<std::string> avgVarNames;
  amrex::Vector<weight_t> avgVarWeights;
  int nAvgVars = nCompOut;
  for (int i = 0; i < nAvgVars; i++) {
    avgVarNames.push_back(outNames[i]);
    avgVarWeights.push_back(w_volume); 
  }

  Vector<Real> dataX(nVars);
  Vector<Real> dataY(nAvgVars);

  // Iterate over input plot files
  for (int iPlot; iPlot < nPlotFiles; ++iPlot) {

    amrex::Print() << "Processing " << iPlot << std::endl;
    AmrData& amrData = dataServicePtrVec[iPlot]->AmrDataRef();

    Vector<Real> probLo = amrData.ProbLo();
    Vector<Real> probHi = amrData.ProbHi();

    int finestLevel = amrData.FinestLevel();
    if (finestLevel > finestLevel_in) finestLevel = finestLevel_in;
    int Nlev = finestLevel + 1;
    amrex::RealBox real_box({AMREX_D_DECL(amrData.ProbLo()[0],
                                          amrData.ProbLo()[1],
                                          amrData.ProbLo()[2])},
                            {AMREX_D_DECL(amrData.ProbHi()[0],
                                          amrData.ProbHi()[1],
                                          amrData.ProbHi()[2])});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 1, 0)};
    Vector<Geometry> geoms(Nlev);
    geoms[0] = amrex::Geometry((amrData.ProbDomain())[0], real_box, amrData.CoordSys(),
                                is_periodic);
    const int nGrow = 1;

    // Vector of MultiFab for output
    Vector<MultiFab> mfv_out(Nlev);

    // Arrays for storing statistics
    int ratio           = 1;
    int bin             = 0;
    bool skip           = false;
    int ii;
    Real vol, rho, m, weight;
    Real vol0           = 1;
    Vector<int> count(nBinsTot, 0);
    Vector<Real> volMean(nBinsTot, 0);
    Vector<Real> volStd(nBinsTot, 0);
    Vector<Real> rhoMean(nBinsTot, 0);
    Vector<Real> rhoStd(nBinsTot, 0);
    Vector<Real> massMean(nBinsTot, 0);
    Vector<Real> massStd(nBinsTot, 0);
    Vector<Vector<Real>> varMeanVal(nAvgVars, Vector<Real>(nBinsTot, 0));
    Vector<Vector<Real>> varStdVal(nAvgVars, Vector<Real>(nBinsTot, 0));
    Vector<Vector<Real>> varMinVal(nAvgVars, Vector<Real>(nBinsTot, 1e30));
    Vector<Vector<Real>> varMaxVal(nAvgVars, Vector<Real>(nBinsTot, -1e30));
    Vector<int> countLvl(Nlev, 0);
    Vector<Real> volLvl(Nlev, 0);

    // Loop over levels
    for (int lev=0; lev<Nlev; ++lev) {  
      // Level info 
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      //const Vector<Real>& dx = amrData.DxLevel()[lev];
      const Vector<Real> dx = amrData.CellSize(lev);
      if ( lev > 0 ) {
        geoms[lev] = amrex::refine(geoms[lev - 1], 2);
      }

      // Input MultiFabs
      MultiFab mf_in(ba, dm, nCompIn, nGrow);
      Print() << "   - Reading data (FillVar) for level " << lev << std::endl;
      amrData.FillVar(mf_in, lev, inNames, destFillComps);
      Print() << "   - Data has been loaded for level " << lev << std::endl;
      // Fix up grow cells.  Use extrap for guess
      const Box& dbox = amrData.ProbDomain()[lev];
      for (amrex::MFIter mfi(mf_in); mfi.isValid(); ++mfi) {
        FArrayBox& fab = mf_in[mfi];
        const Box& box = mfi.validbox();
        pushvtog(BL_TO_FORTRAN_BOX(box),
                 BL_TO_FORTRAN_BOX(dbox),
                BL_TO_FORTRAN_ANYD(fab),
                &nCompIn);
      }
      // Fix up fine-fine and periodic
      mf_in.FillBoundary(geoms[lev].periodicity());

      // Output MultiFabs
      mfv_out[lev] = MultiFab(ba,dm,nCompOut,nGrow);
      //MultiFab mf_out(mfv_out[lev], amrex::make_alias, 0, nCompOut);

      // Intermediate fields MultiFabs
      MultiFab mf_xyz(ba, dm, 3, nGrow);

      MultiFab mf_D(ba, dm, NUM_SPECIES, nGrow);
      MultiFab mf_lam(ba, dm, 1, nGrow);
      MultiFab mf_chi(ba, dm, 1, nGrow);
      MultiFab mf_xi(ba, dm, 1, nGrow);
      MultiFab mf_mu(ba, dm, 1, nGrow);
      MultiFab mf_gradu(ba, dm, BL_SPACEDIM+1, nGrow);
      MultiFab mf_gradv(ba, dm, BL_SPACEDIM+1, nGrow);
      MultiFab mf_gradw(ba, dm, BL_SPACEDIM+1, nGrow);
		  MultiFab mf_su(ba, dm, 6, nGrow); // symmetric velocity tensor 0.5(du_i/dx_j+du_j/dx_i)
		  MultiFab mf_tau(ba, dm, 6, nGrow); // viscous force tensor
      const int L11=0, L22=1, L33=2, L12=3, L13=4, L23=5;
      const int                      L21=3, L31=4, L32=5;
      const int LSUM = 6;

      // Mask if covered by finer level - mf_covered
      MultiFab mf_covered(ba,dm,1,nGrow); 
      mf_covered.setVal(0.0);
      if (lev < finestLevel) {
        auto rref = amrData.RefRatio()[lev];
        BoxArray baf = BoxArray(amrData.boxArray(lev+1)).coarsen(rref);
        for (MFIter mfi(mf_covered); mfi.isValid(); ++mfi) {
          const Box& box = mfi.validbox();
          FArrayBox& fab = mf_covered[mfi];
          std::vector<std::pair<int,Box>> isects = baf.intersections(box);
          for (int ii = 0; ii < isects.size(); ii++)
            fab.setVal(1.0, isects[ii].second);
        }
      }

      // Volume for the current level
      vol = 1;
      for (int i = 0; i < BL_SPACEDIM; i++)
        vol *= dx[i];
      // Normalise cell volume for numerical accuracy
      if (lev == 0)
        vol0 = vol;
      vol /= vol0;
      if (verbose > 1) {
        std::cout << "   - collecting statistics" << std::endl;
        std::cout << "      - dx, dy, dz [m] =  " << dx[0] << ", " << dx[1] << ", " << dx[2] << "\n"
                  << "      - cell volume  at level 0=  " << vol*vol0 << " [m^3] / "
                  << vol << " (norm)."
                  << std::endl;
      }

      //FArrayBox alias_fab(orig_fab, amrex::make_alias, 1, 2);
      // Calculate result variables & collect their statistics

      for (amrex::MFIter mfi(mf_in, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        // Reference FArraybox
        FArrayBox& fab_in = mf_in[mfi];
        FArrayBox& fab_out = mfv_out[lev][mfi];
        FArrayBox& fab_covered = mf_covered[mfi];
        FArrayBox& fab_xyz = mf_xyz[mfi];

        FArrayBox& fab_gradu = mf_gradu[mfi];
        FArrayBox& fab_gradv = mf_gradv[mfi];
        FArrayBox& fab_gradw = mf_gradw[mfi];
        FArrayBox& fab_su    = mf_su[mfi];
        FArrayBox& fab_tau   = mf_tau[mfi];

        // Reference Array4
        Array4<Real> const& rho_a   = mf_in.array(mfi, mi["density"]);
        Array4<Real> const& mixfrac_a  = mf_in.array(mfi, mi["mixture_fraction"]);
        Array4<Real> const& Y_a     = mf_in.array(mfi, mi[spec_names[0]]);
        Array4<Real> const& T_a     = mf_in.array(mfi, mi["temp"]);
        Array4<Real> const& HRR_a   = mf_in.array(mfi, mi["HeatRelease"]);
        Array4<Real> const& U_a     = mf_in.array(mfi, mi["x_velocity"]);
        Array4<Real> const& V_a     = mf_in.array(mfi, mi["y_velocity"]);
        Array4<Real> const& W_a     = mf_in.array(mfi, mi["z_velocity"]);
        Array4<Real> const& x_a     = mf_xyz.array(mfi, 0);
        Array4<Real> const& y_a     = mf_xyz.array(mfi, 1);
        Array4<Real> const& z_a     = mf_xyz.array(mfi, 2);  
        Array4<Real> const& covered_a = mf_covered.array(mfi);

        Array4<Real> const& D_a     = mf_D.array(mfi);
        Array4<Real> const& mu_a    = mf_mu.array(mfi);
        Array4<Real> const& lam_a   = mf_lam.array(mfi);
        Array4<Real> const& chi_a   = mf_chi.array(mfi);
        Array4<Real> const& xi_a    = mf_xi.array(mfi);
        Array4<Real> const& su_a    = mf_su.array(mfi);
        Array4<Real> const& tau_a   = mf_tau.array(mfi);
        Array4<Real> const& gradu_a = mf_gradu.array(mfi);
        Array4<Real> const& gradv_a = mf_gradv.array(mfi);
        Array4<Real> const& gradw_a = mf_gradw.array(mfi);

        Array4<Real> const& rho_out_a = mfv_out[lev].array(mfi, mo["rho"]);
        Array4<Real> const& mixfrac_out_a = mfv_out[lev].array(mfi, mo["mixture_fraction"]);
        Array4<Real> const& T_out_a = mfv_out[lev].array(mfi, mo["temp"]);
        Array4<Real> const& HRR_out_a = mfv_out[lev].array(mfi, mo["HeatRelease"]);
        Array4<Real> const& mu_out_a = mfv_out[lev].array(mfi, mo["mu"]);
        Array4<Real> const& ts_a      = mfv_out[lev].array(mfi, mo["ts11"]);

        // Calculate coordiante
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                  x_a(i,j,k) = amrData.ProbLo()[0] + i * dx[0];
                  y_a(i,j,k) = amrData.ProbLo()[1] + j * dx[1];
                  z_a(i,j,k) = amrData.ProbLo()[2] + k * dx[2];
              }
          }
        }

        // Calculate transport coefficients
        amrex::launch(bx, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
          auto trans = pele::physics::PhysicsType::transport();
          trans.get_transport_coeffs(
            tbx, Y_a, T_a, rho_a, D_a, chi_a, mu_a, xi_a, lam_a, ltransparm);
        });

        // Calculate gradients of velocity
        gradient(BL_TO_FORTRAN_BOX(bx),
                BL_TO_FORTRAN_N_ANYD(fab_in,    mi["x_velocity"]),
                BL_TO_FORTRAN_N_ANYD(fab_gradu, 0),
                &(dx[0]));
        gradient(BL_TO_FORTRAN_BOX(bx),
                BL_TO_FORTRAN_N_ANYD(fab_in,    mi["y_velocity"]),
                BL_TO_FORTRAN_N_ANYD(fab_gradv, 0),
                &(dx[0]));
        gradient(BL_TO_FORTRAN_BOX(bx),
                BL_TO_FORTRAN_N_ANYD(fab_in,    mi["z_velocity"]),
                BL_TO_FORTRAN_N_ANYD(fab_gradw, 0),
                &(dx[0]));

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          rho_out_a(i,j,k) = rho_a(i,j,k);
          mixfrac_out_a(i,j,k) = mixfrac_a(i,j,k);
          T_out_a(i,j,k) = T_a(i,j,k);
          HRR_out_a(i,j,k) = HRR_a(i,j,k);
          mu_out_a(i,j,k) = mu_a(i,j,k);

          su_a(i,j,k,L11) = gradu_a(i,j,k,0);
          su_a(i,j,k,L22) = gradv_a(i,j,k,1);
          su_a(i,j,k,L33) = gradw_a(i,j,k,2);
          su_a(i,j,k,L12) = 0.5 * (gradu_a(i,j,k,1) + gradv_a(i,j,k,0));
          su_a(i,j,k,L13) = 0.5 * (gradu_a(i,j,k,2) + gradw_a(i,j,k,0));
          su_a(i,j,k,L23) = 0.5 * (gradv_a(i,j,k,2) + gradw_a(i,j,k,1));

          Real skk = (1/3.) * (gradu_a(i,j,k,0)+gradv_a(i,j,k,1)+gradw_a(i,j,k,2));
          tau_a(i,j,k,L11) = 2. * mu_a(i,j,k) * (su_a(i,j,k,L11) - skk);
          tau_a(i,j,k,L22) = 2. * mu_a(i,j,k) * (su_a(i,j,k,L22) - skk);
          tau_a(i,j,k,L33) = 2. * mu_a(i,j,k) * (su_a(i,j,k,L33) - skk);
          tau_a(i,j,k,L12) = 2 * mu_a(i,j,k) * su_a(i,j,k,L12);
          tau_a(i,j,k,L13) = 2 * mu_a(i,j,k) * su_a(i,j,k,L13); 
          tau_a(i,j,k,L23) = 2 * mu_a(i,j,k) * su_a(i,j,k,L23); 

          Real coeff = 2 * mu_a(i,j,k) / rho_a(i,j,k);
          ts_a(i,j,k,L11) = su_a(i,j,k,L11);
          ts_a(i,j,k,L22) = su_a(i,j,k,L22);
          ts_a(i,j,k,L33) = su_a(i,j,k,L33);
          ts_a(i,j,k,L12) = su_a(i,j,k,L12);
          ts_a(i,j,k,L13) = su_a(i,j,k,L13);
          ts_a(i,j,k,L23) = su_a(i,k,k,L23);

          ts_a(i,j,k,LSUM) = ts_a(i,j,k,L11) * ts_a(i,j,k,L11) + 
                             ts_a(i,j,k,L12) * ts_a(i,j,k,L12) + 
                             ts_a(i,j,k,L13) * ts_a(i,j,k,L13) +
                             ts_a(i,j,k,L21) * ts_a(i,j,k,L21) + 
                             ts_a(i,j,k,L22) * ts_a(i,j,k,L22) +
                             ts_a(i,j,k,L23) * ts_a(i,j,k,L23) +
                             ts_a(i,j,k,L31) * ts_a(i,j,k,L31) +
                             ts_a(i,j,k,L32) * ts_a(i,j,k,L32) + 
                             ts_a(i,j,k,L33) * ts_a(i,j,k,L33);


        });

        // Collect statistics for mfi
        countLvl[lev] += bx.volume();
        for (IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
          //for (int i = 0; i < nVars; i++) {
          //  dataX[i] = mf_in.array(mfi, mi[varNames[i]])(iv[0],iv[1],iv[2]);
          //}
          //dataX[0] = mf_x.array(mfi)(iv[0], iv[1], iv[2]);
          //dataX[1] = mf_y.array(mfi)(iv[0], iv[1], iv[2]);
          //dataX[2] = mf_z.array(mfi)(iv[0], iv[1], iv[2]);

          fab_xyz.getVal(dataX.dataPtr(), iv);
          fab_out.getVal(dataY.dataPtr(), iv);
          //Skip points overlapping with a finer levels
          skip = false;
          if (covered_a(iv[0],iv[1],iv[2]) > 0.0) skip = true;
          for (int i=0; i<nVars; i++) {
            if ((incOutOfBounds[i] == 0) &&
                ((dataX[i] < varBounds[i][0]) ||
                 (dataX[i] >= varBounds[i][1]))) {
              skip = true;
            }
          }
          if (skip) continue;
          // Compute the bin
          for (int i=0; i<nVars; i++) {
            bins[i] = computeBin(dataX[i], varBounds[i][0], varBounds[i][1], nBins[i], binType[i]);
          }
          bin = 0;
          for (int i=0; i<nVars-1; i++) {
            bin += bins[i];
            bin *= nBins[i+1];
          }
          bin += bins[nVars-1];

          // Gather statistics
          rho = dataY[mo["density"]];
          m = rho * vol;
          count[bin] += 1;
          volMean[bin] += vol;
          volStd[bin] += vol * vol;
          rhoMean[bin] += rho;
          rhoStd[bin] += rho * rho;
          massMean[bin] += m;
          massStd[bin] += m * m;
          for (int i=0; i<nCompOut; i++) {
            varMeanVal[i][bin] += dataY[i] * vol;
            varStdVal[i][bin] += dataY[i] * dataY[i] * vol;
            varMinVal[i][bin] += std::min(varMinVal[i][bin], dataY[i]);
            varMaxVal[i][bin] += std::min(varMaxVal[i][bin], dataY[i]);
          }
        } // iv of mfi to collect statistics
      } // mfi - set variables && collect stats

      // Cell volume for each level
      volLvl[lev] = countLvl[lev] * vol;
    } // lev - set variables && collect stats
    ParallelDescriptor::ReduceIntSum(count.dataPtr(), count.size(), ioproc);
    ParallelDescriptor::ReduceRealSum(volMean.dataPtr(), volMean.size(), ioproc);
    ParallelDescriptor::ReduceRealSum(volStd.dataPtr(), volStd.size(), ioproc);
    ParallelDescriptor::ReduceRealSum(rhoMean.dataPtr(), rhoMean.size(), ioproc);
    ParallelDescriptor::ReduceRealSum(rhoStd.dataPtr(), rhoStd.size(), ioproc);
    ParallelDescriptor::ReduceRealSum(massMean.dataPtr(), massMean.size(), ioproc);
    ParallelDescriptor::ReduceRealSum(massStd.dataPtr(), massStd.size(), ioproc);
    for (int i=0; i<nAvgVars; i++) {
      ParallelDescriptor::ReduceRealSum(varMeanVal[i].dataPtr(), varMeanVal[i].size(), ioproc);
      ParallelDescriptor::ReduceRealSum(varStdVal[i].dataPtr(), varStdVal[i].size(), ioproc);
      ParallelDescriptor::ReduceRealMin(varMinVal[i].dataPtr(), varMinVal[i].size(), ioproc);
      ParallelDescriptor::ReduceRealMax(varMaxVal[i].dataPtr(), varMaxVal[i].size(), ioproc);
    }
    ParallelDescriptor::ReduceIntSum(countLvl.dataPtr(), countLvl.size(), ioproc);
    ParallelDescriptor::ReduceRealSum(volLvl.dataPtr(), volLvl.size(), ioproc);

    // Write statistics
    if (isioproc) {
      std::string h5fileName = outDir + "/" + basename(plotFileNames[iPlot]) + "_" + outputLabel + ".h5";
      if (verbose > 1)
        std::cout << "  Write HDF5 file: " << h5fileName << std::endl;
      H5::H5File h5file(h5fileName, H5F_ACC_TRUNC);

      // Write time
      {
        Real time = amrData.Time();
        hsize_t dims[1];
        dims[0] = 1;
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset(h5file.createDataSet("time", H5T_REAL, dataspace));
        dataset.write(&time, H5T_REAL);
      }

      // Write grid parameter (e.g. for VTK/XDMF)
      if (writeGrid) {
        H5::Group h5grid(h5file.createGroup("GRID"));

        // Write DIMENSIONS
        {
          uint DIMENSIONS[nVars];
          for (int i = 0; i < nVars; i++) {
            DIMENSIONS[nVars - (i + 1)] = nBins[i] + 1;
          }

          hsize_t dims[1];
          dims[0] = nVars;
          H5::DataSpace dataspace(1, dims);

          H5::DataSet dataset(
              h5grid.createDataSet("DIMENSIONS", H5T_UINT, dataspace));
          dataset.write(DIMENSIONS, H5T_UINT);
        }

        // Write ORIGIN
        {
          Real ORIGIN[nVars];
          for (int i = 0; i < nVars; i++) {
            ORIGIN[nVars - (i + 1)] = varBounds[i][0];
          }

          hsize_t dims[1];
          dims[0] = nVars;
          H5::DataSpace dataspace(1, dims);

          H5::DataSet dataset(
              h5grid.createDataSet("ORIGIN", H5T_REAL, dataspace));
          dataset.write(ORIGIN, H5T_REAL);
        }

        // Write SPACING
        {
          Real SPACING[nVars];
          for (int i = 0; i < nVars; i++) {
            SPACING[nVars - (i + 1)] =
                (varBounds[i][1] - varBounds[i][0]) / std::max(nBins[i], 1);
          }

          hsize_t dims[1];
          dims[0] = nVars;
          H5::DataSpace dataspace(1, dims);

          H5::DataSet dataset(
              h5grid.createDataSet("SPACING", H5T_REAL, dataspace));
          dataset.write(SPACING, H5T_REAL);
        }
        h5grid.close();
      } // Write grid

      if (writeCoordinates) {
        H5::Group h5coord(h5file.createGroup("COORDS"));
        for (int i = 0; i < nVars; i++) {
          uint d = nBins[i] + 1;
          Real o = varBounds[i][0];
          Real coord[d];
          switch (binType[i]) {
          case b_lin: {
            Real s =
                (varBounds[i][1] - varBounds[i][0]) / std::max(nBins[i], 1);
            for (int j = 0; j < d; j++) {
              coord[j] = Real(j) * s + o;
            }
            break;
          }
          case b_log: {
            Real s = std::log10(varBounds[i][1] / varBounds[i][0]) /
                     std::max(nBins[i], 1);
            for (int j = 0; j < d; j++) {
              coord[j] = std::pow(10, Real(j) * s + std::log10(o));
            }
            break;
          }
          }

          // if (incOutOfBounds[i] > 0) {
          //     coord[0] = varGMin[i];
          //     coord[d-1] = varGMax[i];
          // }

          hsize_t dims[1];
          dims[0] = d;
          H5::DataSpace dataspace(1, dims);

          char name[10];
          snprintf(name, sizeof(name), "C%d", i);
          H5::DataSet dataset(h5coord.createDataSet(name, H5T_REAL, dataspace));
          dataset.write(coord, H5T_REAL);

          H5::Attribute attribute;
          attribute = dataset.createAttribute("name", H5T_STR, H5S_SCALAR);
          attribute.write(H5T_STR, varNames[i]);

          attribute = dataset.createAttribute("type", H5T_STR, H5S_SCALAR);
          attribute.write(H5T_STR, bin_map[binType[i]]);

          attribute =
              dataset.createAttribute("incOutOfBounds", H5T_INT, H5S_SCALAR);
          attribute.write(H5T_INT, &incOutOfBounds[i]);

        }
        h5coord.close();
      } // Write coordinate

      // Global (unconditional) data
      H5::Group h5domain(h5file.createGroup("DOMAIN"));
      // Domain size
      writeRealData(h5domain, "prob_lo", Vector<int>(1, BL_SPACEDIM), probLo,
                    0.0);
      writeRealData(h5domain, "prob_hi", Vector<int>(1, BL_SPACEDIM), probHi,
                    0.0);

      // Count/ volume per level (unconditional)
      writeIntData(h5domain, "count_lvl", Vector<int>(1, Nlev), countLvl,
                   0.0);
      for (auto &v : volLvl)
        v *= vol0;
      writeRealData(h5domain, "volume_lvl", Vector<int>(1, Nlev), volLvl,
                    0.0);
      h5domain.close();

      // Conditional data
      H5::Group h5data(h5file.createGroup("DATA"));
      for (int i = 0; i < nAvgVars; i++) {
        writeRealData(h5data, avgVarNames[i] + "_mean", nBins, varMeanVal[i],
                      0.0, avgVarWeights[i], 0.0);
        writeRealData(h5data, avgVarNames[i] + "_std", nBins, varStdVal[i], 0.0,
                      avgVarWeights[i], 4);
        writeRealData(h5data, avgVarNames[i] + "_min", nBins, varMinVal[i], 0.0,
                      w_notdefined, 4);
        writeRealData(h5data, avgVarNames[i] + "_max", nBins, varMaxVal[i], 0.0,
                      w_notdefined, 4);
      }

      // Count
      int countTot = 0;
      for (const auto &c : count)
        countTot += c;
      writeIntData(h5data, "count", nBins, count, countTot, w_notdefined, 4);

      // Volume
      Real volTot = 0;
      for (const auto &v : volMean) volTot += v;
      volTot *= vol0;
      for (auto &v : volMean) v *= 1.0; //vol0;
      for (auto &v : volStd) v *= vol0*vol0;
      writeRealData(h5data, "volume_sum", nBins, volMean, volTot, w_notdefined, 4);
      writeRealData(h5data, "volume_mean", nBins, volMean, 0.0, w_none, 4);
      writeRealData(h5data, "volume_std", nBins, volStd, 0.0, w_none, 4);

      h5data.close();
      h5file.close();
    } // statistics

    // Write output fields
    if (true) {
      std::string outfilename(outDir + "/" + basename(plotFileNames[iPlot]) + "_derived"); 
      Print() << "Writing ISRN derived data to " << outfilename << std::endl;
      //bool verb = false;
	    amrex::Vector<int> istep;
      amrex::Vector<IntVect> ref_ratio;
      for (int lev=0; lev<Nlev; ++lev) {
		    istep.emplace_back(1);
		    IntVect iv(AMREX_D_DECL(2,2,2));
        ref_ratio.emplace_back(iv);
      }
	    WriteMultiLevelPlotfile(outfilename, Nlev, 
          GetVecOfConstPtrs(mfv_out), outNames, geoms, 0.0, istep, ref_ratio);
    }
  } // iPlot
} // end of main