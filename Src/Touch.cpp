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

#define NUMMIXF 0
#define NUMAGE 0
#define NUMAGEPV 0
#define NUMAUX (NUMMIXF + NUMAGE + NUMAGEPV)

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

  // If write Derived field
  int writeDerivedField = 0;
  pp.query("writeDerivedField", writeDerivedField);

  // Get species names
  auto eos = pele::physics::PhysicsType::eos();
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
  amrex::Real atwCHON[4] = {0.0};
  pele::physics::eos::atomic_weightsCHON<pele::physics::PhysicsType::eos_type>(
    atwCHON);
  for (int i = 0; i < 4; i++) {
    amrex::Print() << atwCHON[i] << std::endl;
  }
  int ecompCHON[NUM_SPECIES * 4];
  pele::physics::eos::element_compositionCHON<
    pele::physics::PhysicsType::eos_type>(ecompCHON);
  amrex::Print() << "Number of atm (C,H,O,N) in species (molecular weight):" << std::endl;
  amrex::Real mwt[NUM_SPECIES];
  eos.molecular_weight(mwt);
  for (int i = 0; i < NUM_SPECIES; ++i) {
    for (int k = 0; k < 4; k++) {
      amrex::Print() << ecompCHON[i*4+k] << ", ";
    }
    amrex::Print() << mwt[i];
    amrex::Print() << std::endl;
  }

  amrex::Real Beta_mix[4] = {0.0};
  Beta_mix[0] = (atwCHON[0] != 0.0) ? 2.0 / atwCHON[0] : 0.0;
  Beta_mix[1] = (atwCHON[1] != 0.0) ? 1.0 / (2.0 * atwCHON[1]) : 0.0;
  Beta_mix[2] = (atwCHON[2] != 0.0) ? -1.0 / atwCHON[2] : 0.0;
  Beta_mix[3] = 0.0;

  amrex::Real spec_Bilger_fact[NUM_SPECIES] = {0.0};
  amrex::Real YF[NUM_SPECIES] = {0.0}; YF[H2_ID] = 1.0;
  amrex::Real YO[NUM_SPECIES] = {0.0}; YO[O2_ID] = 0.232; YO[N2_ID] = 1-YO[O2_ID];
  amrex::Real Zfu = 0.0;
  amrex::Real Zox = 0.0;
  for (int i = 0; i < NUM_SPECIES; ++i) {
    spec_Bilger_fact[i] = 0.0;
    for (int k = 0; k < 4; k++) {
      spec_Bilger_fact[i] +=
        Beta_mix[k] * (ecompCHON[i * 4 + k] * atwCHON[k] / mwt[i]);
    }
    Zfu += spec_Bilger_fact[i] * YF[i];
    Zox += spec_Bilger_fact[i] * YO[i];
  }
  amrex::Print() << "Zfu: " << Zfu << ", Zox: " << Zox << std::endl;

  // Initialize transport data
  pele::physics::transport::TransportParams<
      pele::physics::PhysicsType::transport_type>
      trans_parms;
  trans_parms.allocate();

  // Get the transport data pointer
  auto const* ltransparm = trans_parms.device_trans_parm();
  amrex::Real atw[4] = {0.0};

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

  // List of input fields
  std::string fn;
  amrex::Vector<std::string> inNames;
  std::map<std::string,int> mi;
  int nCompIn = 0;
  amrex::Vector<int> destFillComps;

  const int IYSP = 0;
  //for (const auto &spn : spec_names) {
  //  fn = "Y(" + spn + ")";
  //  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1;
  //}
  fn = "temp";  
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int ITEMP = mi[fn];
  //fn = "HeatRelease";  
  //inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IHRR = mi[fn];
  //fn = "density";
  //inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IRHO = mi[fn];
  //fn = "x_velocity";
  //inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IU = mi[fn];
  //fn = "y_velocity";
  //inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IV = mi[fn];
  //fn = "z_velocity"; 
  //inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IW = mi[fn];
#if (NUMAUX > 0)
  fn = "mixture_fraction_userdef_0";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IMIXF0 = mi[fn];
  fn = "mixture_fraction_userdef_1";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IMIXF1 = mi[fn];
  fn = "age_0";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IAGE0 = mi[fn];
  fn = "age_1";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IAGE1 = mi[fn];
  fn = "agepv_0";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IAGEPV0 = mi[fn];
  fn = "agepv_1";
  inNames.emplace_back(fn); nCompIn = inNames.size(); mi[fn] = nCompIn - 1; const int IAGEPV1 = mi[fn];
#endif

  for (int i = 0; i < nCompIn; i++) {
    destFillComps.emplace_back(i);
  }
  amrex::Print() << "Read " << nCompIn 
    << " fields (inNames[i] -> destFillComps[i]):" << std::endl;
  for (int i=0; i<inNames.size(); ++i) {
    amrex::Print() << "   " << inNames[i] << " -> " << destFillComps[i] << std::endl;
  }
 
  // List of output fields
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
  fn = "HeatReleaseFI";  
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "Y(H2)";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "pv";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "rhorr(NO)";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "rhorr(N2O)";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "rhorr(NNH)";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "FI";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "R10";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
  fn = "zone";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1; 
#if (NUMAUX > 0)
  fn = "agepv_1";
  outNames.emplace_back(fn); nCompOut = outNames.size(); mo[fn] = nCompOut - 1;
#endif

  // List of fields to be conditioned upon (X of <Y|X>)
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

  // List of variables to be conditionally averaged (Y of <Y|X>)
  amrex::Vector<std::string> avgVarNames;
  amrex::Vector<weight_t> avgVarWeights;
  std::map<std::string, int> mav;
  int nAvgVars = 0;
  //int nAvgVars = nCompOut;
  //for (int i = 0; i < nAvgVars; i++) {
  //  avgVarNames.push_back(outNames[i]);
  //  avgVarWeights.push_back(w_volume); 
  //}
  // rho
  fn = "rho"; avgVarNames.emplace_back(fn);
  int ID_rho = avgVarNames.size()-1; mav[fn] = avgVarNames.size()-1; 
  // Z
  fn = "mixture_fraction"; avgVarNames.emplace_back(fn);
  int ID_Z = avgVarNames.size()-1; mav[fn] = avgVarNames.size()-1; 
  // rhoT
  fn = "rhoT"; avgVarNames.emplace_back(fn);
  int ID_rhoT = avgVarNames.size()-1; mav[fn] = avgVarNames.size()-1; 
  fn = "rhoT2"; avgVarNames.emplace_back(fn);
  int ID_rhoT2 = avgVarNames.size()-1; mav[fn] = avgVarNames.size()-1; 
  // HRR
  fn = "HeatRelease"; avgVarNames.emplace_back(fn); 
  int ID_HRR = avgVarNames.size()-1; mav[fn] = avgVarNames.size()-1; 
  fn = "HeatRelease2"; avgVarNames.emplace_back(fn); 
  int ID_HRR2 = avgVarNames.size()-1; mav[fn] = avgVarNames.size()-1; 
  // pv
  fn = "pv"; avgVarNames.emplace_back(fn); 
  int ID_pv = avgVarNames.size()-1; mav[fn] = avgVarNames.size()-1; 
  // rhoY
  int ID_rhoY = ID_pv + 1;
  for (int isp = 0; isp < NUM_SPECIES; isp++) {
    fn = "rhoY(" + spec_names[isp] + ")"; avgVarNames.emplace_back(fn); 
    mav[fn] = avgVarNames.size() - 1;
  }
  // rhoY2
  int ID_rhoY2 = ID_rhoY + NUM_SPECIES;
  for (int isp = 0; isp < NUM_SPECIES; isp++) {
    fn = "rhoY2(" + spec_names[isp] + ")"; avgVarNames.emplace_back(fn); 
    mav[fn] = avgVarNames.size() - 1;
  }
  // wdot
  int ID_wdot = ID_rhoY2 + NUM_SPECIES;
  for (int isp = 0; isp < NUM_SPECIES; isp++) {
    fn = "wdot(" + spec_names[isp] + ")"; avgVarNames.emplace_back(fn); 
    mav[fn] = avgVarNames.size() - 1;
  }
  // wdot2
  int ID_wdot2 = ID_wdot + NUM_SPECIES;
  for (int isp = 0; isp < NUM_SPECIES; isp++) {
    fn = "wdot2(" + spec_names[isp] + ")"; avgVarNames.emplace_back(fn); 
    mav[fn] = avgVarNames.size() - 1;
  }
  fn = "rhorr(NO)"; 
  avgVarNames.emplace_back(fn); mav[fn] = avgVarNames.size() - 1; 
  fn = "rhorr(N2O)"; 
  avgVarNames.emplace_back(fn); mav[fn] = avgVarNames.size() - 1; 
  fn = "rhorr(NNH)"; 
  avgVarNames.emplace_back(fn); mav[fn] = avgVarNames.size() - 1; 
	// progress rate of reaction
  int ID_PRR = ID_wdot2 + NUM_SPECIES + 3;
  for (int ir = 0; ir < NUM_REACTIONS; ir++) {
    std::string s = std::to_string(ir);
    fn = "R(" + std::string(s) + ")"; avgVarNames.emplace_back(fn); 
    mav[fn] = avgVarNames.size() - 1;
  }
#if (NUMAUX > 0)
  // Residence times
	fn = "agepv_1"; avgVarNames.emplace_back(fn);
  int ID_agepv1 = avgVarNames.size()-1; mav[fn] = avgVarNames.size()-1;
#endif

  nAvgVars = avgVarNames.size();
  for (int i = 0; i < nAvgVars; i++) {
    avgVarWeights.emplace_back(w_volume);
  }
  Vector<Real> dataX(nVars);
  Vector<Real> dataY(nAvgVars);

  // Conditional variable fields index and names
  amrex::Vector<std::string> midNames;
  std::map<std::string, int> mm;
  int nCompMid = 0;
  for (int i=0; i<nVars; i++) {
    fn = varNames[i];
    midNames.emplace_back(fn); nCompMid = midNames.size(); mm[fn] = nCompMid - 1; 
  }
  amrex::Print() << "midNames size: " << midNames.size() << std::endl;

  // Progress variable  
  Vector<Real> pv_min = {-0.23290922, -0.24825103, -0.26359285, -0.27893466, -0.29427648,
       -0.3096183 , -0.32496011, -0.34030193, -0.35564374, -0.37098556,
       -0.38632737, -0.40166919, -0.41701101, -0.43235282, -0.44769464,
       -0.46303645, -0.47837827, -0.49372008, -0.5090619 , -0.52440372,
       -0.53974553, -0.55508735, -0.57042916, -0.58577098, -0.60111279,
       -0.61645461, -0.63179642, -0.64713824, -0.66248006, -0.67782187,
       -0.69316369, -0.7085055 , -0.72384732, -0.73918913, -0.75453095,
       -0.76987277, -0.78521458, -0.8005564 , -0.81589821, -0.83124003,
       -0.84658184, -0.86192366, -0.87726547, -0.89260729, -0.90794911,
       -0.92329092, -0.93863274, -0.95397455, -0.96931637, -0.98465818,
       -1.};
  Vector<Real> pv_max = {-0.23290756,  0.11198707,  0.23760585,  0.21397159,  0.18827926,
        0.16245364,  0.13662336,  0.11079683,  0.08497476,  0.05915878,
        0.03335191,  0.00755889, -0.01821313, -0.04395368, -0.06964796,
       -0.09527449, -0.12080348, -0.14619516, -0.1713993 , -0.19635705,
       -0.22100643, -0.24529162, -0.26917403, -0.29264112, -0.31570967,
       -0.33842297, -0.36084441, -0.38305109, -0.40512986, -0.42717719,
       -0.44930301, -0.47163929, -0.494354  , -0.51766916, -0.54187145,
       -0.56726962, -0.59402377, -0.62193754, -0.65056081, -0.67951156,
       -0.70858909, -0.73771137, -0.76684866, -0.7959908 , -0.82513446,
       -0.85427858, -0.88342283, -0.91256711, -0.94171141, -0.9708557 ,
       -1.};
  amrex::Real ztab_min = 0.0;
  amrex::Real ztab_max = 1.0;
  int nztab = pv_max.size();
  amrex::Real dz_tab = (ztab_max - ztab_min) / amrex::Real(nztab - 1);
  Vector<Real> ztab(nztab, 0.0);
  for (int iz = 0; iz < nztab; iz++) {
    ztab[iz] = ztab_min + Real(iz) * dz_tab;
    //amrex::Print() << "iz: " << iz << ", ztabl[" << iz << "] = " << ztab[iz] << std::endl;
  }
  //amrex::Print() << "nztab = " << nztab << ", dz_tab = " << dz_tab << std::endl;

  // Temporary data that can be repetitively used
  int reaction_map[NUM_REACTIONS];
  GET_RMAP(reaction_map);

  // Iterate over input plot files
  for (int iPlot; iPlot < nPlotFiles; ++iPlot) {

    amrex::Print() << "Processing " << iPlot << std::endl;
    AmrData& amrData = dataServicePtrVec[iPlot]->AmrDataRef();

    Vector<Real> probLo = amrData.ProbLo();
    Vector<Real> probHi = amrData.ProbHi();
    amrex::Real Lx = probHi[0] - probLo[0]; 
    amrex::Real Ly = probHi[1] - probLo[1]; 
    amrex::Real Lz = probHi[2] - probLo[2]; 

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
    } // lev - set variables && collect stats

  } // iPlot
} // end of main
