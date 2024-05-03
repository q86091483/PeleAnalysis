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
    amrex::Print() << "iz: " << iz << ", ztabl[" << iz << "] = " << ztab[iz] << std::endl;
  }
  amrex::Print() << "nztab = " << nztab << ", dz_tab = " << dz_tab << std::endl;

  // Temporary data that can be repetitively used

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

      // Output MultiFabs
      mfv_out[lev] = MultiFab(ba,dm,nCompOut,nGrow);
      //MultiFab mf_out(mfv_out[lev], amrex::make_alias, 0, nCompOut);

      // Mid MultiFabs - conditional variables except x, y, z
      MultiFab mf_mid(ba, dm, nCompMid, nGrow);
      MultiFab mf_xyz(ba, dm, 3, nGrow);

      // Intermediate fields MultiFabs
      MultiFab mf_D(ba, dm, NUM_SPECIES, nGrow);
      MultiFab mf_lam(ba, dm, 1, nGrow);
      MultiFab mf_chi(ba, dm, 1, nGrow);
      MultiFab mf_xi(ba, dm, 1, nGrow);
      MultiFab mf_mu(ba, dm, 1, nGrow);
      MultiFab mf_gradu(ba, dm, BL_SPACEDIM+1, nGrow);
      MultiFab mf_gradv(ba, dm, BL_SPACEDIM+1, nGrow);
      MultiFab mf_gradw(ba, dm, BL_SPACEDIM+1, nGrow);
      MultiFab mf_gradYfu(ba, dm, BL_SPACEDIM+1, nGrow);
      MultiFab mf_gradYox(ba, dm, BL_SPACEDIM+1, nGrow);
		  //MultiFab mf_su(ba, dm, 6, nGrow); // symmetric velocity tensor 0.5(du_i/dx_j+du_j/dx_i)
		  //MultiFab mf_tau(ba, dm, 6, nGrow); // viscous force tensor
      //const int L11=0, L22=1, L33=2, L12=3, L13=4, L23=5;
      //const int                      L21=3, L31=4, L32=5;
      //const int LSUM = 6;

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
      //vol /= vol0; // Commented out Apr 30 2024
      if (verbose > 1) {
        std::cout << "   - collecting statistics" << std::endl;
        std::cout << "      - dx, dy, dz [m] =  " << dx[0] << ", " << dx[1] << ", " << dx[2] << "\n"
                  << "      - cell volume  at level 0=  " << vol*vol0 << " [m^3] / "
                  << vol << " (norm)."
                  << std::endl;
      }

      // Jet trajectory
      const amrex::Real coeff_A = 0.55;
      const amrex::Real coeff_B = 0.2;
      const amrex::Real coeff_J = 6.6;
      const amrex::Real Djet = 4.5E-4;
      const amrex::Real slope_zone = 1.5; 
      const amrex::Real ZONE_LEEWARD = 1.5;
      const amrex::Real ZONE_WINDWARD = 2.5;
      const amrex::Real ZONE_INTERACTION = 3.5;

      //FArrayBox alias_fab(orig_fab, amrex::make_alias, 1, 2);
      // Calculate result variables & collect their statistics

      for (amrex::MFIter mfi(mf_in, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();
        countLvl[lev] += bx.volume();

        // FArraybox reference to mf_in
        FArrayBox& fab_in = mf_in[mfi];

        // FArraybox reference mfv_out
        FArrayBox& fab_out = mfv_out[lev][mfi];

        // FArraybox reference mfv_out
        FArrayBox& fab_mid = mf_mid[mfi];

        // FArraybox reference to intermediate
        FArrayBox& fab_xyz = mf_xyz[mfi];
        FArrayBox& fab_covered = mf_covered[mfi];
        FArrayBox& fab_gradu = mf_gradu[mfi];
        FArrayBox& fab_gradv = mf_gradv[mfi];
        FArrayBox& fab_gradw = mf_gradw[mfi];
        FArrayBox& fab_gradYfu = mf_gradYfu[mfi];
        FArrayBox& fab_gradYox = mf_gradYox[mfi];

        //FArrayBox& fab_su    = mf_su[mfi];
        //FArrayBox& fab_tau   = mf_tau[mfi];

        // Array4 reference to mf_in
        Array4<Real> const& rho_a   = mf_in.array(mfi, mi["density"]);
        Array4<Real> const& mixfrac_a  = mf_in.array(mfi, mi["mixture_fraction"]);
        Array4<Real> const& Y_a     = mf_in.array(mfi, mi[spec_names[0]]);
        Array4<Real> const& T_a     = mf_in.array(mfi, mi["temp"]);
        Array4<Real> const& HRR_a   = mf_in.array(mfi, mi["HeatRelease"]);
        Array4<Real> const& U_a     = mf_in.array(mfi, mi["x_velocity"]);
        Array4<Real> const& V_a     = mf_in.array(mfi, mi["y_velocity"]);
        Array4<Real> const& W_a     = mf_in.array(mfi, mi["z_velocity"]);

        // Array reference to mfv_out
        Array4<Real> const& rho_out_a = mfv_out[lev].array(mfi, mo["rho"]);
        Array4<Real> const& mixfrac_out_a = mfv_out[lev].array(mfi, mo["mixture_fraction"]);
        Array4<Real> const& Y_H2_a = mfv_out[lev].array(mfi, mo["Y(H2)"]);
        Array4<Real> const& T_out_a = mfv_out[lev].array(mfi, mo["temp"]);
        Array4<Real> const& HRR_out_a = mfv_out[lev].array(mfi, mo["HeatRelease"]);
        Array4<Real> const& HRRFI_out_a = mfv_out[lev].array(mfi, mo["HeatReleaseFI"]);
        Array4<Real> const& pv_out_a = mfv_out[lev].array(mfi, mo["pv"]);
        Array4<Real> const& rhorr_NO_out_a = mfv_out[lev].array(mfi, mo["rhorr(NO)"]);
        Array4<Real> const& rhorr_N2O_out_a = mfv_out[lev].array(mfi, mo["rhorr(N2O)"]);
        Array4<Real> const& rhorr_NNH_out_a = mfv_out[lev].array(mfi, mo["rhorr(NNH)"]);
        Array4<Real> const& FI_out_a = mfv_out[lev].array(mfi, mo["FI"]);
        Array4<Real> const& R10_out_a = mfv_out[lev].array(mfi, mo["R10"]);
        Array4<Real> const& zone_out_a = mfv_out[lev].array(mfi, mo["zone"]);

        //Array4<Real> const& mu_out_a = mfv_out[lev].array(mfi, mo["mu"]);
        //Array4<Real> const& ts_a      = mfv_out[lev].array(mfi, mo["ts11"]);

        // Array reference to mf_mid
        Array4<Real> const& mixfrac_mid_a = mf_mid.array(mfi, mm["mixture_fraction"]);
        Array4<Real> const& pv_mid_a = mf_mid.array(mfi, mm["pv"]);
        Array4<Real> const& zone_mid_a = mf_mid.array(mfi, mm["zone"]);
        Array4<Real> const& FI_mid_a = mf_mid.array(mfi, mm["FI"]);

        Array4<Real> const& rhowdot_mid_a = mf_mid.array(mfi, mm["pv"]);

        // Array reference intermediate
        Array4<Real> const& x_a     = mf_xyz.array(mfi, 0);
        Array4<Real> const& y_a     = mf_xyz.array(mfi, 1);
        Array4<Real> const& z_a     = mf_xyz.array(mfi, 2);  
        Array4<Real> const& covered_a = mf_covered.array(mfi); 
        Array4<Real> const& D_a     = mf_D.array(mfi);
        Array4<Real> const& lam_a   = mf_lam.array(mfi);
        Array4<Real> const& chi_a   = mf_chi.array(mfi);
        Array4<Real> const& xi_a    = mf_xi.array(mfi);
        Array4<Real> const& mu_a    = mf_mu.array(mfi);
        //Array4<Real> const& su_a    = mf_su.array(mfi);
        //Array4<Real> const& tau_a   = mf_tau.array(mfi);
        Array4<Real> const& gradu_a = mf_gradu.array(mfi);
        Array4<Real> const& gradv_a = mf_gradv.array(mfi);
        Array4<Real> const& gradw_a = mf_gradw.array(mfi);
        Array4<Real> const& gradYfu_a = mf_gradYfu.array(mfi);
        Array4<Real> const& gradYox_a = mf_gradYox.array(mfi);

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
        gradient(BL_TO_FORTRAN_BOX(bx),
                BL_TO_FORTRAN_N_ANYD(fab_in,    mi["Y(H2)"]),
                BL_TO_FORTRAN_N_ANYD(fab_gradYfu, 0),
                &(dx[0]));
        gradient(BL_TO_FORTRAN_BOX(bx),
                BL_TO_FORTRAN_N_ANYD(fab_in,    mi["Y(O2)"]),
                BL_TO_FORTRAN_N_ANYD(fab_gradYox, 0),
                &(dx[0]));

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          // var_loc(i,j,k)
          amrex::Real wdot_loc[NUM_SPECIES] = {0.0_rt};
          amrex::Real Ci_MKS[NUM_SPECIES] = {0.0_rt};
          amrex::Real Ci_CGS[NUM_SPECIES] = {0.0_rt};
          amrex::Real Y_loc[NUM_SPECIES] = {0.0_rt};
          amrex::Real X_loc[NUM_SPECIES] = {0.0_rt};
          amrex::Real Qf[NUM_REACTIONS] = {0.0_rt};
          amrex::Real Qr[NUM_REACTIONS] = {0.0_rt};
          amrex::Real rho_loc     = 0.0_rt;
          amrex::Real rho_cgs     = 0.0_rt;
          amrex::Real Pcgs        = 0.0_rt;
          amrex::Real T_loc       = 0.0_rt;
          int reaction_map[NUM_REACTIONS];
          GET_RMAP(reaction_map);

          rho_loc = rho_a(i,j,k);
          rho_cgs = rho_loc * 0.001_rt; // kg/m3 to g/cm3
          T_loc = T_a(i,j,k);
          for (int isp = 0; isp < NUM_SPECIES; isp++) {
            Y_loc[isp] = Y_a(i,j,k,isp);
          }

          // mfv_out
          rho_out_a(i,j,k) = rho_a(i,j,k);
          T_out_a(i,j,k) = T_a(i,j,k);
          Y_H2_a(i,j,k) = Y_a(i,j,k,H2_ID);
          // mfv_out - mixture fraction
          amrex::Real Zlocal = 0.0;
          for (int isp = 0; isp < NUM_SPECIES; ++isp) {
            Zlocal += spec_Bilger_fact[isp] * Y_a(i,j,k,isp);
          }
          Zlocal = (Zlocal-Zox) / (Zfu-Zox);
          mixfrac_out_a(i,j,k) = Zlocal;
          // mfv_out - progress variable
          amrex::Real pv0, pv1;
          int iztab0 = 0;
          iztab0 = floor((Zlocal-ztab_min)/(ztab_max-ztab_min)*nztab);
          iztab0 = std::max(std::min(iztab0, nztab - 1), 0);
          if (iztab0 >= (nztab-1)) iztab0 = iztab0 - 1;
          amrex::Real r0 = 1 - (Zlocal - ztab[iztab0]) / dz_tab;
          pv0 = pv_min[iztab0]*r0 + pv_min[iztab0+1]*(1-r0); 
          pv1 = pv_max[iztab0]*r0 + pv_max[iztab0+1]*(1-r0); 
          pv_out_a(i,j,k) = -1.0 * Y_a(i,j,k,H2_ID) + 
                            -1.0 * Y_a(i,j,k,O2_ID) +
                            +1.0 * Y_a(i,j,k,H2O_ID); 
          pv_out_a(i,j,k) = (pv_out_a(i,j,k)-pv0)/(pv1-pv0);
          if (Zlocal < 1E-4) pv_out_a(i,j,k) = 0.0;
          // mfv_out - wdot
          eos.RTY2WDOT(rho_cgs, T_loc, Y_loc, wdot_loc);  // g/cm3
          //for (int n = 0; n < NUM_SPECIES; n++) {
          //  rr_a(i,j,k,n) = wdot[n] * 1000.0_rt; // rhodot, CGS -> MKS conversion
          //}
          rhorr_NO_out_a(i,j,k) = wdot_loc[NO_ID] * 1000.0_rt; // kg/m3
          rhorr_N2O_out_a(i,j,k) = wdot_loc[N2O_ID] * 1000.0_rt;
          rhorr_NNH_out_a(i,j,k) = wdot_loc[NNH_ID] * 1000.0_rt;
          // net rate of reaction progress
          eos.Y2X(Y_loc, X_loc);
          CKPX(rho_cgs, T_a(i,j,k), X_loc, Pcgs);        
          CKYTCR(rho_cgs, T_a(i,j,k), Y_loc, Ci_CGS);
          for (int isp = 0; isp < NUM_SPECIES; isp++) {
            Ci_MKS[isp] = Ci_CGS[isp]*1.0e6_rt;                             // CGS -> MKS conversion
          }
          //CKKFKR(Pcgs, T_a(i,j,k), X_loc, Qf, Qr);
          progressRateFR(Qf, Qr, Ci_MKS, T_a(i,j,k));
          R10_out_a(i,j,k) = Qf[0]; //(Qf[reaction_map[10]] - Qr[reaction_map[10]]) * 1.0E6; // - Qr[9];
          // mfv_out - FI 
          FI_out_a(i,j,k) = 0.0;
          FI_out_a(i,j,k) = gradYfu_a(i,j,k,0) * gradYox_a(i,j,k,0) +
                            gradYfu_a(i,j,k,1) * gradYox_a(i,j,k,1) +
                            gradYfu_a(i,j,k,2) * gradYox_a(i,j,k,2);
          FI_out_a(i,j,k) = FI_out_a(i,j,k) / std::fabs(FI_out_a(i,j,k));
          //if (HRR_a(i,j,k) < 1E3) {
          //  FI_out_a(i,j,k) = 0.0;
          //}
          HRR_out_a(i,j,k) = HRR_a(i,j,k);
          if (FI_out_a(i,j,k) > 0.0) {
            HRRFI_out_a(i,j,k) = HRR_a(i,j,k);
          } else if (FI_out_a(i,j,k) < 0.0) {
            HRRFI_out_a(i,j,k) = -HRR_a(i,j,k);
          }
          // zone
          //ys = J * Djet * A * np.power(xs / (J*Djet), B)
          amrex::Real yp = 0.0_rt; 
          amrex::Real hj = coeff_J * Djet * coeff_A * std::pow(x_a(i,j,k)/(coeff_J*Djet), coeff_B); 
          amrex::Real bl0, bl1, br0, br1, zl0, zl1, zr0, zr1, zp;
          // Use periodicity along y direction
          yp = probLo[1] + std::fmod(y_a(i,j,k)-probLo[1], Ly/2.);
          bl0 = hj - (-slope_zone) * (-2.0*Djet); // Lean towards left or right
          bl1 = hj - (-slope_zone) * (+2.0*Djet);
          br0 = hj - (+slope_zone) * (-2.0*Djet);
          br1 = hj - (+slope_zone) * (+2.0*Djet);
          zl0 = (-slope_zone) * yp + bl0;
          zl1 = (-slope_zone) * yp + bl1;
          zr0 = (+slope_zone) * yp + br0;
          zr1 = (+slope_zone) * yp + br1;
          zp = z_a(i,j,k); 
          zone_out_a(i,j,k) = -1.0;
          if ((zp>zr0) and (zp>zl0)) {
            zone_out_a(i,j,k) = ZONE_WINDWARD;
          } else if ((zp<zl0) and (zp<zr0)) {
            zone_out_a(i,j,k) = ZONE_LEEWARD;
          } else {
            zone_out_a(i,j,k) = ZONE_INTERACTION;
          }
          if (x_a(i,j,k) < 0.0) {
            zone_out_a(i,j,k) = -1.0; // Disregard inlet regions
          }

          // mf_mid
          mixfrac_mid_a(i,j,k) = mixfrac_out_a(i,j,k);
          pv_mid_a(i,j,k) = pv_out_a(i,j,k);
          zone_mid_a(i,j,k) = zone_out_a(i,j,k);
          FI_mid_a(i,j,k) = FI_out_a(i,j,k);

        });

        // Collect statistics for mfi
        amrex::Real rho_loc = 0.0_rt;
        amrex::Real rho_cgs = 0.0_rt;
        amrex::Real T_loc = 0.0_rt;
        amrex::Real Y_loc[NUM_SPECIES] = {0.0_rt};
        amrex::Real wdot_loc[NUM_SPECIES] = {0.0_rt};
        for (IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
          // Index
          int i = iv[0];
          int j = iv[1];
          int k = iv[2];
          // rho, T, Y
          rho_loc = rho_a(i,j,k);
          rho_cgs = rho_loc * 0.001_rt; // kg/m3 to g/cm3
          T_loc   = T_a(i,j,k);
          for (int isp = 0; isp < NUM_SPECIES; isp++) {
            Y_loc[isp] = Y_a(i,j,k,isp);
            wdot_loc[isp] = 0.0;
          }
          // wdot
          eos.RTY2WDOT(rho_cgs, T_loc, Y_loc, wdot_loc);  // g/cm3
          for (int isp = 0; isp < NUM_SPECIES; isp++) {
            wdot_loc[isp] = wdot_loc[isp] * 1000.0_rt; // kg/m3
          }
          
          // Get dataX such as x, y, z and variables in mf_mid
          for (int ivar = 0; ivar < nVars; ivar++) {
            std::string vn = varNames[ivar]; 
            if (vn == "x") {
              dataX[ivar] = mf_xyz.array(mfi)(i, j, k, 0);
            } else if (vn == "y") {
              dataX[ivar] = mf_xyz.array(mfi)(i, j, k, 1);
            } else if (vn == "z") {
              dataX[ivar] = mf_xyz.array(mfi)(i, j, k, 2); 
            } else {
              dataX[ivar] = mf_mid.array(mfi)(i, j, k, mm[vn]);
            } 
          }
          //dataX[0] = mf_xyz.array(mfi)(iv[0], iv[1], iv[2], 0);
          //dataX[1] = mf_xyz.array(mfi)(iv[0], iv[1], iv[2], 1);
          //dataX[2] = mf_xyz.array(mfi)(iv[0], iv[1], iv[2], 2);

          //Skip points overlapping with a finer levels
          skip = false;
          if (covered_a(i, j, k) > 0.0) skip = true;
          for (int ivar=0; ivar<nVars; ivar++) {
            if ((incOutOfBounds[ivar] == 0) &&
                ((dataX[ivar] < varBounds[ivar][0]) ||
                 (dataX[ivar] >= varBounds[ivar][1]))) {
              skip = true;
            }
          }
          if (skip) continue;

          // Get dataY
          //fab_out.getVal(dataY.dataPtr(), iv); // original way
          dataY[ID_rho]   = rho_a(i,j,k);
          dataY[ID_Z]     = mixfrac_out_a(i,j,k);
          dataY[ID_rhoT]  = rho_a(i,j,k) * T_a(i,j,k); 
          dataY[ID_rhoT2] = rho_a(i,j,k) * T_a(i,j,k) * T_a(i,j,k);
          dataY[ID_HRR]   = HRR_a(i,j,k);
          dataY[ID_HRR2]  = HRR_a(i,j,k) * HRR_a(i,j,k);
          for (int isp = 0; isp < NUM_SPECIES; isp++) {
            dataY[ID_rhoY + isp] = rho_a(i,j,k) * Y_a(i,j,k,isp);
            dataY[ID_rhoY2 + isp] = rho_a(i,j,k) * Y_a(i,j,k,isp) * Y_a(i,j,k,isp);
            dataY[ID_wdot + isp] = rho_a(i,j,k) * wdot_loc[isp];
            dataY[ID_wdot2 + isp] = rho_a(i,j,k) * wdot_loc[isp] * wdot_loc[isp];
          }
          dataY[ID_pv] = pv_out_a(i,j,k);

          dataY[mav["rhorr(NO)"]] = rhorr_NO_out_a(i,j,k);
          dataY[mav["rhorr(N2O)"]] = rhorr_N2O_out_a(i,j,k);
          dataY[mav["rhorr(NNH)"]] = rhorr_NNH_out_a(i,j,k);

          // Compute the bin in X space (<Y|X>)
          for (int ivar=0; ivar<nVars; ivar++) {
            bins[ivar] = computeBin(dataX[ivar], varBounds[ivar][0], varBounds[ivar][1], nBins[ivar], binType[ivar]);
          }
          bin = 0;
          for (int ivar=0; ivar<nVars-1; ivar++) {
            bin += bins[ivar];
            bin *= nBins[ivar+1];
          }
          bin += bins[nVars-1];

          // Gather statistics
          rho = dataY[mav["density"]];
          m = rho * vol;
          count[bin] += 1;
          volMean[bin] += vol;
          volStd[bin] += vol * vol;
          rhoMean[bin] += rho;
          rhoStd[bin] += rho * rho;
          massMean[bin] += m;
          massStd[bin] += m * m;
          for (int ivar=0; ivar<nAvgVars; ivar++) {
            varMeanVal[ivar][bin] += dataY[ivar] * vol;
            varStdVal[ivar][bin] += dataY[ivar] * dataY[ivar] * vol;
            varMinVal[ivar][bin] += std::min(varMinVal[ivar][bin], dataY[ivar]);
            varMaxVal[ivar][bin] += std::min(varMaxVal[ivar][bin], dataY[ivar]);
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
    if (writeDerivedField) {
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
          GetVecOfConstPtrs(mfv_out), outNames, geoms, amrData.Time(), istep, ref_ratio);
    }
  } // iPlot
} // end of main
