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
      amrex::Print(ioproc) << "   " << plotFileNames[i] << std::endl;
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
  int finestLevel(-1);
  pp.query("finestLevel", finestLevel);
  int nLevels = finestLevel + 1;

  // Variables to be conditioned upon
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

  // Get species names
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);

  // Initialize transport data
  auto eos = pele::physics::PhysicsType::eos();
  pele::physics::transport::TransportParams<
      pele::physics::PhysicsType::transport_type>
      trans_parms;
  trans_parms.allocate();

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

  // Initialize input index and names
  amrex::Vector<std::string> inNames;
  amrex::Vector<int> destFillComps;
  std::map<std::string,int> mi;
  int nCompIn = 0;

  auto declareIn = [&inNames, &destFillComps, &nCompIn](const std::string &fn) {
    inNames.emplace_back(fn);
    destFillComps.emplace_back(nCompIn);
    nCompIn++;
    return 0;
  };
  // 1.species
  for (const auto &spn : spec_names) {
    std::string fn = "Y(" + spn + ")";
    declareIn(fn);
  }
  std::string spName= spec_names[0];
  // 2. temperature
  std::string Tname = "temp";
  declareIn(Tname); mi[Tname] = nCompIn - 1;
  std::string Rname = "density";
  declareIn(Rname); mi[Rname] = nCompIn - 1;
  std::string Uname = "x_velocity";
  declareIn(Uname); mi[Uname] = nCompIn - 1;
  std::string Vname = "y_velocity";
  declareIn(Vname); mi[Vname] = nCompIn - 1;
  std::string Wname = "w_velocity";
  declareIn(Wname); mi[Wname] = nCompIn - 1;
  amrex::Print() << "Read " << nCompIn 
    << " fields (inNames[i] -> destFillComps[i]):" << std::endl;
  for (int i=0; i<inNames.size(); ++i) {
    amrex::Print() << "   " << inNames[i] << " -> " << destFillComps[i] << std::endl;
  }

  // Initialize output index and names
  amrex::Vector<std::string> outNames;
  outNames.emplace_back("rhoU"); 
  outNames.emplace_back("rhoV");
  outNames.emplace_back("rhoW");

  std::map<std::string, int> mo;
  int nCompOut = outNames.size();

  //int nCompIn = 0;

} // end of main