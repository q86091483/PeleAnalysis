#include <string>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>

#include <H5Cpp.h>

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

using namespace amrex;

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
}

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
}

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

static void print_usage (int argc, char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0]
        << "\n"
        << "   help                            : Print this help\n"
        << "   verbose = 1                     : Verbosity level (default = 0)\n"
        << "   finestLevel = 2                 : Finest level to include in the statistics (default = -1)\n"
        << "   infile = plt1 plt2              : List of plot files\n"
        << "   outputDir = cond_stats          : Name of the output directory for the result files (default = cond_stats)\n"
        << "   writeGrid = 0                   : Write uniform grid properties (e.g. for XDMF/VTK) to the result files (default = 0)\n"
        << "   writeCoordinates = 1            : Write coordinates (i.e. the bins) to the result files (default = 1)\n"
        << "   vars = temp density             : List of conditioning variables\n"
        << "   nBins = 11 12                   : List of bins per conditioning variables\n"
        << "   varBounds = 600 1200 15.0 24.0  : List of bounds for the conditioning variable bins (2x nVars)\n"
        << "   incOutOfBounds = 0 1            : List of switches to indicate whether values outside the specified bounds should be included in the lower/upper bins (default = 0)\n"
        << "   binType = lin log               : List of types to indicate whether the bins are equally space in linear- or log-space (default = lin)\n"
        << "   avgVars = HeatRelease Y(N2)     : List of variables for statistics\n"
        << "   avgVarWeights = volume density  : List of weights for statistics\n"
        << std::endl;
    exit(1);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    ParmParse pp;
    const bool isioproc = ParallelDescriptor::IOProcessor();
    const int ioproc = ParallelDescriptor::IOProcessorNumber();

    if (argc < 2 && isioproc) {
        print_usage(argc, argv);
    }

    if (pp.contains("help")) {
        print_usage(argc, argv);
    }

    int verbose = 0;
    pp.query("verbose", verbose);
    if (verbose > 2) {
        AmrData::SetVerbose(true);
    }
    if (! isioproc) {
        verbose = 0;
    }

    const int nPlotFiles = pp.countval("infile");
    if (nPlotFiles < 1) {
        Print(ioproc) << "Bad nPlotFiles, exiting ..." << std::endl; 
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    const int nPlotFilesDerived = pp.countval("infile_derived");
    if (nPlotFilesDerived < 1) {
        Print(ioproc) << "No derived field is needed." << std::endl; 
    }
    if (verbose > 0) {
        Print(ioproc) << "Processing " << nPlotFiles + nPlotFilesDerived << " plotfiles ..." << std::endl;
    }

    // Settings of IO
    std::string outDir = "cond_ISRN";
    pp.query("outputDir", outDir);
    UtilCreateDirectory(outDir, 0755);
    bool writeGrid = false;
    pp.query("writeGrid", writeGrid);
    bool writeCoordinates = false;
    pp.query("writeCoordinates", writeCoordinates);
    std::string outputLabel;
    pp.query("outputLabel", outputLabel);

    // Plot file names
    amrex::Vector<std::string> plotFileNames(nPlotFiles);
    for (int i = 0; i < nPlotFiles; ++i) {
        pp.get("infile", plotFileNames[i], i);
        if (verbose > 1) {
            amrex::Print(ioproc) << "    " << plotFileNames[i] << std::endl;
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
            std::cout << "  " << varNames[i]
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

    // Variables to be average
    amrex::Vector<std::string> avgVarNames;
    amrex::Vector<weight_t> avgVarWeights;
    for (int i = 0; i < pp.countval("avgVars"); i++) {
        std::string tmp;
        pp.get("avgVars", tmp, i);
        if (tmp == "density") {
            continue;
        }
        avgVarNames.push_back(tmp);
        pp.get("avgVarWeights", tmp, i);
        if (tmp == "none") {
            avgVarWeights.push_back(w_none);
        } else if (tmp == "volume") {
            avgVarWeights.push_back(w_volume);
        } else if (tmp == "density") {
            avgVarWeights.push_back(w_density);
        } else if (tmp == "mass") {
            avgVarWeights.push_back(w_mass);
        } else {
            std::string msg = "Unknown weight factor: " + tmp;
            Error(msg.c_str());
        }
    }
    const int nAvgVars(avgVarNames.size());
    if (verbose > 0) {
        std::cout << "Variable to be averaged:" << std::endl;
        for (int i=0; i<nAvgVars; i++) {
            std::cout << "  " << avgVarNames[i]
                      << " (weight = " << avgVarWeights[i]
                      << ")" << std::endl;
        }
    }
    // All variables to be loaded
    Vector<std::string> destFillNames;
    for (int i=0; i<nVars; i++) {
        destFillNames.push_back(varNames[i]);
    }
    for (int i=0; i<nAvgVars; i++) {
        destFillNames.push_back(avgVarNames[i]);
    }
    destFillNames.push_back("density");
    const int nDest = destFillNames.size();
    Vector<int> destFills(nDest);
    std::iota(destFills.begin(), destFills.end(), 0);
    Vector<Real> dataIV(nDest + 1);
    if (verbose > 0) {
        std::cout << "Variables to be loaded:" << std::endl;
        for (int i=0; i<nDest; i++) {
            std::cout << "  " << destFillNames[i]
                      << " (id = " << destFills[i]
                      << ")" << std::endl;
        }
    }
    int idx = -1, idy = -1, idz = -1;
    for (int iDest = 0; iDest < nDest; ++iDest) {
        const std::string& fn = destFillNames[iDest];
        if (fn == "x") idx = iDest;
        if (fn == "y") idy = iDest;
        if (fn == "z") idz = iDest;
    }

    // Loop over input files
    for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {

        // Open file and get an amrData pointer
        const std::string& infile = plotFileNames[iPlot];
        if (verbose > 0) 
            std::cout << "\nOpening " << infile << " ..." << std::endl;
        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(infile, fileType);
        if (!dataServices.AmrDataOk())
            DataServices::Dispatch(DataServices::ExitRequest, NULL);
        AmrData& amrData = dataServices.AmrDataRef();
        if (verbose > 0) 
            std::cout << "Loading " << infile << " done." << std::endl;

        // If needed, open file for derived fields and get an amrData pointer
        std::string infile_derived;
        if (verbose > 0) 
            std::cout << "\nOpening " << infile_derived << " ..." << std::endl;
        infile_derived = plotFileNamesDerived[iPlot];
        DataServices dataServicesDerived(infile_derived, fileType);
        if (!dataServicesDerived.AmrDataOk()) {
            DataServices::Dispatch(DataServices::ExitRequest, NULL);
        }
        if (verbose > 0) 
            std::cout << "Loading " << infile_derived << " done." << std::endl;
        AmrData& amrData_derived = dataServicesDerived.AmrDataRef();

        // Field names in plot files
        const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
        amrex::Print() << "Fields in the original plot files: " << std::endl;
        for (auto it = plotVarNames.begin(); it < plotVarNames.end(); ++it) {
            amrex::Print() << *it << std::endl;
        }

        // Derived field names 
        amrex::Vector<std::string> plotVarNamesDerived; 
        if (nPlotFilesDerived >= 1) {
            amrex::Print() << "Derived fields plot files: " << std::endl;
            const Vector<std::string> ptn = amrData_derived.PlotVarNames();
            for (auto it = ptn.begin(); it < ptn.end(); ++it) {
                plotVarNamesDerived.push_back(*it);
                amrex::Print() << *it << std::endl;
            }
        }

        // Check the names of the variables are present in the plotfile
        amrex::Vector<int> iPlotGroup(nDest, -1);
        for (int i=0; i<nDest; i++) {
            const std::string& fn = destFillNames[i];
            if (fn== "x" || fn=="y" || fn=="z") {
                iPlotGroup[i] = 0;
            } else if (std::find(plotVarNames.begin(), plotVarNames.end(), fn) != plotVarNames.end()) {
                iPlotGroup[i] = 0;
            } else if (std::find(plotVarNamesDerived.begin(), plotVarNamesDerived.end(), fn) != plotVarNamesDerived.end()){
                iPlotGroup[i] = 1;
            } else {
                std::string msg = "Bad variable name: " + destFillNames[i];
                Error(msg.c_str());
            }
        }

        // Get finest level
        if (finestLevel < 0)
            finestLevel = 0; //amrData.FinestLevel();
            nLevels = finestLevel+1;
        if (verbose > 0)
            amrex::Print(ioproc) << "  Finest level: " << finestLevel << std::endl;
        
        // Get domain size
        Vector<Real> probLo = amrData.ProbLo();
        Vector<Real> probHi = amrData.ProbHi();
        Real domainVol = 1;
        for (int i=0; i<BL_SPACEDIM; i++) {
            domainVol *= probHi[i] - probLo[i];
        }
        if (verbose > 0)
            amrex::Print(ioproc) << "  Domain size:\n"
                      << "    x0, y0, z0 [m] =  " << probLo[0] << ", " << probLo[1] << ", " << probLo[2] << "\n"
                      << "    x1, y1, z1 [m] =  " << probHi[0] << ", " << probHi[1] << ", " << probHi[2] << "\n"
                      << "  Domain volume [m^3] = " << domainVol 
                      << std::endl;
                      // Get min/max for each component
        if (verbose > 1)
            std::cout << "  Global min/max values:" << std::endl;
        Vector<Box> probDomain = amrData.ProbDomain();
        Vector<Box> probDomain_derived = amrData_derived.ProbDomain();
        Vector<Real> varGMin(nVars, 1e30);
        Vector<Real> varGMax(nVars, -1e30);
        for (int i=0; i<nVars; i++) {
            const std::string& fn = destFillNames[i];
            if (fn=="x" || fn=="y" || fn=="z") {
                continue;
            } else if (iPlotGroup[i] == 0) {
                Real mn, mx;
                for (int iLevel=0; iLevel<nLevels; iLevel++) {
                    amrData.MinMax(probDomain[iLevel], varNames[i], iLevel, mn, mx);
                    varGMin[i] = std::min(varGMin[i], mn);
                    varGMax[i] = std::max(varGMax[i], mx);
                }
            } else if (iPlotGroup[i] == 1) {
                Real mn, mx;
                for (int iLevel=0; iLevel<nLevels; iLevel++) {
                    amrData_derived.MinMax(probDomain_derived[iLevel], varNames[i], iLevel, mn, mx);
                    varGMin[i] = std::min(varGMin[i], mn);
                    varGMax[i] = std::max(varGMax[i], mx);
                }
            }
        }
        for (int i=0; i<nDest; i++) {
            Real minVal = 1e30;
            Real maxVal = -1e30;
            Real mn, mx;
            const std::string& fn = destFillNames[i];
            if (fn=="x" || fn=="y" || fn=="z") {
                continue;
            } else if (iPlotGroup[i] == 0) {
                for (int iLevel=0; iLevel<nLevels; iLevel++) {
                    amrData.MinMax(probDomain[iLevel], destFillNames[i], iLevel, mn, mx);
                    minVal = std::min(minVal, mn);
                    maxVal = std::max(maxVal, mx);
                }
            } else if (iPlotGroup[i] == 1) {
                for (int iLevel=0; iLevel<nLevels; iLevel++) {
                    // amrData_derived.MinMax(probDomain_derived[iLevel], destFillNames[i], iLevel, mn, mx);
                    // minVal = std::min(minVal, mn);
                    // maxVal = std::max(maxVal, mx);
                }
            }
            if (verbose > 1)
                std::cout << "  " << destFillNames[i] << " min/max: " << minVal << "/" << maxVal << std::endl;
        }

        const int ngrow     = 0;
        int ratio           = 1;
        int bin             = 0;
        bool skip           = false;
        int ii;
        Real vol, rho, m, weight;
        Real vol0           = 1;

        // Initialise arrays for conditional data
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
        Vector<int> countLvl(nLevels, 0);
        Vector<Real> volLvl(nLevels, 0);

        for (int iLevel=0; iLevel<nLevels; iLevel++) {
            // Cell size for the current level
            const Vector<Real> dx = amrData.CellSize(iLevel);
            // Load data into MultiFab
            if (verbose > 1)
                std::cout << "  Loading data for level " << iLevel << " ... ";
            DistributionMapping dm(amrData.boxArray(iLevel));
            MultiFab mf(amrData.boxArray(iLevel), dm, nDest+1, ngrow);
            // MultiFab mf_derived(amrData_derived.boxArray(iLevel), dm, nDest+1, ngrow);
            // amrData.FillVar(mf, iLevel, destFillNames, destFills);
            for (int iDest = 0; iDest < nDest; iDest++) {
                const std::string& fn = destFillNames[iDest];
                if (fn=="x" || fn=="y" || fn=="z") {
                    amrex::Print() << probLo[0] << std::endl;
                    // Set coordinates x, y, z
                    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                        const Box& box = mfi.validbox();
                        FArrayBox& fab = mf[mfi];
                        Array4<Real> const& phi = fab.array();
                        IntVect iv_lo = box.smallEnd();
                        Real x_lo = probLo[0] + iv_lo[0] * dx[0];
                        Real y_lo = probLo[1] + iv_lo[1] * dx[1];
                        Real z_lo = probLo[2] + iv_lo[2] * dx[2];
                        const auto lo = lbound(box);
                        const auto hi = ubound(box);
                        for (int k = lo.z; k <= hi.z; ++k) {
                            for (int j = lo.y; j <= hi.y; ++j) {
                                for (int i = lo.x; i <= hi.x; ++i) {
                                    if (idx > 0) phi(i,j,k,idx) = probLo[0] + i * dx[0];
                                    if (idy > 0) phi(i,j,k,idy) = probLo[1] + j * dx[1];
                                    if (idz > 0) phi(i,j,k,idz) = probLo[2] + k * dx[2];
                                }
                            }
                        }
                    } // mf - set x y z
                } else if (iPlotGroup[iDest] == 0) {
                    amrData.FillVar(mf, iLevel, fn, iDest);
                } else if (iPlotGroup[iDest] == 1) {
                    amrData_derived.FillVar(mf, iLevel, fn, iDest);
                }
            } // iDest
            mf.setVal(1, nDest, 1);
            if (verbose > 1)
                std::cout << "done." << std::endl;
            
            // Zero out overlapping data
            if (iLevel < finestLevel)
            {
                ratio = amrData.RefRatio()[iLevel];
                BoxArray baf = BoxArray(amrData.boxArray(iLevel+1)).coarsen(ratio);
                for (MFIter mfi(mf); mfi.isValid(); ++mfi)
                {
                    const Box& box = mfi.validbox();
                    FArrayBox& fab = mf[mfi];
                    std::vector<std::pair<int,Box>> isects = baf.intersections(box);
                    for (int ii = 0; ii < isects.size(); ii++)
                        fab.setVal(0,isects[ii].second,nDest,1);
                }
            }

            // Volume for the current level
            vol = 1;
            for (int i = 0; i < BL_SPACEDIM; i++)
                vol *= dx[i];
            // Normalise cell volume for numerical accuracy
            if (iLevel == 0)
                vol0 = vol;
            vol /= vol0;
            if (verbose > 1) {
                std::cout << "    dx, dy, dz [m] =  " << dx[0] << ", " << dx[1] << ", " << dx[2] << "\n"
                          << "    cell volume =  " << vol*vol0 << " [m^3] / "
                          << vol << " (norm)."
                          << std::endl;
            }

            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                FArrayBox& fab = mf[mfi];
                const Box& box = mfi.validbox();

                // Cell count for each level
                countLvl[iLevel] += box.volume();

                for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv)) {
                    fab.getVal(dataIV.dataPtr(), iv);

                    // Skip points overlapping with a finer levels
                    if (dataIV[nDest] == 0) {
                        continue;
                    }

                    // Skip points outside the specified bounds
                    skip = false;
                    for (int i=0; i<nVars-1; i++) {
                        if ((incOutOfBounds[i] == 0) &&
                            ((dataIV[i] <  varBounds[i][0]) ||
                             (dataIV[i] >= varBounds[i][1]))) {
                            skip = true;
                        }
                    }
                    if (skip) {
                        continue;
                    }

                    // Compute the bin
                    for (int i=0; i<nVars; i++) {
                        bins[i] = computeBin(dataIV[i], varBounds[i][0], varBounds[i][1], nBins[i], binType[i]);
                    }
                    bin = 0;
                    for (int i=0; i<nVars-1; i++) {
                        bin += bins[i];
                        bin *= nBins[i+1];
                    }
                    bin += bins[nVars-1];

                    // Gather statistics
                    rho = dataIV[nDest-1];
                    m = rho*vol;
                    count[bin] += 1;
                    volMean[bin] += vol;
                    volStd[bin] += vol*vol;
                    rhoMean[bin] += rho;
                    rhoStd[bin] += rho*rho;
                    massMean[bin] += m;
                    massStd[bin] += m*m;

                    for (int i=0; i<nAvgVars; i++) {
                        switch (avgVarWeights[i])
                        {
                            case w_none:
                            {
                                weight = 1.0;
                                break;
                            }
                            case w_volume:
                            {
                                weight = vol;
                                break;
                            }
                            case w_density:
                            {
                                weight = rho;
                                break;
                            }
                            case w_mass:
                            {
                                weight = m;
                                break;
                            }
                        }
                        ii = nVars + i;
                        varMeanVal[i][bin] += dataIV[ii]*weight;
                        varStdVal[i][bin] += dataIV[ii]*dataIV[ii]*weight;
                        varMinVal[i][bin] = std::min(varMinVal[i][bin], dataIV[ii]);
                        varMaxVal[i][bin] = std::max(varMaxVal[i][bin], dataIV[ii]);
                    } // i - nAvgVars
                } // iv - box
            } // mfi - mf

            // Cell volume for each level
            volLvl[iLevel] = countLvl[iLevel]*vol;
        } // iLevel - nLevels

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

        // Output result
        if (isioproc) {
            std::string h5fileName = outDir + "/" + basename(infile) + "_" + outputLabel + ".h5";
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
                    for (int i=0; i<nVars; i++) {
                        DIMENSIONS[nVars-(i+1)] = nBins[i] + 1;
                    }

                    hsize_t dims[1];
                    dims[0] = nVars;
                    H5::DataSpace dataspace(1, dims);

                    H5::DataSet dataset(h5grid.createDataSet("DIMENSIONS", H5T_UINT, dataspace));
                    dataset.write(DIMENSIONS, H5T_UINT);
                }

                // Write ORIGIN
                {
                    Real ORIGIN[nVars];
                    for (int i=0; i<nVars; i++) {
                        ORIGIN[nVars-(i+1)] = varBounds[i][0];
                    }

                    hsize_t dims[1];
                    dims[0] = nVars;
                    H5::DataSpace dataspace(1, dims);

                    H5::DataSet dataset(h5grid.createDataSet("ORIGIN", H5T_REAL, dataspace));
                    dataset.write(ORIGIN, H5T_REAL);
                }

                // Write SPACING
                {
                    Real SPACING[nVars];
                    for (int i=0; i<nVars; i++) {
                        SPACING[nVars-(i+1)] = (varBounds[i][1] - varBounds[i][0])/std::max(nBins[i], 1);
                    }

                    hsize_t dims[1];
                    dims[0] = nVars;
                    H5::DataSpace dataspace(1, dims);

                    H5::DataSet dataset(h5grid.createDataSet("SPACING", H5T_REAL, dataspace));
                    dataset.write(SPACING, H5T_REAL);
                }
                h5grid.close();
            } // Write grid

            if (writeCoordinates) {
                H5::Group h5coord(h5file.createGroup("COORDS"));
                for (int i=0; i<nVars; i++) {
                    uint d = nBins[i] + 1;
                    Real o = varBounds[i][0];
                    Real coord[d];
                    switch (binType[i])
                    {
                        case b_lin:
                        {
                            Real s = (varBounds[i][1] - varBounds[i][0])/std::max(nBins[i], 1);
                            for (int j = 0; j < d; j++)
                            {
                                coord[j] = Real(j)*s + o;
                            }
                            break;
                        }
                        case b_log:
                        {
                            Real s = std::log10(varBounds[i][1]/varBounds[i][0])/std::max(nBins[i], 1);
                            for (int j = 0; j < d; j++)
                            {
                                coord[j] = std::pow(10, Real(j)*s + std::log10(o));
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

                    attribute = dataset.createAttribute("incOutOfBounds", H5T_INT, H5S_SCALAR);
                    attribute.write(H5T_INT, &incOutOfBounds[i]);

                    attribute = dataset.createAttribute("min", H5T_REAL, H5S_SCALAR);
                    attribute.write(H5T_REAL, &varGMin[i]);

                    attribute = dataset.createAttribute("max", H5T_REAL, H5S_SCALAR);
                    attribute.write(H5T_REAL, &varGMax[i]);
                }
                h5coord.close();
            } // Write coordinate

            // Global (unconditional) data
            H5::Group h5domain(h5file.createGroup("DOMAIN"));
            // Domain size
            writeRealData(h5domain, "prob_lo", Vector<int>(1, BL_SPACEDIM), probLo, 0.0);
            writeRealData(h5domain, "prob_hi", Vector<int>(1, BL_SPACEDIM), probHi, 0.0);
            writeRealData(h5domain, "volume", Vector<int>(1, 1), Vector<Real>(1, domainVol), 0.0);

            // Count/ volume per level (unconditional)
            writeIntData(h5domain, "count_lvl", Vector<int>(1, nLevels), countLvl, 0.0);
            for (auto& v : volLvl) v *= vol0;
            writeRealData(h5domain, "volume_lvl", Vector<int>(1, nLevels), volLvl, 0.0);
            h5domain.close();

            // Conditional data
            H5::Group h5data(h5file.createGroup("DATA"));

            Vector<Real> invCount(count.begin(), count.end());
            for (auto& t : invCount) t = 1.0/(t + 1e-15);

            // Others
            for (int i=0; i<nAvgVars; i++) {
                auto& a = varMeanVal[i];
                auto& v = varStdVal[i];
                switch (avgVarWeights[i])
                {
                    case w_none:
                    {
                        for (int j=0; j<a.size(); j++) {
                            a[j] *= invCount[j];
                            v[j] *= invCount[j];
                        }
                        break;
                    }
                    case w_volume:
                    {
                        for (int j=0; j<a.size(); j++) {
                            a[j] /= (volMean[j] + 1e-15);
                            v[j] /= (volMean[j] + 1e-15);
                        }
                        break;
                    }
                    case w_density:
                    {
                        for (int j=0; j<a.size(); j++) {
                            a[j] /= (rhoMean[j] + 1e-15);
                            v[j] /= (rhoMean[j] + 1e-15);
                        }
                        break;
                    }
                    case w_mass:
                    {
                        for (int j=0; j<a.size(); j++) {
                            a[j] /= (massMean[j] + 1e-15);
                            v[j] /= (massMean[j] + 1e-15);
                        }
                        break;
                    }
                }
                for (int j=0; j<v.size(); j++) {
                    v[j] = std::sqrt(v[j] - a[j]*a[j]);
                }
                writeRealData(h5data, avgVarNames[i] + "_mean", nBins, varMeanVal[i], 0.0, avgVarWeights[i], 0.0);
                writeRealData(h5data, avgVarNames[i] + "_std", nBins, varStdVal[i], 0.0, avgVarWeights[i], 4);
                writeRealData(h5data, avgVarNames[i] + "_min", nBins, varMinVal[i], 0.0, w_notdefined, 4);
                writeRealData(h5data, avgVarNames[i] + "_max", nBins, varMaxVal[i], 0.0, w_notdefined, 4);
            } // Others

            // Count
            int countTot = 0;
            for (const auto& c : count) countTot += c;
            writeIntData(h5data, "count", nBins, count, countTot, w_notdefined, 4);

            // Volume
            Real volTot = 0;
            for (const auto& v : volMean) volTot += v;
            volTot *= vol0;
            for (auto& v : volMean) v *= vol0;
            for (auto& v : volStd) v *= vol0*vol0;
            writeRealData(h5data, "volume_sum", nBins, volMean, volTot, w_notdefined, 4);
            for (int i=0; i<volMean.size(); i++) {
                volMean[i] *= invCount[i];
            }
            for (int i=0; i<volStd.size(); i++) {
                volStd[i] = std::sqrt(volStd[i]*invCount[i] - volMean[i]*volMean[i]);
            }
            writeRealData(h5data, "volume_mean", nBins, volMean, 0.0, w_none, 4);
            writeRealData(h5data, "volume_std", nBins, volStd, 0.0, w_none, 4);

            // Density
            writeRealData(h5data, "density_sum", nBins, rhoMean, 0.0, w_notdefined, 4);
            for (int i=0; i<rhoMean.size(); i++) {
                rhoMean[i] *= invCount[i];
            }
            for (int i=0; i<rhoStd.size(); i++) {
                rhoStd[i] = std::sqrt(rhoStd[i]*invCount[i] - rhoMean[i]*rhoMean[i]);
            }
            writeRealData(h5data, "density_mean", nBins, rhoMean, 0.0, w_none, 4);
            writeRealData(h5data, "density_std", nBins, rhoStd, 0.0, w_none, 4);

            // Mass
            for (auto& m : massMean) m *= vol0;
            for (auto& m : massStd) m *= vol0*vol0;
            Real massTot = 0;
            for (const auto& m : massMean) massTot += m;
            writeRealData(h5data, "mass_sum", nBins, massMean, massTot, w_notdefined, 4);
            for (int i=0; i<massMean.size(); i++) {
                massMean[i] *= invCount[i];
            }
            for (int i=0; i<massStd.size(); i++) {
                massStd[i] = std::sqrt(massStd[i]*invCount[i] - massMean[i]*massMean[i]);
            }
            writeRealData(h5data, "mass_mean", nBins, massMean, 0.0, w_none, 4);
            writeRealData(h5data, "mass_std", nBins, massStd, 0.0, w_none, 4);

            h5data.close();
            h5file.close();

        } // IOprocessor
    } // iPlot - plotFileNames

    amrex::Finalize();
    return 0;
}
