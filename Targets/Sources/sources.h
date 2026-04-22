#include <iostream>
#include <fstream>
#include <filesystem>
#include <json.hpp>
#include <math.h>

using json = nlohmann::json;

#include "pzcmesh.h"
#include "TPZElementMatrixT.h"
#include "pzmanvector.h"
#include "TPZGeoMeshTools.h"
#include "TPZCompMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzvec.h"
#include "TPZVTKGenerator.h"
#include <fstream>
#include "TPZMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZNullMaterialCS.h"
#include "TPZNullMaterial.h"
#include "TPZAnalysis.h"
//#include "TPZCreateMultiphysicsSpace.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
//#include "TPZStepSolver.h"
//using std::cout;
//using std::endl;
//using std::cin;
#include <set>
#include "TPZAnalyticSolution.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHDivApproxCreator.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzintel.h"



enum MatID{
    EMatId = 1,
    EVugId = 600,
    EVugBcId = 100,
    EFracId = 700,
    EFracBcId = 200,
    EbcInletId = 2,
    EbcOutletId = 3,
    EbcTop = 4,
    EbcBottom = 5,
    ELagrange = 6
};

struct DomData {
    std::string name = "none"; // name of the gmesh tag
    int matId = 0;
    double permeability = 0.0;
    int dim = 2;
};

struct BcData {
    std::string name = "none"; // name of the bc
    int matId = 0;
    int type = 0;
    TPZManVector<double, 3>  value = {0.0, 0.0, 0.0};
};

class ReadJson
{
public:

    //Constructor
    //@fileName path to the input .json file
    //@This constructor reads and parses the json file,
    //@storing all simulation parameters internally
    ReadJson(std::string fileName);
    //Methods
    //

    //@Returns the name of the mesh (used for identification/output)
    std::string MeshName();
    //@Returns the directory or full path to the mesh file
    std::string MeshFile();
    //@Returns a const reference to the vector containing domain data
    //@Each DomData entry represents a material region in the domain
    const std::vector<DomData> &DomainData();
    //@Returns the approximation type used in the simulation
    //@(e.g., mixed, continuous, etc.)
    int approxType();
    //@Returns the problem type identifier
    //@(used to select physics/solver strategy)
    int problemType();
    //@Returns the polynomial order for pressure approximation
    int pressOrder();
    //@Returns the spatial dimension of the problem
    int dim();
    //@Returns the mesh resolution level (refinement indicator)

    int resolution();
    //@Returns a vector with boundary condition data
    //@Each BcData represents one boundary condition entry
    //
    
    std::vector<BcData> BCInput();
    //@Returns boundary condition data specifically for fractures
    std::vector<BcData> FracBCInput();
    //
    //    MATERIAL SETS

    //
    //@Set of material IDs associated with vug boundary conditions
    //
    //
    std::set<int> fvugBcIds;
    //@Set of material IDs associated with vug regions

    std::set<int> fvugIds;
    //@Set of material IDs associated with fracture boundary conditions

    std::set<int> ffracBcIds;
    //@Set of material IDs associated with fracture regions

    std::set<int> ffracIds;


private:
    //@Raw json object storing the parsed input file
        json fInputFile;
        
        //@Mesh name (used for identification)
        std::string fMeshName;
        
        //@Mesh directory or file path
        std::string fMeshDirectory;

        //@Vector storing domain/material data
        std::vector<DomData> fDomainDataVec;
        
        //@Approximation type identifier
        int fApproxType;

        //@Problem type identifier
        int fProblemType;
        
        //@Polynomial order for pressure space
        int fPresspOrder;

        //@Problem dimension (2D or 3D)
        int fDim;
        
        //@Mesh resolution level
        int fResolution;
        
        //@Permeability for standard domain
        double fPerm;

        //@Permeability for vug regions
        double fVugPerm;

        //@Permeability for fracture regions
        double fFracPerm;

        //@Fluid viscosity
        double fVisc;

        //@Boundary condition data (standard)
        std::vector<BcData> fBcDataVec;

        //@Boundary condition data for fractures
        std::vector<BcData> fFracBcVec;
    
};

class VugsApproxSpaceGenerator
{
public:
    //Constructor
    VugsApproxSpaceGenerator();
    //Create TPZGeomesh using the gmsh mesh and material id/dimension information
    //@inputData
    //@filename
    //@meshName
    TPZGeoMesh* generateGMeshWithPhysTagVec(ReadJson inputData, std::string filename, std::string meshName);
    
    //Create Hdiv mesh giving the specifics mat ids and TPZ mesh
    //@A TPZGeoMesh mesh is needed as a base.
    //@volId contains a set with  material ids of volume type.
    //"orderp"??? why not just "order" as in "CreatePressureMesh"
    //@bcId contains a set with  material ids of volume type -1.
    //@input data have all data from .json input file.
    TPZCompMesh *CreateFluxMesh(TPZGeoMesh *, std::set<int> &volId, std::set<int> &bcId, int &orderp, ReadJson inputData);
    //Create pressure mesh giving the specifics mat ids and TPZ mesh
    //@A TPZGeoMesh mesh is needed as a base.
    //@volId contains a set with  material ids of volume type.
    //order??? why not just orderp as in "CreatePressureMesh"
    //@bcId contains a set with  material ids of volume type -1.
    //@input data have all data from .json input file.
    TPZCompMesh *CreatePressureMesh(TPZGeoMesh *, std::set<int> &volId, std::set<int> &bcId,int order, ReadJson inputData);
    //Create computational mesh of multiphysicis type
    //@Meshvec contains both the Pressure and Hdiv mesh
    //@input data have all data from .json input file.
    TPZMultiphysicsCompMesh *CreateMultiMesh(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> meshvec, ReadJson inputData);
    //Create computational mesh based on a geometric mesh
    //@gmesh is the TPZGeoMesh
    //@input data have all data from .json input file.
    TPZCompMesh *CreateMesh(TPZGeoMesh* gmesh, ReadJson inputData);
    //Extracts material ids information from a geometric mesh
    //@volId contains a set with  material ids of volume type.
    //@bcId contains a set with  material ids of volume type -1.
    void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId);
    //Insert some objects in a computational mesh of type 0 "pressure", 1 "flux" based on a list of matids.
    //@matIdsEls contains a set with  material ids of volume type
    //@matIdsBcs contains a set with  material ids of volume type -1.
    //@typeMesh for flux or pr meshe
    //@cmesh is the computational mesh
    //@input data have all data from .json input file.
    void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsEls, std::set<int> matIdsBcs, int typeMesh, ReadJson inputData);
    //Set the same connect index for a group of elements with the same mat id.
    //@cmesh is the computational mesh
    //@gmesh is the geometric mesh
    void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);
    //Identify vug and fracture entities and give a specific mat id for each one
    //Vugs contourns are in the range of and volume elements in the range of
    //Fracture contourns(points) are in the range of  XXXX and volume entities in the range of
    //@gmesh is the TPZGeoMesh
    //@input data have all data from .json input file.
    void MeshWithSegmentPhil(ReadJson inputData, TPZGeoMesh *gmesh);
    //Identify vug and fracture entities and give a specific mat id for each one
    //Vugs contourns are in the range of and volume elements in the range of
    //Fracture contourns(points) are in the range of  XXXX and volume entities in the range of
    //@gmesh is the TPZGeoMesh
    //@input data have all data from .json input file.
    void MeshWithSegment(ReadJson inputData, TPZGeoMesh *gmesh);
    //Assign value of 1 for the normal side orientation
    //@input data have all data from .json input file.
    void SideOrientation(TPZCompMesh *cmesh, ReadJson inputData);
    //Assing value of -1 times  the original gel orientation for fracture normal side orientation
    //@cmesh is the computational mesh.
    void SideOrientation1D(TPZCompMesh *cmesh);
    //*
    //@gmesh is the geometric mesh
    //@cmesh is the computational mesh.
    void DuplicateConnectFracture(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);
    //Set condensed "true" for points of fracture limits.
    //@cmesh is the computational mesh.
    void CondenseEndFrac(TPZCompMesh* cmesh);
    //Create interface between Hdiv(boundary) elements and pressure elements for vugs/fractures
    //@gmesh is the geometric mesh
    void CreateInterfaceGeoEls(TPZGeoMesh *gmesh);
    //*
    //@cmesh is the computational mesh.
    //@gmesh is the geometric mesh
    //input Data is not used?!
    void InsertInterfaceEls(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh, ReadJson inputData);
    //Matrix assemble and solve is done here.
    //@an is the analysis of the problem
    //@cmesh is the computational mesh.
    //input Data is not used?!
    void Solve(TPZLinearAnalysis* an, TPZCompMesh* cmesh, ReadJson inputData);
    //This function is not used.
    void PostProcess(ReadJson inputData);
    //Print computational mesh in both .txt and .vtk format,
    void PrintCompMesh(TPZCompMesh *cmesh);
    //Print geometric mesh in both .txt and .vtk format.
    void PrintGeoMesh(TPZGeoMesh *gmesh);
};


