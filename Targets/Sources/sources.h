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

struct BcData {
    std::string name = "none"; // name of the bc
    int matId = 0;
    int type = 0; // bc type 
    TPZManVector<double, 3>  value = {0.0, 0.0, 0.0}; // bc value
};

extern std::set<int> vugBcIds;
extern std::set<int> vugIds;

extern std::set<int> fracBcIds;
extern std::set<int> fracIds;

class ReadJson
{
public:

    //Constructor
    ReadJson(std::string fileName);
    
    //Methods
    std::string MeshName();

    std::string MeshFile();

    std::map<std::string, int> DomainData();

    std::map<std::string, int> VugData();

    std::map<std::string, int> FracData();

    int approxType();

    int problemType();

    int pressOrder();

    int dim();

    int resolution();

    double perm();

    double permVug();

    double permFrac();

    double visc();

    std::vector<BcData> BCInput();


private:

    json fInputFile;
    
    std::string fMeshName;
    
    std::string fMeshDirectory;

    std::map<std::string, int> fDomainData;

    std::map<std::string, int> fVugData;

    std::map<std::string, int> fFracData;
    
    int fApproxType;

    int fProblemType;
    
    int fPresspOrder;

    int fDim;
    
    int fResolution;
    
    double fPerm;

    double fVugPerm;

    double fFracPerm;

    double fVisc;

    std::vector<BcData> fBcDataVec;

};

TPZGeoMesh* generateGMeshWithPhysTagVec(ReadJson inputData, std::string filename, std::string meshName);

TPZCompMesh *CreateFluxMesh(TPZGeoMesh *, std::set<int> &volId, std::set<int> &bcId, int &orderp, ReadJson inputData);

TPZCompMesh *CreatePressureMesh(TPZGeoMesh *, std::set<int> &volId, std::set<int> &bcId,int order, ReadJson inputData);

TPZMultiphysicsCompMesh *CreateMultiMesh(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> meshvec, ReadJson inputData);

TPZCompMesh *CreateMesh(TPZGeoMesh* gmesh, ReadJson inputData);

void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId);

void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsEls, std::set<int> matIdsBcs, int typeMesh, ReadJson inputData);

void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);

void MeshWithSegmentPhil(ReadJson inputData, TPZGeoMesh *gmesh);

void MeshWithSegment(ReadJson inputData, TPZGeoMesh *gmesh);

void SideOrientation(TPZCompMesh *cmesh);

void DuplicateConnectFracture(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);

void CreateInterfaceGeoEls(TPZGeoMesh *gmesh);

//void InsertInterfaceEls(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh, std::set<int> &vugIds,std::set<int> &vugBcIds);
void InsertInterfaceEls(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh);

void Solve(TPZLinearAnalysis* an, TPZCompMesh* cmesh, ReadJson inputData);

void PostProcess(ReadJson inputData);

void PrintCompMesh(TPZCompMesh *cmesh);

void PrintGeoMesh(TPZGeoMesh *gmesh);