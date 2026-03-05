#include <iostream>
#include <filesystem>
#include <math.h>
#include "pzcmesh.h"
#include "TPZElementMatrixT.h"

//#include <opencv2/opencv.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
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
    EbcInletId = 2,
    EbcOutletId = 3,
    EbcNoFlux = 4,
    ELagrange = 5
};

//TODO
std::set<int> vugIds;
std::set<int> vugBcIds;

//TPZCompMesh* HdivMesh(TPZGeoMesh *);
//TPZCompMesh* Pressuremesh(TPZGeoMesh *, int order);

TPZCompMesh *CreateFluxMesh(TPZGeoMesh *, std::set<int> &volId, std::set<int> &bcId, int &orderp);

TPZCompMesh *CreatePressureMesh(TPZGeoMesh *,std::set<int> &volId, std::set<int> &bcId,int order);

void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId);

void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs,int typeMesh);

void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);

void MeshWithSegmentVugs(TPZGeoMesh *gmesh);

void PrintCompMesh(TPZCompMesh *cmesh);

void CreateInterfaceGeoEls(TPZGeoMesh *gmesh);

void InsertInterfaceEls(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh);

TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);

void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);

void PrintCompMesh(TPZCompMesh *cmesh);

void SideOrientation(TPZCompMesh *cmesh);

void insertAtomicMaterialsf(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs);

void insertAtomicMaterialsp(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs);
