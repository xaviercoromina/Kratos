//Author: Miguel Angel Celigueta. maceli@cimne.upc.edu

#if !defined(KRATOS_RENUMBERING_NODES_UTILITY_FOR_PLASMA_DYNAMICS)
#define KRATOS_RENUMBERING_NODES_UTILITY_FOR_PLASMA_DYNAMICS

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/model_part.h"

#include "pybind11/stl.h"


namespace Kratos
{
class RenumberingNodesUtilityForPlasmaDynamics
{

public:

KRATOS_CLASS_POINTER_DEFINITION(RenumberingNodesUtilityForPlasmaDynamics);

RenumberingNodesUtilityForPlasmaDynamics(ModelPart& mp1) {
    mListOfModelParts.push_back(&mp1);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}

RenumberingNodesUtilityForPlasmaDynamics(ModelPart& mp1, ModelPart& mp2): RenumberingNodesUtilityForPlasmaDynamics(mp1){
    mListOfModelParts.push_back(&mp2);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}

RenumberingNodesUtilityForPlasmaDynamics(ModelPart& mp1, ModelPart& mp2, ModelPart& mp3): RenumberingNodesUtilityForPlasmaDynamics(mp1, mp2) {
    mListOfModelParts.push_back(&mp3);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}

RenumberingNodesUtilityForPlasmaDynamics(ModelPart& mp1, ModelPart& mp2, ModelPart& mp3, ModelPart& mp4): RenumberingNodesUtilityForPlasmaDynamics(mp1, mp2, mp3) {
    mListOfModelParts.push_back(&mp4);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}

RenumberingNodesUtilityForPlasmaDynamics(ModelPart& mp1, ModelPart& mp2, ModelPart& mp3, ModelPart& mp4, ModelPart& mp5): RenumberingNodesUtilityForPlasmaDynamics(mp1, mp2, mp3, mp4) {
    mListOfModelParts.push_back(&mp5);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}


virtual ~RenumberingNodesUtilityForPlasmaDynamics(){}

void Renumber() {
    int id = 1;
    for (int i=0; i<(int)mListOfModelParts.size(); i++){
        ModelPart& mp = *mListOfModelParts[i];
        std::map<int,int>& new_to_old = mListOfMapsOfIdsNewToOld[i];

        for (int j = 0; j < (int)mp.Nodes().size(); j++){
            auto it = mp.NodesBegin() + j;
            new_to_old[id] = it->Id();
            it->SetId(id);
            id++;
        }
    }
}

void UndoRenumber() {
    for (int i=0; i<(int)mListOfModelParts.size(); i++){
        ModelPart& mp = *mListOfModelParts[i];
        std::map<int,int>& new_to_old = mListOfMapsOfIdsNewToOld[i];

        for (int j = 0; j < (int)mp.Nodes().size(); j++){
            auto it = mp.NodesBegin() + j;
            it->SetId(new_to_old[it->Id()]);
        }
    }
}



private:

std::vector<ModelPart*> mListOfModelParts;
std::vector<std::map<int,int> > mListOfMapsOfIdsNewToOld;

}; // Class RenumberingNodesUtilityForPlasmaDynamics

} // namespace Kratos.

#endif // KRATOS_RENUMBERING_NODES_UTILITY_FOR_PLASMA_DYNAMICS  defined

