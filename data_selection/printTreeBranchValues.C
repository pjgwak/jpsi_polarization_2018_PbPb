#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>

void printTreeBranchValues() {
    const char *fileName = "/mnt/Oniatree/Oniatree_miniAOD_MC_Prompt.root"; // ROOT 파일 경로
    const char *treeName = "hionia/myTree";                                 // TTree 이름
    int eventNumber = 2876801;                                              // 원하는 이벤트 번호

    // 파일 열기
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsOpen() == false) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return;
    }

    // TTree 가져오기
    TTree *tree = (TTree*)file->Get(treeName);
    if (!tree) {
        std::cerr << "Error: Cannot find tree " << treeName << " in the file." << std::endl;
        file->Close();
        return;
    }

    // 특정 이벤트 가져오기
    tree->GetEntry(eventNumber);

    // 브랜치 정보 출력
    TObjArray *branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i)
    {
        TBranch *branch = (TBranch *)branches->At(i);
        std::cout << "Branch name: " << branch->GetName() << std::endl;

        // 브랜치에 저장된 리프(Leaf) 가져오기
        TLeaf *leaf = branch->GetLeaf(branch->GetName());
        if (leaf)
        {
            int nLeaves = leaf->GetNdata();
            for (int j = 0; j < nLeaves; ++j)
            {
                std::cout << "  Value " << j << ": " << leaf->GetValue(j) << std::endl;
            }
        }
    }

    file->Close();
}