// Merge only evt_tree from an explicit list of ROOT files into one output ROOT file.
//
// Usage:
//   root -l -b -q 'macros/merge_evt_tree.C("prepped_data","prepped_data/merged_evt_tree.root")'
//
// Optional arguments:
//   recursive = false  -> search subdirectories when true
//   forceOverwrite = false -> refuse to overwrite output file unless true

#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TSystem.h>
#include <TSystemDirectory.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

bool ends_with(const std::string &value, const std::string &suffix) {
    if (value.size() < suffix.size()) return false;
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string join_path(const std::string &a, const std::string &b) {
    if (a.empty()) return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

void collect_root_files(const std::string &dirPath, bool recursive, std::vector<std::string> &files) {
    TSystemDirectory dir(dirPath.c_str(), dirPath.c_str());
    TList *entries = dir.GetListOfFiles();
    if (!entries) {
        return;
    }

    std::vector<std::string> local;
    TIter next(entries);
    while (TObject *obj = next()) {
        const char *nameC = obj->GetName();
        if (!nameC) continue;
        std::string name(nameC);
        if (name == "." || name == "..") continue;

        std::string fullPath = join_path(dirPath, name);
        if (gSystem->AccessPathName(fullPath.c_str(), kReadPermission) == true) {
            continue;
        }

        FileStat_t st;
        if (gSystem->GetPathInfo(fullPath.c_str(), st) != 0) {
            continue;
        }

        if (R_ISDIR(st.fMode)) {
            if (recursive) {
                collect_root_files(fullPath, recursive, files);
            }
            continue;
        }

        if (ends_with(name, ".root")) {
            local.push_back(fullPath);
        }
    }

    std::sort(local.begin(), local.end());
    files.insert(files.end(), local.begin(), local.end());
}

void collect_root_files_from_manifest(const std::string &manifestPath, std::vector<std::string> &files) {
    std::ifstream in(manifestPath);
    if (!in.is_open()) {
        std::cerr << "ERROR: Could not open file list: " << manifestPath << std::endl;
        return;
    }

    std::string path;
    while (std::getline(in, path)) {
        if (!path.empty()) {
            files.push_back(path);
        }
    }
}

bool file_has_evt_tree(const std::string &filePath) {
    std::unique_ptr<TFile> file(TFile::Open(filePath.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Could not open ROOT file: " << filePath << std::endl;
        return false;
    }
    TObject *obj = file->Get("evt_tree");
    if (!obj) {
        std::cerr << "ERROR: Missing evt_tree in file: " << filePath << std::endl;
        return false;
    }
    if (!obj->InheritsFrom("TTree")) {
        std::cerr << "ERROR: evt_tree exists but is not a TTree in file: " << filePath << std::endl;
        return false;
    }
    return true;
}

}  // namespace

int merge_evt_tree(const char *inputSource,
                   const char *outputFile,
                   bool recursive = false,
                   bool forceOverwrite = false) {
    if (!inputSource || !outputFile) {
        std::cerr << "ERROR: inputSource and outputFile are required." << std::endl;
        return 1;
    }

    const std::string inDir(inputSource);
    const std::string outFile(outputFile);

    if (inDir.empty() || outFile.empty()) {
        std::cerr << "ERROR: inputSource and outputFile cannot be empty." << std::endl;
        return 1;
    }

    if (gSystem->AccessPathName(inDir.c_str(), kReadPermission) != false) {
        std::cerr << "ERROR: Input path is not readable or does not exist: " << inDir << std::endl;
        return 1;
    }

    FileStat_t inStat;
    if (!forceOverwrite && gSystem->AccessPathName(outFile.c_str(), kFileExists) == false) {
        std::cerr << "ERROR: Output file already exists. Refusing to overwrite without force flag: " << outFile << std::endl;
        return 1;
    }

    std::vector<std::string> files;
    if (gSystem->GetPathInfo(inDir.c_str(), inStat) != 0) {
        std::cerr << "ERROR: Could not inspect input path: " << inDir << std::endl;
        return 1;
    }

    if (R_ISDIR(inStat.fMode)) {
        collect_root_files(inDir, recursive, files);
    } else {
        collect_root_files_from_manifest(inDir, files);
    }

    if (files.empty()) {
        std::cerr << "ERROR: No input ROOT files found: " << inDir << std::endl;
        return 1;
    }

    // Prevent accidentally feeding the output file back into the chain.
    files.erase(std::remove(files.begin(), files.end(), outFile), files.end());

    if (files.empty()) {
        std::cerr << "ERROR: No valid input files remain after excluding output path." << std::endl;
        return 1;
    }

    for (const auto &path : files) {
        if (!file_has_evt_tree(path)) {
            std::cerr << "ERROR: Validation failed. Aborting merge." << std::endl;
            return 1;
        }
    }

    TChain chain("evt_tree");
    for (const auto &path : files) {
        chain.Add(path.c_str());
    }

    if (chain.GetEntries() < 0) {
        std::cerr << "ERROR: Failed to read entries from input evt_tree chain." << std::endl;
        return 1;
    }

    TFile out(outFile.c_str(), "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "ERROR: Could not create output ROOT file: " << outFile << std::endl;
        return 1;
    }

    out.cd();
    std::cout << "Merging " << files.size() << " files into: " << outFile << std::endl;
    TTree *merged = chain.CloneTree(-1, "fast");
    if (!merged) {
        std::cerr << "ERROR: CloneTree returned null. Merge failed." << std::endl;
        return 1;
    }

    merged->SetName("evt_tree");
    merged->Write("evt_tree", TObject::kOverwrite);
    int num_entries = merged->GetEntries();
    out.Close();

    std::cout << "Merged " << files.size() << " files into: " << outFile
              << " with " << num_entries << " evt_tree entries." << std::endl;

    return 0;
}
