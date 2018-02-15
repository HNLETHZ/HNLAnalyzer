/*
 * Creating an Analysis Tree Class with MakeClass from one input Root file.
 *
 * Execute the function like this from the command line:
 *
 * echo input.root | root -b -l MakeAnalysisTreeClass.C
 *
 * This macro assumes that TTree is named 'T' and analysis class is named 'AnalysisTreeClass'
 * Change the corresponding lines below if this is different for your case.
 */

#include <string>
#include <iostream>

#include "TList.h"
#include "TFile.h"
#include "TTree.h"

// Change two following lines if your TTree and AnalysisTreeClass have different names
const char *TTreeName = "T";
const char *ClassName = "AnalysisTreeClass";

void MakeAnalysisTreeClass() {

	using namespace std;

	string inputRootFile;
	cin >> inputRootFile;
	
	TFile *f = new TFile(inputRootFile.c_str());
	
	cout << "Input Root file: "     << inputRootFile << endl;
	cout << "TTree name: "          << TTreeName << endl;
	cout << "Analysis tree class name: " << ClassName << endl;
	
	TTree *Tree = (TTree*)f->Get(TTreeName);
	Tree->MakeClass(ClassName);

}
