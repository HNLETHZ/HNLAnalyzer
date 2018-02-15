/*
 * Merge a number of .root files into one output file.
 * Inspired by: http://www-d0.fnal.gov/Run2Physics/cs/howto/tmb_analyze/TMBAnalyze-8.html
 *
 * Execute the function like this from the command line:
 *
 * echo outputfile.root input1.root input2.root ... | root -b -l merge_root.C+
 *
 * This macro assumes that the TTree is named 'T'. Change
 * the corresponding line below if this is different for your case.
 */

#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

// Change this line if your TTree has a different name
const char *TTreeName = "T";

void MergeRootFiles() {
	using namespace std;

	TChain chain(TTreeName);

	string outfile;
	cin >> outfile;

	Int_t total_events = 0;
	std::string filename;
	while(cin >> filename) {
        	TFile *f = new TFile(filename.c_str());
		
		TTree *tree = (TTree *)f->Get(TTreeName);
		Int_t nEvents = (Int_t )tree->GetEntries();
		cout << "Adding "<< nEvents << " events from file: " << filename << endl;

		chain.Add(filename.c_str());
		
		total_events += nEvents;
	}

	cout << "Output file: " << outfile << endl;

	cout << "Merging trees...patience..." << endl;
	chain.SetMaxTreeSize(10000000000); // Set the maximum size of a Tree file = 8 Gb
	chain.Merge(outfile.c_str());

	cout << "Total Events: " << total_events << endl;
}
