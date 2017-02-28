#include<iostream>
#include<TChain.h>
#include<TFile.h>
#include<TTree.h>
#include"feio.h"
#include<TStyle.h>
#include<TCanvas.h>
#include<TH1I.h>
#include<TF1.h>


int main(int argc, char** argv)
{
 TH1I* h1[RICHfrontend::NCHANNELS];
 for(int ich=0;ich<RICHfrontend::NCHANNELS;ich++)
	h1[ich] = new TH1I(Form("h1_%03d",ich), Form("channel %03d",ich), 400,300.5,700.5);
 TF1* f1 = new TF1("f1","gaus",0,1000);
 UShort_t m3sig[RICHfrontend::NCHANNELS]={};

 UShort_t fadc[RICHfrontend::NCHANNELS];
 TChain* h22 = new TChain("h22");
 for(int iarg=1;iarg<argc;iarg++)
	h22->AddFile(argv[iarg]);
 h22->SetBranchAddress("fadc", fadc);

 TFile* ff = new TFile(Form("skim_%s", basename(argv[1])), "recreate");
 TTree* tt = (TTree*) h22->CloneTree(0);

 unsigned long int nen = h22->GetEntries();
 for(int ien=0;ien<std::min(1000000.0, (double) nen);ien++){
	h22->GetEntry(ien);
	for(int ich=0;ich<RICHfrontend::NCHANNELS;ich++)
		h1[ich]->Fill(fadc[ich]);
 }

 for(int ich=0;ich<RICHfrontend::NCHANNELS;ich++){
	double ymax = h1[ich]->GetMaximum();
	h1[ich]->GetXaxis()->SetRange(h1[ich]->FindFirstBinAbove(0.001*ymax)-5, h1[ich]->FindLastBinAbove(0.001*ymax)+5);
	h1[ich]->Fit(f1,"Q");

	m3sig[ich] = (UShort_t) (f1->GetParameter(1) + 4*fabs(f1->GetParameter(2)));
	delete h1[ich];
 }

 for(int ien=0;ien<nen;ien++){
	h22->GetEntry(ien);

	for(int ich=0;ich<RICHfrontend::NCHANNELS;ich++)
	if(fadc[ich] > m3sig[ich]){
		tt->Fill();
		break;
	}
 }
 ff->Write();

 delete h22, tt, ff;
 return 0;
}
