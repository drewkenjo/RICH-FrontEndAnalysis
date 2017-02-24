#include<iostream>
#include<list>
#include<TFile.h>
#include<TTree.h>
#include<THnSparse.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent();
	~goodRICHEvent();
	void Fill(rawEvent&);

  private:
	THnSparseI *hh[NCHANNELS];
	int ipix[64];
};


goodRICHEvent::goodRICHEvent(){
	int ipixels[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};
	std::copy(ipixels, ipixels+64, ipix);

	int nbins[3] = {3000, 3000, 64};
	double vmin[3] = {0.5, 0.5, 0.5};
	double vmax[3] = {3000.5, 3000.5, 64.5};

	for(int ich=0;ich<NCHANNELS;ich++)
		hh[ich] = new THnSparseI(Form("hh%03d",ich), Form("channel %03d",ich), 3, nbins, vmin, vmax);
}


void goodRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);

	for(int ichan=0; ichan<NCHANNELS; ichan++){
		int iasic = ichan/64;
		for(int ich=iasic*64;ich<(iasic+1)*64;ich++){
			double vals[] = {(double) fadc[ichan], (double) fadc[ich], (double) ipix[ich%64]};
			hh[ichan]->Fill(vals);
		}
	}
}


goodRICHEvent::~goodRICHEvent(){
}


//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev;

 for(int iarg=1;iarg<argc;iarg++){
	if(!TString(argv[iarg]).Contains(TRegexp(".root$"))) continue;
	TString bname(TString(argv[iarg])(TRegexp("thr[0-9]*")));

	TFile* ff = new TFile(argv[iarg]);
	TTree* tt = (TTree*) ff->Get("h22");

	tt->SetBranchAddress("trigID", &rawEv.trigID);
	tt->SetBranchAddress("timeStamp", &rawEv.timeStamp);

	tt->SetBranchAddress("fadc", rawEv.fadc);

	tt->SetBranchAddress("fnedge", &rawEv.fnedge);
	tt->SetBranchAddress("ftdc", rawEv.ftdc);
	tt->SetBranchAddress("fchan", rawEv.fchan);
	tt->SetBranchAddress("fpolar", rawEv.fpolar);

	int nen = tt->GetEntries();
	for(int ien=0; ien<nen; ien++){
		tt->GetEntry(ien);
		ev.Fill(rawEv);
	}

	delete tt, ff;
 }

 return 0;
}

