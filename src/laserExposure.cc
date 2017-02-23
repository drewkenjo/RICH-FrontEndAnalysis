#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<TH2I.h>
#include<TGraph.h>
#include<TCanvas.h>
#include<TStyle.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(int);
	~goodRICHEvent();
	void Fill(rawEvent&);
	void Print();

  private:
	TH1I* h1[NCHANNELS];
	TH2I* htot;
	TCanvas* c1;
	int nasic;
};

goodRICHEvent::goodRICHEvent(int _nasic):nasic(_nasic){
	c1 = new TCanvas("c1", "c1", nasic*400, 400);
	c1->Print("laser_exposure.pdf[");
     for(int ich=0;ich<NCHANNELS;ich++)
		h1[ich]  = new TH1I(Form("h1_%03d",ich), Form("Channel %d; ADC",ich), 4100, -0.5, 4099.5 );

	htot = new TH2I("htot", "laser exposure", nasic*8, 0.5,8*nasic+.5, 8, 0.5, 8.5);
}

void goodRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);
	for(int ich=0;ich<NCHANNELS;ich++)
		h1[ich]->Fill(fadc[ich]);
}

void goodRICHEvent::Print()
{
     int ipixel[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};
	TH2I* hperfile = new TH2I("hperfile", "laser exposure", nasic*8, 0.5,8*nasic+.5, 8, 0.5, 8.5);
	for(int ichan=0;ichan<NCHANNELS;ichan++){
		int ipmt = ichan/64;

		if(nasic==2 && ipmt==1) continue;
		if(nasic==2 && ipmt==2) ipmt--;

		int ipx = ipixel[ichan%64]-1;
		int icol = ipx%8+1 + (nasic-ipmt-1)*8;
		int irow = ipx/8+1;

		double ww = h1[ichan]->Integral(h1[ichan]->GetMaximumBin()+50, h1[ichan]->GetNbinsX());
		hperfile->Fill(icol, irow, ww);
		h1[ichan]->Reset();
	}
	htot->Add(hperfile);
	hperfile->Draw("colz");
	c1->Print("laser_exposure.pdf");

	delete hperfile;
}

goodRICHEvent::~goodRICHEvent(){
	htot->Draw("colz");
	c1->Print("laser_exposure.pdf");
	c1->Print("laser_exposure.pdf]");

	delete htot, c1;
}



//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev(TString(argv[1]).Contains("2ASIC") ? 2:3);

 for(int iarg=1;iarg<argc;iarg++){
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
	ev.Print();

	delete tt, ff;
 }

 return 0;
}

