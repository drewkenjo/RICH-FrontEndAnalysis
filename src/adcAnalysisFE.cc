#include<iostream>
#include<TChain.h>
#include<TH1I.h>
#include<TH2I.h>
#include<TMultiGraph.h>
#include<TGraph.h>
#include<TCanvas.h>
#include<TString.h>
#include<TStyle.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(TString);
	~goodRICHEvent();
	void Fill(rawEvent&);
  private:
	TH1I* h1[NCHANNELS];
	TString fname;
	int ipix[64];
};

goodRICHEvent::goodRICHEvent(TString _fname):fname(_fname){
     int ipixel[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};
	std::copy(ipixel, ipixel+64, ipix);

     for(int ich=0;ich<NCHANNELS;ich++)
		h1[ich]  = new TH1I(Form("h1_%03d",ich), Form("Channel %d, pixel %d; ADC",ich,ipix[ich%64]), 4100, -0.5, 4099.5 );
}

void goodRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);
	for(int ich=0;ich<NCHANNELS;ich++)
		h1[ich]->Fill(fadc[ich]);
}

goodRICHEvent::~goodRICHEvent(){
	TCanvas* c1 = new TCanvas("c1","c1",1100,400);
	c1->SetLogy();
	c1->Print("adc_analysis.pdf[");

	TH2I* h2 = new TH2I("h2","laser intensity;iASIC [0-2] + MAPMT's iCOLUMN [1-8];MAPMT's iROW [0-8]",24,0.5,24.5,8,.5,8.5);

	for(int ich=0;ich<NCHANNELS;ich++){
	     h1[ich]->GetXaxis()->SetRange(h1[ich]->FindFirstBinAbove(1), h1[ich]->FindLastBinAbove(1));
		TH1I* hh = (TH1I*) h1[ich]->Clone("hh");
	     hh->GetXaxis()->SetRange(hh->GetMaximumBin()+20, std::max(hh->GetMaximumBin()+21, h1[ich]->FindLastBinAbove(1)));
		hh->SetFillColor(kCyan);
		h1[ich]->Draw();
		hh->Draw("same");
		c1->Print("adc_analysis.pdf");

		int iy = (ipix[ich%64]-1)/8;
		int ix = ipix[ich%64]-iy*8 + (2-ich/64)*8;
		h2->Fill(ix, 8-iy, hh->Integral());
		delete hh;
	}
	gStyle->SetOptStat(0);
	c1->SetLogy(0);
	h2->Draw("colz");
	c1->Print("adc_analysis.pdf");
	c1->Print("adc_analysis.pdf]");

	for(int ich=0;ich<NCHANNELS;ich++)
		delete h1[ich];
}



//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;

 TChain *tt = new TChain("h22");
 for(int iarg=1;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);
 tt->SetBranchAddress("trigID", &rawEv.trigID);
 tt->SetBranchAddress("timeStamp", &rawEv.timeStamp);

 tt->SetBranchAddress("fadc", rawEv.fadc);

 tt->SetBranchAddress("fnedge", &rawEv.fnedge);
 tt->SetBranchAddress("ftdc", rawEv.ftdc);
 tt->SetBranchAddress("fchan", rawEv.fchan);
 tt->SetBranchAddress("fpolar", rawEv.fpolar);

 goodRICHEvent ev(argv[1]);

 int nen = tt->GetEntries();
 for(int ien=0; ien<nen; ien++){
	tt->GetEntry(ien);
	ev.Fill(rawEv);
 }

 return 0;
}

