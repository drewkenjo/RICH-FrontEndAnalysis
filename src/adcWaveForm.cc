#include<iostream>
#include<TChain.h>
#include<TH2I.h>
#include<TMultiGraph.h>
#include<TGraph.h>
#include<TCanvas.h>
#include<TString.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent();
	~goodRICHEvent();
	void Fill(rawEvent&);

  private:
	TH2I* h2[NCHANNELS];
	int maxhold = 100;
	double maxz = 0;
};

goodRICHEvent::goodRICHEvent(){
	for(int ich=0;ich<NCHANNELS;ich++)
		h2[ich]  = new TH2I(Form("h2_%03d",ich), Form("Channel %d, pixel %d; ADC",ich,chan2pix[ich%64]), maxhold, 0.5, maxhold+.5, 4300, -200.5, 4099.5 );
}

void goodRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);

	for(int ichan=0; ichan<NCHANNELS; ichan++)
	for(int iedge=0; iedge<ftdc[ichan].size(); iedge++)
	if(fpolar[ichan][iedge]==fLeadingEdge
		&& iedge<ftdc[ichan].size()-1
		&& fpolar[ichan][iedge+1]==fTrailingEdge
		){
			ftime[ichan].push_back(ftdc[ichan][iedge]);
			fdur[ichan].push_back(ftdc[ichan][iedge+1] - ftdc[ichan][iedge]);
			iedge++;
	}

	for(int ich=0;ich<NCHANNELS;ich++)
		h2[ich]->Fill(hold, fadc[ich]);
}

goodRICHEvent::~goodRICHEvent(){
	TCanvas* c1 = new TCanvas("c1","c1",800,600);
	c1->SetLogz();
	c1->Print("hold_vs_adc.pdf[");

//	for(int ich=0;ich<NCHANNELS;ich++){
	int ich=22;{
		int ixmax, iymax, izmax;
		h2[ich]->GetMaximumBin(ixmax, iymax, izmax);
		h2[ich]->GetYaxis()->SetRange(iymax+10, 4000);
		maxz = std::max(h2[ich]->GetMaximum(), maxz);
		h2[ich]->GetYaxis()->SetRange(1, iymax-10);
		maxz = std::max(h2[ich]->GetMaximum(), maxz);

		h2[ich]->SetMaximum(maxz);
		h2[ich]->GetYaxis()->SetRange(h2[ich]->FindFirstBinAbove(2, 2)-50, h2[ich]->FindLastBinAbove(2,2)+50);
		h2[ich]->Draw("colz");
		c1->Print("hold_vs_adc.pdf");
	}
	c1->Print("hold_vs_adc.pdf]");
	
	for(int ich=0;ich<NCHANNELS;ich++)
		delete h2[ich];
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
 tt->SetBranchAddress("hold", &rawEv.hold);

 tt->SetBranchAddress("fadc", rawEv.fadc);

 tt->SetBranchAddress("fnedge", &rawEv.fnedge);
 tt->SetBranchAddress("ftdc", rawEv.ftdc);
 tt->SetBranchAddress("fchan", rawEv.fchan);
 tt->SetBranchAddress("fpolar", rawEv.fpolar);

 goodRICHEvent ev;

 int nen = tt->GetEntries();
 for(int ien=0; ien<nen; ien++){
	tt->GetEntry(ien);
	ev.Fill(rawEv);
 }

 return 0;
}

