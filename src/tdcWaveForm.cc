#include<iostream>
#include<TFile.h>
#include<TTree.h>
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
	void Fill(rawEvent&, int);

  private:
	TH2I* h2[2][3][NCHANNELS];
};

goodRICHEvent::goodRICHEvent(){
	int ipix[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};

	for(int idelay=0; idelay<2; idelay++)
	for(int ipolar=0; ipolar<3; ipolar++)
	for(int ich=0;ich<NCHANNELS;ich++)
		h2[idelay][ipolar][ich]  = new TH2I(Form("h2_%d_%d_%03d",idelay,ipolar,ich),
			Form("Channel %d, pixel %d; TDC; Threshold",ich,ipix[ich%64]),
			1500, 0.5, 1500.5,
			81, 145, 955);
}

void goodRICHEvent::Fill(rawEvent &rev, int thr)
{
	RICHEvent::Fill(rev);

	for(int ichan=0; ichan<NCHANNELS; ichan++)
	for(int iedge=0; iedge<ftdc[ichan].size(); iedge++)
	if(fpolar[ichan][iedge]==fLeadingEdge
		&& iedge<ftdc[ichan].size()-1
		&& fpolar[ichan][iedge+1]==fTrailingEdge
		){
			int hittime = ftdc[ichan][iedge];
			int hitduration = ftdc[ichan][iedge+1] - ftdc[ichan][iedge];

			h2[0][1-fpolar[ichan][iedge]][ichan]->Fill(ftdc[ichan][iedge], thr);
			h2[0][1-fpolar[ichan][iedge+1]][ichan]->Fill(ftdc[ichan][iedge+1], thr);
			if(hittime>120){
				h2[1][1-fpolar[ichan][iedge]][ichan]->Fill(ftdc[ichan][iedge], thr);
				h2[1][1-fpolar[ichan][iedge+1]][ichan]->Fill(ftdc[ichan][iedge+1], thr);
			}

			iedge++;
	}
}

goodRICHEvent::~goodRICHEvent(){
	TCanvas* c1 = new TCanvas("c1","c1",800,600);
	c1->SetLogz();
	c1->Print("thr_vs_tdc.pdf[");

//	for(int ich=0;ich<NCHANNELS;ich++){
	int ich=22;{
		for(int idelay=0;idelay<2;idelay++){
			h2[idelay][2][ich] = (TH2I*) h2[idelay][0][ich]->Clone("htot");
			h2[idelay][2][ich]->Add(h2[idelay][1][ich]);
/*
			for(int zmax=0;zmax<5;zmax+=4)
			for(int ipolar=0; ipolar<3; ipolar++){
				h2[idelay][ipolar][ich]->GetXaxis()->SetRange(h2[0][2][ich]->FindFirstBinAbove(zmax)-10, h2[0][2][ich]->FindLastBinAbove(zmax)+20);
				h2[idelay][ipolar][ich]->GetYaxis()->SetRange(h2[0][2][ich]->FindFirstBinAbove(2,2)-2, h2[0][2][ich]->FindLastBinAbove(2,2)+2);
				h2[idelay][ipolar][ich]->Draw("colz");
				if(!(idelay>0))
					c1->Print("thr_vs_tdc.pdf");
			}
*/
		}
		h2[0][2][ich]->GetXaxis()->SetRange(h2[0][2][ich]->FindFirstBinAbove(2)-10, 190);
		h2[0][2][ich]->GetYaxis()->SetRange(h2[0][2][ich]->FindFirstBinAbove(2,2)-2, h2[0][2][ich]->FindLastBinAbove(2,2)+2);
		h2[0][2][ich]->Draw("colz");
		c1->Print("thr_vs_tdc.pdf");
	}
	c1->Print("thr_vs_tdc.pdf]");
	
	for(int idelay=0; idelay<2; idelay++)
	for(int ipolar=0; ipolar<3; ipolar++)
	for(int ich=0;ich<NCHANNELS;ich++)
		delete h2[idelay][ipolar][ich];
}



//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev;

for(int iarg=1;iarg<argc;iarg++){
 TFile* ff = new TFile(argv[iarg]);
 TTree* tt = (TTree*) ff->Get("h22");

 tt->SetBranchAddress("trigID", &rawEv.trigID);
 tt->SetBranchAddress("timeStamp", &rawEv.timeStamp);
 tt->SetBranchAddress("hold", &rawEv.hold);

 tt->SetBranchAddress("fadc", rawEv.fadc);

 tt->SetBranchAddress("fnedge", &rawEv.fnedge);
 tt->SetBranchAddress("ftdc", rawEv.ftdc);
 tt->SetBranchAddress("fchan", rawEv.fchan);
 tt->SetBranchAddress("fpolar", rawEv.fpolar);

 int nen = tt->GetEntries();
 for(int ien=0; ien<nen; ien++){
	tt->GetEntry(ien);
	ev.Fill(rawEv, iarg*10+190);
 }

 delete tt, ff;
}
 return 0;
}

