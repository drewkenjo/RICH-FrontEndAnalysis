#include<iostream>
#include<list>
#include<TFile.h>
#include<TTree.h>
#include<TH2F.h>
#include<TPad.h>
#include<TCanvas.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(){};
	~goodRICHEvent(){};
	void Fill(rawEvent&);
	double getADC(int ich){return fadc[ich];};
	double hasTDC(int ich){return ftime[ich].size()>0 && ftime[ich][0]<150;};
};

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
}


//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 TH1F* h0 = new TH1F("h0", "", 1,0,1);
 for(int iarg=1;iarg<argc;iarg++){
	if(!TString(argv[iarg]).Contains(TRegexp(".root$"))) continue;
	h0->Fill(TString(TString(argv[iarg])(TRegexp("thr[0-9]*"))), 1);
	h0->Fill("all", 1);
 }
 h0->LabelsDeflate();
 h0->LabelsOption("a");

 int nasic = TString(argv[1]).Contains("2ASIC") ? 2 : 3;
 TH2F* h2[NCHANNELS];
 for(int ich=0;ich<NCHANNELS;ich++){
		int nconf = h0->GetNbinsX();
		h2[ich]  = new TH2F(Form("h2_%03d",ich), Form("Channel %d, pixel %d; ADC",ich,chan2pix[ich%64]), 4100,0.5,4100.5, nconf, .5, nconf+.5);
		for(int iconf=1;iconf<=nconf;iconf++)
			h2[ich]->GetYaxis()->SetBinLabel(iconf, h0->GetXaxis()->GetBinLabel(iconf));
 }

 rawEvent rawEv;
 goodRICHEvent ev;

 double v0 = h0->GetXaxis()->FindBin("all");
 double w0 = 1.0/h0->GetBinContent(h0->GetXaxis()->FindBin("all"));
 for(int iarg=1;iarg<argc;iarg++){
	if(!TString(argv[iarg]).Contains(TRegexp(".root$"))) continue;
	TString bname(TString(argv[iarg])(TRegexp("thr[0-9]*")));
	double v1 = h0->GetXaxis()->FindBin(bname);
	double w1 = 1.0/h0->GetBinContent(h0->GetXaxis()->FindBin(bname));

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
		for(int ich=0;ich<NCHANNELS;ich++){
			if(ev.hasTDC(ich))
				h2[ich]->Fill(ev.getADC(ich), v1, w1);

			h2[ich]->Fill(ev.getADC(ich), v0, w0);
		}
	}
	delete tt, ff;
 }

 TCanvas* c1 = new TCanvas("c1","c1",800,400);
 c1->SetLogz();
 c1->Print("adc_tdc_corr.pdf[");
 for(int ich=0;ich<NCHANNELS;ich++){
	h2[ich]->GetXaxis()->SetRange(h2[ich]->FindFirstBinAbove(10), h2[ich]->FindLastBinAbove(1));
	h2[ich]->Draw("colz");
	c1->Print("adc_tdc_corr.pdf");

	if(ich==63 && nasic==2) ich+=64;
 }
 c1->Print("adc_tdc_corr.pdf]");

 return 0;
}

