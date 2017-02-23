#include<iostream>
#include<list>
#include<TFile.h>
#include<TTree.h>
#include<TF1.h>
#include<TH2F.h>
#include<TH3F.h>
#include<TPad.h>
#include<TCanvas.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(int);
	~goodRICHEvent();
	void Fill(rawEvent&);

  private:
	TH3F* h_dt_t0_adc[NCHANNELS];
	int nasic;
};


goodRICHEvent::goodRICHEvent(int _nasic):nasic(_nasic){
	int ipix[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};
	for(int ich=0;ich<NCHANNELS;ich++){
		h_dt_t0_adc[ich]  = new TH3F(Form("h_dt_t0_adc_%03d",ich), Form("Channel %d, pixel %d;TDC hit, duration;TDC hit, leading time;ADC",ich,ipix[ich%64]),
			100,0.5,100.5,		50, 100.5, 150.5,		340,450.5,2150.5);
//			100,0.5,100.5,		275, 100.5, 1200.5,		340,450.5,2150.5);
	}
}


void goodRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);

	for(int ichan=0; ichan<NCHANNELS; ichan++){
		for(int iedge=0; iedge<ftdc[ichan].size(); iedge++)
		if(fpolar[ichan][iedge]==fLeadingEdge
			&& iedge<ftdc[ichan].size()-1
			&& fpolar[ichan][iedge+1]==fFallingEdge
			){
				ftime[ichan].push_back(ftdc[ichan][iedge]);
				fdur[ichan].push_back(ftdc[ichan][iedge+1] - ftdc[ichan][iedge]);
				iedge++;
		}

		if(ftime[ichan].size()>0)
			h_dt_t0_adc[ichan]->Fill(fdur[ichan][0], ftime[ichan][0], fadc[ichan]);
		else
			h_dt_t0_adc[ichan]->Fill(-1, -1, fadc[ichan]);
	}
}


goodRICHEvent::~goodRICHEvent(){
	TString pdfname("adc_vs_t0.pdf");
	TCanvas* c1 = new TCanvas("c1","c1",1000,1200);
	c1->Divide(1,3,.0001,.0001);
	c1->Print(pdfname+"[");
	for(int ich=0; ich<NCHANNELS; ich++){
		TH3F* h3 = h_dt_t0_adc[ich];

		c1->cd(1);
		gPad->SetTopMargin(0);
		gPad->SetGrid();

		TH1* hyz = h3->Project3D("yz");
		hyz->Draw("colz");

//////////////
		c1->cd(2);
		gPad->SetTopMargin(0);
		gPad->SetGrid();

		TH1* hyx = h3->Project3D("yx");
		TF1* f1 = new TF1("f1", "[0]+expo(1)",0,100);
		f1->SetParameters(100,4,-0.04);
		hyx->Fit(f1,"Q0");
		f1->SetParameter(0,f1->GetParameter(0)+7);
		hyx->Draw("colz");
		f1->Draw("same");

//////////////
		c1->cd(3);
		gPad->SetTopMargin(0);
		gPad->SetGrid();

		TH1* hadc0 = (TH1*) h3->Project3D("hadc_z0");
		hadc0->SetLineStyle(2);
		hadc0->SetLineWidth(2);
		TH1* hadc1 = (TH1*) h3->Project3D("z1 NUF NOF");
		TH1* hadc2 = (TH1*) h3->Project3D("z2 NUF NOF");
		hadc2->Reset();
		hadc2->SetLineColor(kRed);
		for(int ibx=1;ibx<=h3->GetNbinsX();ibx++)
		for(int iby=1;iby<=h3->GetNbinsY();iby++)
			if(h3->GetYaxis()->GetBinCenter(iby) > f1->Eval(h3->GetXaxis()->GetBinCenter(ibx)))
			for(int ibz=1;ibz<=h3->GetNbinsZ();ibz++)
			for(int ien=0;ien<h3->GetBinContent(ibx,iby,ibz);ien++)
				hadc2->Fill(hadc2->GetBinCenter(ibz));
		double hmax = hadc1->GetMaximum();
		hadc1->GetXaxis()->SetRange(hadc0->FindFirstBinAbove(0)-10, hadc1->FindLastBinAbove(hmax/50));
		hadc1->Draw();
		hadc0->Draw("same");
		hadc1->Draw("same");
		hadc2->Scale(hmax/hadc2->GetMaximum());
		hadc2->Draw("same hist");

//////////////

		c1->Print(pdfname);

		if(ich==63 && nasic==2) ich+=64;

		delete hyz, hyx, hadc0, hadc1, hadc2, f1, h3;
	}
	c1->Print(pdfname+"]");

	for(int ich=0; ich<NCHANNELS; ich++)
	if(h_dt_t0_adc[ich])
		delete h_dt_t0_adc[ich];
}


//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev(TString(argv[1]).Contains("2ASIC") ? 2 : 3);

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

Some other changes
