#include<iostream>
#include<libgen.h>
#include<stdint.h>
#include<TCanvas.h>
#include<TFile.h>
#include<TTree.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TLegend.h>
#include<TBox.h>
#include<TH1I.h>
#include<TH2I.h>
#include<TF1.h>
#include"feio.h"


using RICHfrontend::NPIXELS;
using RICHfrontend::NCHANNELS;
using RICHfrontend::chan2pix;
using RICHfrontend::pix2chan;


int main(int argc, char** argv)
{
 if(argc<3 || TString(argv[1]).Atoi()<1 || TString(argv[1]).Atoi()>3){
	std::cerr<<"USAGE: "<<argv[0]<<" asic#[1-3] data file[s]"<<std::endl;
	std::cerr<<"SPE spectra for different conditions: no adjacent signal, no signal in PMT, no conditions"<<std::endl;
	exit(222);
 }

 int iasic = TString(argv[1]).Atoi()-1;


//////////////////////////////////////////////////////////////
 std::map<TString, TH1I*> hadc[NPIXELS];
 UShort_t _fadc[NCHANNELS];


//////////////////////////////////////////////////////////////
 for(int iarg=2;iarg<argc;iarg++){
	if(!TString(argv[iarg]).Contains(TRegexp(".root$")))
		continue;

	UShort_t *fadc = &_fadc[NPIXELS*iasic];
	if(TString(argv[iarg]).Contains("skim") || TString(argv[iarg]).Contains("richtree"))
		fadc = _fadc;

	std::cerr<<"\r"<<iarg<<" / "<<argc<<" processed....";

	TFile* ff = new TFile(argv[iarg]);
	TTree* tt = (TTree*) ff->Get("h22");
	tt->SetBranchAddress("fadc", _fadc);

	TString dname(dirname(argv[iarg]));
	TH1I* tadc[NPIXELS];
	for(int ipix=0;ipix<NPIXELS;ipix++){
		if(hadc[ipix].find(dname) == hadc[ipix].end())
			hadc[ipix][dname] = new TH1I(Form("hadc_%d%02d",ipix/64,ipix%64), Form("adc pix %02d;ADC",ipix%64+1), 3000,.5,3000.5);
		tadc[ipix] = hadc[ipix][dname];
	}

	long unsigned int nen = tt->GetEntries();
	for(long unsigned int ien=0;ien<nen;ien++){
		tt->GetEntry(ien);

		for(int ich=0;ich<NPIXELS;ich++){
			int ipix = chan2pix[ich]-1;
			tadc[ipix]->Fill(fadc[ich]);
		}
	}

	delete tt, ff;
 }
 std::cerr<<std::endl;


//////////////////////////////////////////////////////////////
 TF1* f1 = new TF1("f1", "gaus", 0,1);
 int dfmax=0;
 double pedm4[NPIXELS], pedm0[NPIXELS];
 for(int ipix=0;ipix<NPIXELS;ipix++){
	TH1I* hspe = hadc[ipix].begin()->second;

	double mm = hspe->GetBinCenter(hspe->GetMaximumBin());
	hspe->GetXaxis()->SetRange(1,hspe->GetMaximumBin());
	double ss = hspe->GetRMS();

	hspe->GetXaxis()->SetRangeUser(mm-9*ss, mm+9*ss);
	f1->SetRange(mm-5*ss, mm+5*ss);
	f1->SetParameters(hspe->GetEntries(), mm, ss);
	hspe->Fit(f1,"QRN");
	mm = f1->GetParameter(1);
	ss = fabs(f1->GetParameter(2));

	pedm0[ipix] = mm;
	pedm4[ipix] = mm+4*ss;

	hspe->GetXaxis()->SetRangeUser(mm+6*ss, 3000);
	double amp = (hspe->GetMean()-mm);
	hspe->GetXaxis()->SetRange(0,-1);

	dfmax = std::max((double)dfmax, 2.5*amp+6*ss);
 }
 delete f1;


//////////////////////////////////////////////////////////////
 gStyle->SetLineScalePS(1.2);
 gStyle->SetOptStat(0);
 gStyle->SetPadLeftMargin(0.05);
 gStyle->SetPadRightMargin(0);

 TLegend* leg = new TLegend(0.8,0.7,1.0,0.9);
 for(std::map<TString, TH1I*>::iterator it=hadc[0].begin(); it!=hadc[0].end(); it++)
	leg->AddEntry(it->second, it->first);

 TBox bb;
 bb.SetFillStyle(0);
 bb.SetLineWidth(3);
 bb.SetLineColor(kBlue);

 TString pdfname("adcComparison.pdf");
 TCanvas* c1 = new TCanvas("c1","c1", 1100,800);
 c1->Divide(1,2,0.00001,0.00001);
 c1->Print(pdfname+"[");

 for(int ipix=0;ipix<NPIXELS;ipix++){
	c1->cd(1)->SetLogy();
	gPad->SetGrid();

	gPad->DrawFrame(pedm0[ipix]-0.1*dfmax, 0.15, pedm0[ipix]+.9*dfmax, hadc[ipix].begin()->second->GetMaximum()*.5)->SetTitle(Form("pixel %d",ipix+1));
	int ifile = 1;
	for(std::map<TString, TH1I*>::iterator it=hadc[ipix].begin(); it!=hadc[ipix].end(); it++){
		it->second->SetLineColor(ifile++);
		it->second->Draw("same");
	}

	leg->Draw();

	/////////////////////////////
	TH1I* zadc = hadc[ipix].begin()->second;
	int adcthr = pedm4[ipix]+0.1*dfmax;
	adcthr = pedm4[ipix];

	zadc->GetXaxis()->SetRangeUser(adcthr, pedm0[ipix]+dfmax*0.8);
	double ymax = zadc->GetMaximum();
	zadc->GetXaxis()->SetRange(0,-1);

	bb.DrawBox(pedm0[ipix]-0.03*dfmax, 0.15, pedm0[ipix]+dfmax*0.06, 1.3*ymax);

	/////////////////////////////
	c1->cd(2)->SetLogy();
	gPad->SetTopMargin(0);
	gPad->SetGrid();

	gPad->DrawFrame(pedm0[ipix]-0.03*dfmax, 0.15, pedm0[ipix]+dfmax*0.06, 1.3*ymax)->SetTitle(";ADC");
	ifile = 1;
	for(std::map<TString, TH1I*>::iterator it=hadc[ipix].begin(); it!=hadc[ipix].end(); it++){
		it->second->SetLineColor(ifile++);
		it->second->Draw("same");
	}

	c1->Print(pdfname);
 }

 c1->Print(pdfname+"]");

 for(int ipix=0;ipix<NPIXELS;ipix++){
//	for(std::map<TString, TH1I*>::iterator it=hadc[ipix].begin(); it!=hadc[ipix].end(); it++)
//		delete it->second;
 }

 return 0;
}

