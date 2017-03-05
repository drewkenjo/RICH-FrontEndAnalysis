#include<iostream>
#include<stdint.h>
#include<TCanvas.h>
#include<TChain.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TH2I.h>
#include<TF1.h>
#include<TLine.h>
#include"feio.h"

using RICHfrontend::NPIXELS;
using RICHfrontend::chan2pix;
using RICHfrontend::pix2chan;

int main(int argc, char** argv)
{
 if(argc<3 || TString(argv[1]).Atoi()<1 || TString(argv[1]).Atoi()>64){
	std::cerr<<"USAGE: "<<argv[0]<<" pix#[1-64] data file[s]"<<std::endl;
	exit(222);
 }

 int ipixel = TString(argv[1]).Atoi()-1;
 int ichannel = RICHfrontend::pix2chan[ipixel];

 UShort_t fadc[NPIXELS];
 TChain *tt = new TChain("h22");
 for(int iarg=2;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);
 tt->SetBranchAddress("fadc", fadc);

//////////////////////////////////////////////////////////////
 TH2I *hadc2adc[NPIXELS], *zadc2adc[NPIXELS], *uadc2adc[NPIXELS];
 for(int ipix=0;ipix<NPIXELS;ipix++){
	hadc2adc[ipix] = new TH2I(Form("hadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix+1), 300,0.5,3000.5, 300,0.5,3000.5);
	zadc2adc[ipix] = new TH2I(Form("zadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix+1), 650,200.5,1500.5, 650,200.5,1500.5);
 }

 long unsigned int nen = tt->GetEntries();
 for(int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);

	for(int ich=0;ich<NPIXELS;ich++){
		int ipix = chan2pix[ich]-1;

		hadc2adc[ipix]->Fill(fadc[ichannel], fadc[ich]);
		zadc2adc[ipix]->Fill(fadc[ichannel], fadc[ich]);
	}
 }


//////////////////////////////////////////////////////////////
 gStyle->SetLineScalePS(.7);
 gStyle->SetOptStat(0);
 gStyle->SetPadBottomMargin(0.0001);
 gStyle->SetPadTopMargin(0.0001);
 gStyle->SetPadLeftMargin(0.0001);
 gStyle->SetPadRightMargin(0.0001);
 gStyle->SetOptLogz();

 TString pdfname(Form("adc2adc_ondata_4pixel%02d.pdf",ipixel+1));
 TCanvas* c1 = new TCanvas("c1","c1", 2000,2000);
 c1->Divide(8,8,0.00001,0.00001);
 c1->Print(pdfname+"[");

 for(int ipix=0;ipix<NPIXELS;ipix++){
	c1->cd(ipix+1);
	hadc2adc[ipix]->Draw("col");
 }
 c1->Print(pdfname);

 for(int ipix=0;ipix<NPIXELS;ipix++){
	c1->cd(ipix+1);
	zadc2adc[ipix]->Draw("col");
 }
 c1->Print(pdfname);

 c1->cd()->Clear();
 c1->Divide(3,3,0.00001,0.00001);
 for(int icol=0;icol<3;icol++)
 for(int irow=0;irow<3;irow++){
	c1->cd(icol*3+irow+1)->SetGrid();

	int ipix = ipixel+irow-1 + (icol-1)*8;
	if(ipix<0 || ipix>63) continue;

	TH1* zx = zadc2adc[ipix]->ProjectionX();
	TH1* zy = zadc2adc[ipix]->ProjectionY();
	double x0 = zx->GetBinCenter(zx->GetMaximumBin());
	double y0 = zy->GetBinCenter(zy->GetMaximumBin());
	delete zx, zy;

	zadc2adc[ipix]->Draw("col");
 }
 c1->Print(pdfname);

 c1->Print(pdfname+"]");

 for(int ipix=0;ipix<NPIXELS;ipix++)
	delete hadc2adc[ipix], zadc2adc[ipix];
 delete tt;

 return 0;
}

