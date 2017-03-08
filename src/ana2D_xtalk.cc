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

using RICHfrontend::NCHANNELS;
using RICHfrontend::NPIXELS;
using RICHfrontend::chan2pix;
using RICHfrontend::pix2chan;

int main(int argc, char** argv)
{
 if(argc<3 || TString(argv[1]).Atoi()<1 || TString(argv[1]).Atoi()>264){
	std::cerr<<"USAGE: "<<argv[0]<<" pix#[1-264] data file[s]"<<std::endl;
	exit(222);
 }

 int ipixel = TString(argv[1]).Atoi()-1;
 int iasic = ipixel/100;
 ipixel = ipixel%100;
 int ichannel = RICHfrontend::pix2chan[ipixel];

 UShort_t _fadc[NCHANNELS];
 UShort_t *fadc = &_fadc[64*iasic];

 TChain *tt = new TChain("h22");
 for(int iarg=2;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);
 tt->SetBranchAddress("fadc", _fadc);
 long unsigned int nen = tt->GetEntries();

//////////////////////////////////////////////////////////////
 TH1I* hspe[NPIXELS];
 for(int ipix=0;ipix<NPIXELS;ipix++)
	hspe[ipix] = new TH1I(Form("hspe%02d", ipix), "", 3000,0.5,3000.5);

 for(long unsigned int ien=0;ien<std::min(100000.,(double)nen);ien++){
	tt->GetEntry(ien);
	for(int ich=0;ich<NPIXELS;ich++){
		int ipix = chan2pix[ich]-1;
		hspe[ipix]->Fill(fadc[ich]);
	}
 }

 int dfmax=0, dfzoom=0;
 int fmin[NPIXELS];
 for(int ipix=0;ipix<NPIXELS;ipix++){
	double mm = hspe[ipix]->GetMean();
	double ss = hspe[ipix]->GetRMS();

	hspe[ipix]->GetXaxis()->SetRangeUser(mm+6*ss, hspe[ipix]->GetBinCenter(hspe[ipix]->FindLastBinAbove(0))+1.0);
	double amp = (hspe[ipix]->GetMean()-mm);

	fmin[ipix] = mm-3*ss;
	
	int nbins = (2*amp + 6*ss)/16+1;

	if(nbins*16>dfmax)
		dfmax = nbins*16;

	if(ipix==ipixel)
		dfzoom = ((int) ((0.3*amp + 6*ss)/16+1))*16;

	delete hspe[ipix];
 }

//////////////////////////////////////////////////////////////
 TH2I *hadc2adc[NPIXELS], *zadc2adc[NPIXELS];
 for(int ipix=0;ipix<NPIXELS;ipix++){
	hadc2adc[ipix] = new TH2I(Form("hadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix+1),
			dfmax/4,0.5+fmin[ipixel]-dfmax/16,0.5+fmin[ipixel]+dfmax*15/16,
			dfmax/4,0.5+fmin[ipix]-dfmax/16,0.5+fmin[ipix]+dfmax*15/16);
	zadc2adc[ipix] = new TH2I(Form("zadc2adc_%02d",ipix), Form("adc vs adc for pix %02d",ipix+1),
//			dfzoom,	0.5+fmin[ipixel]-dfzoom/16,	0.5+fmin[ipixel]+dfzoom*15/16,
			dfzoom,	0.5+fmin[ipixel],	0.5+fmin[ipixel]+dfzoom,
			dfmax/2,	0.5+fmin[ipix]-dfmax/16,		0.5+fmin[ipix]+dfmax*15/16);
 }

 for(long unsigned int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);
	if(ien%(nen/100)==0 && ien>0)
		std::cerr<<"\r"<<ien*100/nen<<"% processed....";

	for(int ich=0;ich<NPIXELS;ich++){
		int ipix = chan2pix[ich]-1;

		hadc2adc[ipix]->Fill(fadc[ichannel], fadc[ich]);
		zadc2adc[ipix]->Fill(fadc[ichannel], fadc[ich]);
	}
 }
 std::cerr<<std::endl;


//////////////////////////////////////////////////////////////
 gStyle->SetLineScalePS(.7);
 gStyle->SetOptStat(0);
 gStyle->SetPadBottomMargin(0.0);
 gStyle->SetPadLeftMargin(0.05);
 gStyle->SetPadTopMargin(0.0);
 gStyle->SetPadRightMargin(0.0);
 gStyle->SetOptLogz();

 TString pdfname(Form("adc2adc_ondata_4pixel%02d.pdf",ipixel+1));
 TCanvas* c1 = new TCanvas("c1","c1", 2000,2000);
 c1->Divide(8,8,0.00001,0.00001);
 c1->Print(pdfname+"[");

 for(int ipix=0;ipix<NPIXELS;ipix++){
	c1->cd(ipix+1);
	if(ipix/8==7) gPad->SetBottomMargin(.05);
	hadc2adc[ipix]->Draw("col");
 }
 c1->Print(pdfname);

 for(int ipix=0;ipix<NPIXELS;ipix++){
	c1->cd(ipix+1);
	if(ipix/8==7) gPad->SetBottomMargin(.05);
	zadc2adc[ipix]->Draw("col");
 }
 c1->Print(pdfname);

 c1->cd()->Clear();
 c1->Divide(3,3,0.00001,0.00001);
 for(int irow=0;irow<3;irow++)
 for(int icol=0;icol<3;icol++){
	c1->cd(irow*3+icol+1)->SetGrid();
	if(irow==2) gPad->SetBottomMargin(.05);

	int ipix = ipixel+icol-1 + (irow-1)*8;
	if(ipix<0 || ipix>63) continue;

	zadc2adc[ipix]->Draw("col");
 }
 c1->Print(pdfname);

 c1->Print(pdfname+"]");

 for(int ipix=0;ipix<NPIXELS;ipix++)
	delete hadc2adc[ipix], zadc2adc[ipix];
 delete tt;

 return 0;
}

