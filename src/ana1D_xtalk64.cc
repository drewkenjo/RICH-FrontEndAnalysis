#include<iostream>
#include<stdint.h>
#include<TCanvas.h>
#include<TChain.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TLegend.h>
#include<TH1I.h>
#include<TH2I.h>
#include<TF1.h>
#include"feio.h"


using RICHfrontend::NCHANNELS;
using RICHfrontend::NPIXELS;
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
 if(TString(argv[2]).Contains("skim") || TString(argv[2]).Contains("richtree"))
	iasic = 0;

 UShort_t _fadc[NCHANNELS];
 UShort_t *fadc = &_fadc[NPIXELS*iasic];
 TChain *tt = new TChain("h22");
 for(int iarg=2;iarg<argc;iarg++)
 if(TString(argv[iarg]).Contains(TRegexp(".root$")))
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

 TF1* f1 = new TF1("f1", "gaus", 0,1);
 int dfmax=0;
 double pedm4[NPIXELS], pedm0[NPIXELS];

 for(int ipix=0;ipix<NPIXELS;ipix++){
	double mm = hspe[ipix]->GetBinCenter(hspe[ipix]->GetMaximumBin());
	hspe[ipix]->GetXaxis()->SetRange(1,hspe[ipix]->GetMaximumBin());
	double ss = hspe[ipix]->GetRMS();

	hspe[ipix]->GetXaxis()->SetRangeUser(mm-9*ss, mm+9*ss);
	f1->SetRange(mm-5*ss, mm+5*ss);
	f1->SetParameters(hspe[ipix]->GetEntries(), mm, ss);
	hspe[ipix]->Fit(f1,"QNR");
	mm = f1->GetParameter(1);
	ss = fabs(f1->GetParameter(2));

	pedm0[ipix] = mm;
	pedm4[ipix] = mm+3*ss;

	hspe[ipix]->GetXaxis()->SetRangeUser(mm+6*ss, 3000);
	double amp = (hspe[ipix]->GetMean()-mm);

	dfmax = std::max((double)dfmax, 2.5*amp+6*ss);

	delete hspe[ipix];
 }


//////////////////////////////////////////////////////////////
 TH1I *zadc[NPIXELS], *xadc9[NPIXELS], *xadc[NPIXELS][5];
 for(int ipix=0;ipix<NPIXELS;ipix++){
	int xoffset = 0.1*dfmax;

	zadc[ipix] = new TH1I(Form("zadc_%02d",ipix), Form("adc pix %02d;ADC",ipix+1), dfmax, 0.5+pedm0[ipix]-xoffset, 0.5+pedm0[ipix]+dfmax-xoffset);

	xadc9[ipix] = new TH1I(Form("xadc9_%02d",ipix), Form("adc pix %02d",ipix+1), dfmax, 0.5+pedm0[ipix]-xoffset, 0.5+pedm0[ipix]+dfmax-xoffset);
	xadc9[ipix]->SetLineColor(kRed);

	for(int isig=0;isig<5;isig++){
		xadc[ipix][isig] = new TH1I(Form("xadc%02d%d",ipix,isig), Form("adc pix %02d, %d signals in PMT",ipix+1, isig), dfmax, 0.5+pedm0[ipix]-xoffset, 0.5+pedm0[ipix]+dfmax-xoffset);
		xadc[ipix][isig]->SetLineColor(1+isig);
		xadc[ipix][isig]->SetLineStyle(3);
		xadc[ipix][isig]->SetLineWidth(2);
	}
 }


//////////////////////////////////////////////////////////////
 for(long unsigned int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);

	if(ien%(nen/100)==0 && ien>0)
		std::cerr<<"\r"<<ien*100/nen<<"% processed....";

	for(int ich=0;ich<NPIXELS;ich++){
		int ipix = chan2pix[ich]-1;
		zadc[ipix]->Fill(fadc[ich]);

		int icol = ipix%8;
		int irow = ipix/8;

		/////////////////////////////////////////////////
		bool noAdjacentSignal = true;
		for(int iadj=0;iadj<9;iadj++){
			int jrow = irow+iadj/3-1;
			int jcol = icol+iadj%3-1;
			if(jrow<0 || jcol<0 || jrow>7 || jcol>7) continue;

			int jpix = jrow*8 + jcol;
			int jch = pix2chan[jpix];

			if(jpix==ipix) continue;

			if(fadc[jch] > pedm4[jpix]){
//			if((fadc[jch] > pedm4[jpix] && fadc[ich] < pedm4[ipix]) ||
//				(fadc[ich] > pedm4[ipix] && fadc[jch] > ((fadc[ich]-pedm0[ipix])+pedm4[jpix]))){
					noAdjacentSignal = false;
					break;
			}
		}

		int signalsInMAPMT=0;
		if(noAdjacentSignal){
			xadc9[ipix]->Fill(fadc[ich]);

			for(int jpix=0;jpix<NPIXELS;jpix++){
				if(jpix==ipix) continue;
				int jch = pix2chan[jpix];

/*
				if((fadc[jch] > pedm4[jpix] && fadc[ich] < pedm4[ipix]) ||
					(fadc[ich] > pedm4[ipix] && fadc[jch] > ((fadc[ich]-pedm0[ipix])+pedm4[jpix]))){
						signalsInMAPMT++;
						break;
				}
*/

				if(fadc[jch] > pedm0[jpix] + 5*(pedm4[jpix]-pedm0[jpix]))
						signalsInMAPMT++;
			}

			if(signalsInMAPMT<5)
				xadc[ipix][signalsInMAPMT]->Fill(fadc[ich]);
		}
	}
 }
 std::cerr<<std::endl;


//////////////////////////////////////////////////////////////
 gStyle->SetLineScalePS(1.2);
 gStyle->SetOptStat(0);
 gStyle->SetPadLeftMargin(0.05);
 gStyle->SetPadRightMargin(0);

 TLegend* leg = new TLegend(0.65,0.7,1.0,0.9);
 leg->AddEntry(zadc[0],"NO conditions","l");
 leg->AddEntry(xadc9[0],"NO signal in Adjacent Pixels","l");
 for(int isig=0;isig<5;isig++)
	leg->AddEntry(xadc[0][isig],Form("%d signals in MAPMT",isig),"l");

 TString pdfname("adcXTalk64.pdf");
 TCanvas* c1 = new TCanvas("c1","c1", 1100,800);
 c1->Divide(1,2,0.00001,0.00001);
 c1->Print(pdfname+"[");

 for(int ipix=0;ipix<64;ipix++){
	c1->cd(1)->SetLogy();
	gPad->SetBottomMargin(0);
	gPad->SetGrid();

	gPad->DrawFrame(pedm0[ipix]-0.1*dfmax, 0.15, pedm0[ipix]+.9*dfmax, zadc[ipix]->GetMaximum())->SetTitle(Form("pixel %d",ipix+1));
	zadc[ipix]->Draw("same");
	xadc9[ipix]->Draw("same");
	for(int isig=0;isig<5;isig++)
		xadc[ipix][isig]->Draw("same");

	leg->Draw();

	/////////////////////////////
	c1->cd(2);
	gPad->SetTopMargin(0);
	gPad->SetGrid();

	zadc[ipix]->GetXaxis()->SetRangeUser(pedm4[ipix]+0.2*dfmax, pedm0[ipix]+dfmax*.85);
	double ymax = zadc[ipix]->GetMaximum();
	zadc[ipix]->GetXaxis()->SetRange(0,-1);

	gPad->DrawFrame(pedm0[ipix]-0.1*dfmax, 0.0, pedm0[ipix]+dfmax*.9, 1.3*ymax)->SetTitle(";ADC");
	zadc[ipix]->Draw("same");
	xadc9[ipix]->Draw("same");
	for(int isig=0;isig<5;isig++)
		xadc[ipix][isig]->Draw("same");

	c1->Print(pdfname);
 }

 c1->Print(pdfname+"]");

 for(int ipix=0;ipix<NPIXELS;ipix++){
	delete zadc[ipix], xadc9[ipix];
	for(int isig=0;isig<5;isig++)
		delete xadc[ipix][isig];
 }
 delete tt;

 return 0;
}

