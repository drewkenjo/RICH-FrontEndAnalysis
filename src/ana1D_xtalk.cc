#include<iostream>
#include<stdint.h>
#include<TCanvas.h>
#include<TChain.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TH1I.h>
#include<TH2I.h>
#include<TF1.h>
#include<TText.h>
#include<TLine.h>
#include<TBox.h>
#include"feio.h"

using RICHfrontend::NPIXELS;
using RICHfrontend::chan2pix;
using RICHfrontend::pix2chan;

int main(int argc, char** argv)
{
 if(argc<2){
	std::cerr<<"USAGE: "<<argv[0]<<" pix# data file[s]"<<std::endl;
	exit(222);
 }

 UShort_t fadc[NPIXELS];
 TChain *tt = new TChain("h22");
 for(int iarg=1;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);
 tt->SetBranchAddress("fadc", fadc);

//////////////////////////////////////////////////////////////
 TH2I *hadc2adc[NPIXELS][9] = {};
 TH1I *zadc[NPIXELS], *xadc0[NPIXELS], *xadc1[NPIXELS], *hped[NPIXELS];
 for(int ipix=0;ipix<NPIXELS;ipix++){
	hped[ipix] = new TH1I(Form("hped_%02d",ipix), Form("adc pix %02d",ipix+1), 400,200.5,600.5);
	zadc[ipix] = new TH1I(Form("zadc_%02d",ipix), Form("adc pix %02d;ADC",ipix+1), 1300,200.5,1500.5);

	xadc1[ipix] = new TH1I(Form("xadc1_%02d",ipix), Form("adc pix %02d",ipix+1), 1300,200.5,1500.5);
	xadc1[ipix]->SetFillColor(17);
	xadc1[ipix]->SetLineWidth(1);

	xadc0[ipix] = new TH1I(Form("xadc0_%02d",ipix), Form("adc pix %02d",ipix+1), 1300,200.5,1500.5);
	xadc0[ipix]->SetLineColor(kRed);

	for(int iadj=0; iadj<9; iadj++){
		int irow = iadj/3-1;
		int icol = iadj%3-1;
		int adjpix = ipix + irow*8 + icol;

		if(adjpix == ipix || adjpix<0 || adjpix>63) continue;
		if(ipix/8==0 && irow<0) continue;
		if(ipix/8==7 && irow>0) continue;
		if(ipix%8==0 && icol<0) continue;
		if(ipix%8==7 && icol>0) continue;

		hadc2adc[ipix][iadj] = new TH2I(Form("hadc2adc_%02d_%d",ipix,iadj), Form("Y (pix %02d) vs X (pix %02d);ADC;ADC",adjpix+1, ipix+1), 325,200.5,1500.5, 325,200.5,1500.5);
	}
 }

 long unsigned int nen = tt->GetEntries();
 for(long unsigned int ien=0;ien<std::min(100000.0,(double)nen);ien++){
	tt->GetEntry(ien);

	for(int ich=0;ich<NPIXELS;ich++){
		int ipix = chan2pix[ich]-1;
		hped[ipix]->Fill(fadc[ich]);
	}
 }

 TF1* f1 = new TF1("f1", "gaus", 0,1);
 double m4sig[NPIXELS], mm[NPIXELS];
 for(int ipix=0;ipix<NPIXELS;ipix++){
	double xped = hped[ipix]->GetMaximumBin();
	hped[ipix]->GetXaxis()->SetRange(xped-10, xped+10);
	f1->SetParameters(hped[ipix]->GetEntries(), hped[ipix]->GetMean(), hped[ipix]->GetRMS());
	hped[ipix]->Fit(f1,"QN");

	mm[ipix] = f1->GetParameter(1);
	m4sig[ipix] = f1->GetParameter(1)+10*fabs(f1->GetParameter(2));
 }

 for(long unsigned int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);

	for(int ich=0;ich<NPIXELS;ich++){
		int ipix = chan2pix[ich]-1;
		zadc[ipix]->Fill(fadc[ich]);

		/////////////////////////////////////////////////
		bool noAdjacentSignal = true;
		for(int iadj=0;iadj<9;iadj++){
			int irow = iadj/3-1;
			int icol = iadj%3-1;
			int adjpix = ipix + irow*8 + icol;
			if(hadc2adc[ipix][iadj] == 0) continue;

			int adjch = pix2chan[adjpix];

			if(fadc[adjch] > m4sig[adjpix] && fadc[adjch] > ((fadc[ich]-mm[ipix])+m4sig[adjpix])){
				hadc2adc[ipix][iadj]->Fill(fadc[ich], fadc[adjch]);
				xadc1[ipix]->Fill(fadc[ich]);

				noAdjacentSignal = false;
				break;
			}
		}
		if(noAdjacentSignal){
			xadc0[ipix]->Fill(fadc[ich]);
		}
	}
 }


//////////////////////////////////////////////////////////////
 gStyle->SetOptStat(0);
 gStyle->SetLineScalePS(1);
 gStyle->SetPadTopMargin(0);
 gStyle->SetPadRightMargin(0);
 gStyle->SetTitleColor(kMagenta);

 TBox bb;
 bb.SetFillColor(38);
 TLine ll;
 TText txt;
 txt.SetTextAlign(22);
 txt.SetTextColor(kRed);
 txt.SetTextSize(0.09);

 TString pdfname("adcXTalk.pdf");
 TCanvas* c1 = new TCanvas("c1","c1", 4000,2000);
 c1->Divide(2,1,0.00001,0.00001);
 c1->cd(1)->Divide(1,2,0.00001,0.00001);
 c1->cd(2)->Divide(3,3,0.00001,0.00001);
 c1->Print(pdfname+"[");

 for(int ipix=0;ipix<64;ipix++){
	c1->cd(1)->cd(1)->SetLogy();
	gPad->SetBottomMargin(0);
	gPad->SetGrid();
	zadc[ipix]->GetXaxis()->SetRangeUser(mm[ipix]-30, mm[ipix]+1400);
	TH1* hnew = (TH1*) zadc[ipix]->Clone("hnew");
	hnew->Draw();
	hnew->SetMaximum(xadc1[ipix]->GetMaximum()*.5);
	xadc1[ipix]->Draw("same");
	xadc0[ipix]->Draw("same");

	c1->cd(1)->cd(2);
	gPad->SetGrid();
	zadc[ipix]->GetXaxis()->SetRangeUser(m4sig[ipix]+100, mm[ipix]+1400);
	double ymax = zadc[ipix]->GetMaximum();
	zadc[ipix]->GetXaxis()->SetRangeUser(mm[ipix]-30, mm[ipix]+1400);
	zadc[ipix]->SetMaximum(1.5*ymax);
	zadc[ipix]->Draw();
	xadc1[ipix]->Draw("same");
	xadc0[ipix]->Draw("same");


	for(int iadj=0; iadj<9; iadj++){
		int irow = iadj/3-1;
		int icol = iadj%3-1;
		int adjpix = ipix + irow*8 + icol;

		c1->cd(2)->cd(iadj+1)->Clear();

		if(iadj==4){
			gPad->Range(-0.88888888,-8.8888888,8,0);
			for(int isq=0;isq<9;isq++){
				ll.DrawLine(isq,0,isq,-8);
				ll.DrawLine(0,-isq,8,-isq);
			}
			for(int iadj=0; iadj<9; iadj++){
				if(hadc2adc[ipix][iadj] == 0) continue;

				double x1 = (ipix%8) + (iadj%3-1);
				double y1 = -((ipix/8) + (iadj/3-1));
				bb.DrawBox(x1+.05,y1-.05,x1+.95,y1-.95);
			}
			txt.DrawText(ipix%8+0.5, -ipix/8-0.5, Form("%d",ipix+1));
		}
		if(hadc2adc[ipix][iadj] == 0) continue;

		hadc2adc[ipix][iadj]->Draw("colz");
	}
	c1->Print(pdfname);

	delete hnew;
 }


 c1->Print(pdfname+"]");

 for(int ipix=0;ipix<NPIXELS;ipix++){
	delete zadc[ipix], xadc0[ipix], xadc1[ipix], hped[ipix];
	for(int jpix=0;jpix<9;jpix++)
		delete hadc2adc[ipix][jpix];
 }
 delete tt;

 return 0;
}

