#include<iostream>
#include<stdint.h>
#include<TCanvas.h>
#include<TChain.h>
#include<TTree.h>
#include<TString.h>
#include<THnSparse.h>
#include<TH1D.h>
#include"feio.h"

int main(int argc, char** argv)
{
 UInt_t key, val;
 UChar_t chan;

 TChain* tt = new TChain("tmap");
 for(int iarg=1;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);

 tt->SetBranchAddress("key", &key);
 tt->SetBranchAddress("val", &val);

 Int_t vbin[] = {3000, 3000, 64};
 Double_t vmin[] = {.5, .5, 0.5};
 Double_t vmax[] = {3000.5, 3000.5, 64.5};
 THnSparseI* hmap = new THnSparseI("hmap", tt->GetTitle(), 3, vbin, vmin, vmax);

 long unsigned int nen = tt->GetEntries();
 for(int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);
	UShort_t adc1 = key & 0xfff;
	UShort_t adc0 = (key>>12) & 0xfff;
	UChar_t ipix = (key>>24) & 0xff;

	double vv[] = {(double) adc0, (double) adc1, (double) ipix};

	hmap->Fill(vv, val);
 }
 delete tt;

 TCanvas* c1 = new TCanvas("c1","c1", 1100,800);
 c1->SetLogy();
 hmap->Projection(0)->Draw();
 c1->Print("1.pdf");

 delete hmap;

 return 0;
}

