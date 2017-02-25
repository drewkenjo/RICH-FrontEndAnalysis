#include<iostream>
#include<stdint.h>
#include<TChain.h>
#include<TCanvas.h>
#include<THnSparse.h>
#include<TH1D.h>
#include"feio.h"


int main(int argc, char** argv)
{
 uint8_t ipixels[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};

 TChain* tt = new TChain("tmap");
 for(int iarg=1;iarg<argc;iarg++)
	tt->AddFile(argv[iarg]);

 int channel, pixel = 22;
 for(channel=0;channel<64;channel++)
	if(ipixels[channel]==pixel)
		break;

 UInt_t key, val;
 UChar_t chan;
 tt->SetBranchAddress("chan", &chan);
 tt->SetBranchAddress("key", &key);
 tt->SetBranchAddress("val", &val);

 Int_t bins[] = {3000, 3000, 64};
 Double_t xmin[] = {.5, .5, 0.5};
 Double_t xmax[] = {3000.5, 3000.5, 64.5};
 THnSparseI hmap("hmap", "hmap", 3, bins, xmin, xmax);

 int nen = tt->GetEntries();
 for(int ien=0;ien<nen;ien++){
	tt->GetEntry(ien);
	if(chan==channel){
		UShort_t adc1 = key & 0xfff;
		UShort_t adc0 = (key>>12) & 0xfff;
		UChar_t ipix = (key>>24) & 0xff;

		double vv[] = {(double) adc0, (double) adc1, (double) ipix};

		std::cout<<channel<<" "<<adc0<<" "<<adc1<<" "<<ipix<<std::endl;
		hmap.Fill(vv, val);
	}
 }
 delete tt;

 TCanvas* c1 = new TCanvas("c1","c1", 1100, 800);
 TH1D* h1 = hmap.Projection(0);
 h1->Draw();
 c1->Print("1.pdf");
 

 return 0;
}

