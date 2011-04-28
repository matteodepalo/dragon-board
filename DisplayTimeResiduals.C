#define anaChannel 5
#define DOMINO_NCELL 1024
#define DOMINO_NCH 8
#define DOMINO_DEPTH 12  // number of bit
#include "datastructure.h"
#include "io.h"
#include "functions.h"

void DisplayTimeResiduals(char *TimeResGraphFile="TimeResidualsGraph1000events.root",
		char *TimeResHistoFile="TimeResidualsHisto1000events.root") {

	TGraphErrors *grTimeRes = OpenTimeResFileGraph(TimeResGraphFile);

	/*TList *hCellTimeResList = OpenTimeResFileHList(TimeResHistoFile);

	TH1 *hCellTmp;

	TCanvas *c = new TCanvas("c", "ChannelTest", 1000, 1000);

	hCellTmp = ((TH1 *) hCellTimeResList->At(50));

	Double_t mean = hCellTmp->GetMean();
	Double_t rms = hCellTmp->GetRMS();
	//gaussian fit for each plotted cell
	//hCellTmp->Fit("gaus", "Q");
	hCellTmp->DrawCopy();*/


	cout << "Draw grTimeRes" << endl;
	TString Title = "grTimeRes";
	TCanvas *grTimeResCanvas = new TCanvas("grTimeRes", Title, 1200, 780);
	grTimeRes->SetMarkerStyle(20);
	grTimeRes->SetMarkerSize(0.3);
	grTimeRes->GetYaxis()->SetLabelSize(0.10);
	grTimeRes->GetXaxis()->SetLabelSize(0.10);
	grTimeRes->Draw("APE");
}
