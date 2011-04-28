#include <vector>
#include <algorithm>

//============== Signal Function Sin =============//
Double_t sigSin(Double_t* x, Double_t* par) {
	Double_t A = par[0], periodo = par[1], fase = par[2], offset = par[3];
	Double_t arg = x[0] + fase;
	return offset + A * sin(2 * Pi * arg / periodo);
}

//============== Function Gauss =============//
Double_t Gauss(Double_t* x, Double_t* par) {
	Double_t A = par[0], mean = par[1], rms = par[2];
	Double_t arg = x[0] - mean;
	return A * exp(-0.5 * pow(arg / rms, 2)) / (rms* sqrt(2 * Pi));
}

//============== Function Square =============//
Double_t Square(Double_t* x, Double_t* par) {
	Double_t A = par[0], periodo = par[1], fase = par[2], offset = par[3];
	Double_t arg = x[0] + fase;
	return offset + A*(pow(-1, floor(2*arg/periodo)));
}

//============== Function Triangle =============//
Double_t Triangle(Double_t* x, Double_t* par) {
	Double_t A = par[0], periodo = par[1], fase = par[2], offset = par[3];
	Double_t arg = x[0] + fase;
	//return (2/A) * (x[0] - A*floor(x[0]/A-0.5)) * pow(-1,floor(x[0]/A-0.5));
	return offset + A*(2/Pi)*asin(sin(2*Pi*arg/periodo));
	//return A*(1 - 4 * TMath::Abs(0.5-(0.5*x[0] + 0.25 - floor(0.5*x[0] + 0.25))));
}

//============== Function Sawtooth =============//
Double_t Sawtooth(Double_t* x, Double_t* par) {
	Double_t A = par[0], periodo = par[1], fase = par[2], offset = par[3];
	Double_t arg = x[0] + fase;
	return offset + A*2*(arg/periodo - floor(arg/periodo+0.5));
}

//============== Create Progress Bar =============//

TGHProgressBar *ProgressBar(char *title) {
	TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(),500,50);
	frmMain->SetWindowName(title);
	frmMain->SetCleanup(kDeepCleanup);
	TGHProgressBar *gProgress = new TGHProgressBar(frmMain, TGHProgressBar::kFancy, 500);
	gProgress->ShowPosition(kTRUE, kTRUE, "%.2f");
	gProgress->SetBarColor("green");
	frmMain->MapSubwindows();
	frmMain->MapWindow();
	return gProgress;
}

//Get indexes of intersection points between an array of domino points and a constant value //
vector<int> *GetIntersections(Double_t *points, Double_t refval) {
	Double_t lastpoint = points[0];
	int lastpush;
	Double_t diff = TMath::Abs(points[0]) - TMath::Abs(refval), lastdiff = diff;
	vector<int> *IdxArray = new vector<int>;

	for(int i = 10; i < DOMINO_NCELL; i++) {
		diff = TMath::Abs(points[i]) - TMath::Abs(refval);

		if(diff * lastdiff < 0)
			if(TMath::Abs(diff) < TMath::Abs(lastdiff))
				IdxArray->push_back(i);
			else
				IdxArray->push_back(i - 1);

		lastpoint = points[i];
		lastdiff = diff;
	}

	return IdxArray;
}
