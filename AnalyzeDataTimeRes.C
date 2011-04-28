#define anaChannel 5
#define DOMINO_NCELL 1024
#define DOMINO_NCH 8
#define DOMINO_DEPTH 12  // number of bit
#include "datastructure.h"
#include "io.h"
#include "functions.h"

const Double_t Pi = TMath::Pi();

void AnalyzeDataTimeRes(char *DataFile = "drs4_peds_5buffers.dat", Int_t nevt,
		Int_t startEv = 1, char *PedFile) {


	// Redefine DOMINO Depth in ADC counts
	const Float_t DominoDepthADC = pow(2, DOMINO_DEPTH);

	// open file

	FILE *fdata = OpenDataFile(DataFile);
	struct channel_struct *p;
	struct channel_struct *dep;

	// create list of histograms for pedestals
	TList *grPedList = new TList();
	TGraphErrors *grPed;

	// create list of histograms for channels
	TList *hCellTimeResList = new TList();
	TH1F *hCellTimeRes;

	//create square function
	TF1 *fsquare = new TF1("fsquare", Square, 0., 1024., 4);
	fsquare->SetParameters(600., 100., 10., 10.);
	fsquare->SetParNames("amplitude", "period", "phase", "DC-Offset");


	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		//
		TString title = "Time residuals cell";
		title += ch;
		hCellTimeRes = new TH1F(title,title, 2*((Int_t) DOMINO_NCELL), (Double_t) -DOMINO_NCELL, (Double_t) DOMINO_NCELL);
		hCellTimeResList->Add(hCellTimeRes);
	}


	// calculate or read pedestals from file
	grPedList = OpenPedestals(PedFile);
	grPed = (TGraphErrors *) grPedList->At(anaChannel);


	// Count number of events in data file
	int nevtDataMax = 0;
	while (!feof(fdata)) {
		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		nevtDataMax++;
	}
	printf("nevtDataMax: %d\n", nevtDataMax);

	if (nevt > (nevtDataMax - startEv) || nevt == 0)
		nevt = nevtDataMax - startEv;
	cout << endl << "==>> Processing " << nevt << " events from file "
			<< DataFile << endl;

	rewind(fdata);


	Int_t ievt = 1;
	// go to first event (startEv)
	while (ievt < startEv) {
		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		if (feof(fdata))
			break;
		ievt++;
	}

	ievt = 1;
	Int_t flagEnd = 0;
	Int_t iMax, iMin;
	Double_t chtmp, difftmp;
	Double_t PedVal, itmp;
	Double_t fitMax, fitMin;
	Double_t mean, rms, RMSmean, RMSsum;
	Double_t xtest, ytest;
	Double_t ratio;
	Double_t ProgrMax;
	Double_t DominoXval[DOMINO_NCELL];
	Double_t DominoYval[DOMINO_NCELL];
	Double_t ProgrDiff[DOMINO_NCELL];
	Double_t ADCerr = 1.; // error on ADC counts read by DOMINO

	// loop on events

	while (ievt <= nevt && !flagEnd) {
		if ((10*(ievt/nevt))%10){
			ratio = (Double_t) ievt / nevt;
			cout << (Int_t) (ratio*100) <<"% " << "(" << ievt << "/" << nevt << ")" << endl;
		}
		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		if (feof(fdata))
			flagEnd = 1;

		p = (struct channel_struct *) &event_data.ch[0]; // read bunch of data
		dep = (struct channel_struct *) &event_data.ch[1]; // read bunch of data
		// goes to channel to analyze
		p += anaChannel;

		// read data and subtract pedestals values

		for (int ch = 0; ch < DOMINO_NCELL; ch++) {
			// Read pedestal value for this cell
			grPed->GetPoint(ch, itmp, PedVal);
			chtmp = (Double_t)(p->data[ch]); // data value
			chtmp = chtmp - PedVal;
			//to add: calibration factor
			DominoXval[ch]= ch;
			DominoYval[ch]= chtmp;
		}

		TGraphErrors *grData = new TGraphErrors(DOMINO_NCELL,DominoXval,DominoYval);
		TGraphErrors *grProgrDiff = new TGraphErrors(DOMINO_NCELL);

		for(ch = 1; ch < DOMINO_NCELL - 1; ch++) {
			difftmp = TMath::Abs(DominoYval[ch + 1] - DominoYval[ch - 1]);
			ProgrDiff[ch] = difftmp;
			grProgrDiff->SetPoint(ch,(Double_t) ch,difftmp);
			grProgrDiff->SetPointError(ch,0.,ADCerr*sqrt(2.));
		}

		Int_t periodo = 127;

		iFirstMax = 0;
		for (ch = 0; ch < 0.95*(periodo/2); ch++) {
			if(ProgrDiff[ch] > ProgrDiff[iFirstMax]) {
				iFirstMax = ch;
			}
		}
		cout << "iFirstMax: " << iFirstMax << endl;

		if(ievt == 1) {
			cout << "Draw Progressive Difference" << endl;
			TString Title = "Progressive Difference";
			TCanvas *cdiff = new TCanvas("cdiff", Title, 1200, 780);
			grProgrDiff->SetMarkerStyle(20);
			grProgrDiff->SetMarkerSize(0.3);
			grProgrDiff->GetYaxis()->SetLabelSize(0.12);
			grProgrDiff->GetXaxis()->SetLabelSize(0.12);
			grProgrDiff->Draw("APEL");
		}
		Int_t iseq=0;
		const Int_t Ngroup = 3; // number of cells grouped together
		Int_t Tcell[Ngroup]; // t1, t2, t3, ... tn
		Double_t deltaT[Ngroup];
		Axis_t XminFit = (Axis_t) (iFirstMax - periodo/4);
		Axis_t XmaxFit = (Axis_t) (iFirstMax + periodo/4);
		while (XmaxFit < DOMINO_NCELL){
			grProgrDiff->Fit("gaus","Q","",XminFit,XmaxFit);
			TF1 *fitfun = grProgrDiff->GetFunction("gaus");
			Double_t mean = fitfun->GetParameter(1); // retrieve the value for the parameter 1
			Tcell[iseq] = (int)(mean+0.5);
			iseq++;
			XminFit += (Axis_t)periodo;
			XmaxFit += (Axis_t)periodo;
			if(ievt == 1) fitfun->DrawCopy("same"); //just for check purposes
			if (iseq == Ngroup) {
				// compute the time residuals for the 3 cells
				deltaT[Tcell]=1.;
				Tcell[0] = Tcell[Ngroup-1];
				iseq = 1;
			}
		}
		cout << "mean: " << mean << endl;
		ievt++;
	}



	/*cout << "Draw Data" << endl;
	TString Title = "Data";
	TCanvas *cdata = new TCanvas("cdata", Title, 1200, 780);
	grData->SetMarkerStyle(20);
	grData->SetMarkerSize(0.3);
	grData->GetYaxis()->SetLabelSize(0.12);
	grData->GetXaxis()->SetLabelSize(0.12);
	grData->Draw("APEL");*/

	hCellTimeResList->Delete();

	fclose(fdata);
}

