#define anaChannel 5
#define DOMINO_NCELL 1024
#define DOMINO_NCH 8
#define DOMINO_DEPTH 12  // number of bit
#include "datastructure.h"
#include "io.h"
#include "functions.h"

const Double_t Pi = TMath::Pi();

void AnalyzeDataTimeResSin(char *DataFile = "drs4_peds_5buffers.dat", Int_t nevt,
		Int_t startEv = 1, char *PedFile, char *PeriodFile="Period.dat") {

	// create progress bar
	TGHProgressBar *gProgress = ProgressBar("Analisi residui temporali segnale SIN");

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

	Double_t histMin = -25;
	Double_t histMax = 25;
	Int_t histNbins = 50;


	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		//
		TString title = "Time residuals cell";
		title += ch;
		hCellTimeRes = new TH1F(title,title,histNbins, histMin, histMax);
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
	printf("nevtDataMax: %d\n", nevtDataMax - 1);

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
	Double_t chtmp;
	Double_t PedVal, itmp;
	Int_t cht, incCht, decCht, xflag, iMax, iMin, iRMS;
	Int_t fitusati = 0;
	Double_t xdiff, incDiff, decDiff, incDiffTemp, decDiffTemp, incXDiff, decXDiff;
	Double_t fitMax, fitMin;
	Double_t RMSmean, RMSsum;
	Double_t xtest, ytest;
	Double_t ratio;
	Double_t DominoXval[DOMINO_NCELL];
	Double_t DominoYval[DOMINO_NCELL];
	Double_t FitXval[DOMINO_NCELL];
	Double_t FitYval[DOMINO_NCELL];

	// create sawtooth function
	TF1 *fsaw = new TF1("fsaw", Sawtooth, 0., 1024, 4);
	fsaw->SetParameters(600., 255., 150., 150.);
	fsaw->SetParNames("amplitude", "period", "phase", "DC-Offset");

	//create square function
	TF1 *fsquare = new TF1("fsquare", Square, 0., 1024., 4);
	fsquare->SetParameters(300., 255., 150., 150.);
	fsquare->SetParNames("amplitude", "period", "phase", "DC-Offset");

	//create triangle function
	TF1 *ftriangle = new TF1("ftriangle", Triangle, 0., 1024., 4);
	ftriangle->SetParameters(300., 255., 150., 150.);
	ftriangle->SetParNames("amplitude", "period", "phase", "DC-Offset");

	Double_t *arr = OpenPeriodFile(PeriodFile);
	Double_t fixedPeriod = arr[0];
	Double_t fixedPeriodRMS = arr[1];
	Double_t fixedAmplitude = 597.2;

	cout << arr[0] << " " << arr[1] << endl;

	gProgress->Reset();
	gProgress->SetMax(nevt);

	gSystem->ProcessEvents();

	// loop on events

	//debug canvas
	TCanvas *cfitTest = new TCanvas("cfitTest", "fit tests", 1200, 780);
	cfitTest->Divide(1,nevt);

	while (ievt <= nevt && !flagEnd) {

		fread((void *) &event_data, 1, sizeof(event_data), fdata);
		if (feof(fdata))
			flagEnd = 1;

		p = (struct channel_struct *) &event_data.ch[0]; // read bunch of data
		dep = (struct channel_struct *) &event_data.ch[1]; // read bunch of data
		// goes to channel to analyze
		p += anaChannel;

		// read data, subtract pedestals values and save results in grAnaChDataTemp graph with
		// fixed error for each point (x = 0.5 and y = 2.1). Also generate an array with Domino
		// X and Y values

		TGraphErrors *grAnaChDataTemp = new TGraphErrors(DOMINO_NCELL);

		for (int ch = 0; ch < DOMINO_NCELL; ch++) {
			// Read pedestal value for this cell
			grPed->GetPoint(ch, itmp, PedVal);
			chtmp = (Double_t)(p->data[ch]); // data value
			chtmp = chtmp - PedVal;
			grAnaChDataTemp->SetPoint(ch, (Double_t) ch, chtmp);
			grAnaChDataTemp->SetPointError(ch, 0.5, 2.1);
			DominoXval[ch]= ch;
			DominoYval[ch]= chtmp;
		}

		// create fit functions
		TF1 *fsin = new TF1("fsin", sigSin, 0., 1024., 4);
		fsin->SetParameters(600., 255., 150., 150.);
		fsin->SetParNames("amplitude", "Period", "Phase", "DC-Offset");
		fsin->FixParameter(1, fixedPeriod);
		fsin->FixParameter(0, fixedAmplitude/2);

		grAnaChDataTemp->Fit("fsin", "Q");
		TF1 *fsinFit = grAnaChDataTemp->GetFunction("fsin");
		fsinFit->SetParNames("amplitude", "Period", "Phase", "DC-Offset");

		// debug draw
		cfitTest->cd(ievt);
		grAnaChDataTemp->SetMarkerStyle(20);
		grAnaChDataTemp->SetMarkerSize(0.3);
		grAnaChDataTemp->GetYaxis()->SetLabelSize(0.12);
		grAnaChDataTemp->GetXaxis()->SetLabelSize(0.12);
		grAnaChDataTemp->Draw("APE");

		Double_t fitAmplitude = TMath::Abs(fsinFit->GetParameter("amplitude"));
		Double_t fitPeriod = fsinFit->GetParameter("Period");
		Double_t chisquare = fsinFit->GetChisquare();
		Double_t fitOffset = fsinFit->GetParameter("DC-Offset");

		if(chisquare > 0.1e+06) {
			gProgress->Increment(1);
			gSystem->DispatchOneEvent(kTRUE);
			ievt++;
			continue;
		}

		fitusati++;

		fitMax = fsinFit->GetMaximum();
		fitMin = fsinFit->GetMinimum();

		// generate array with fit X and Y values

		for (int ch = 0; ch < DOMINO_NCELL; ch++) {
			// get Fit-value and save in array
			FitXval[ch] = DominoXval[ch];
			FitYval[ch] = fsinFit->Eval(FitXval[ch]);
		}

		//reset temp vars

		incCht = 0;
		decCht = 0;
		incDiff = 0;
		incDiffTemp = 0;
		decDiff = 0;
		decDiffTemp = 0;
		xflag = 1;

		// channel 0 and 1023 are computed out of the cicle

		/*// channel 0

		incDiff = DominoYval[0] - FitYval[0];
		incDiffTemp = DominoYval[0] - FitYval[0 + 1];

		if(incDiff < 0) {
			for (int i = 0; i < DOMINO_NCELL; i++, incDiff = incDiffTemp) {
				incDiffTemp = DominoYval[0] - FitYval[i];

				if (incDiff*incDiffTemp < 0) {
					if(abs(incDiffTemp) < abs(incDiff))
						incCht = i;
					else
						incCht = i - 1;
					break;
				}
			}
		}
		else incCht = 0;

		xdiff = FitXval[incCht] - DominoXval[0];
		((TH1F *) hCellTimeResList->At(0))->Fill(xdiff);*/

		/*// channel 1023

		decDiff = DominoYval[DOMINO_NCELL-1] - FitYval[DOMINO_NCELL-1];
		decDiffTemp = DominoYval[DOMINO_NCELL-1] - FitYval[DOMINO_NCELL - 2];

		if(decDiff < 0) {
			for (int j = DOMINO_NCELL -1 ; j >= 0 ; j--, decDiff = decDiffTemp) {
				decDiffTemp = DominoYval[DOMINO_NCELL-1] - FitYval[j];

				if (decDiff*decDiffTemp < 0) {
					if(abs(decDiffTemp) < abs(decDiff))
						decCht = j;
					else
						decCht = j + 1;
					break;
				}
			}
		}
		else decCht = DOMINO_NCELL-1;

		xdiff = FitXval[decCht] - DominoXval[DOMINO_NCELL-1];
		((TH1F *) hCellTimeResList->At(DOMINO_NCELL-1))->Fill(xdiff); */

		// now DOMINO_NCELL - 2 points

		//reset temp vars

		incCht = 0;
		decCht = 0;
		incDiff = 0;
		incDiffTemp = 0;
		decDiff = 0;
		decDiffTemp = 0;
		xflag = 1;

		for (int ch = 1; ch < DOMINO_NCELL - 1; ch++) {

			if(FitYval[ch] > fitOffset + fitAmplitude/2 || FitYval[ch] < fitOffset - fitAmplitude/2) {
				continue;
			}

			incDiff = DominoYval[ch] - FitYval[ch];
			incDiffTemp = DominoYval[ch] - FitYval[ch + 1];

			decDiff = DominoYval[ch] - FitYval[ch];
			decDiffTemp = DominoYval[ch] - FitYval[ch - 1];

			if(abs(incDiffTemp) < abs(incDiff) || (incDiff*incDiffTemp < 0 && abs(decDiffTemp) > abs(decDiff))) {
				for (int i = ch; i < DOMINO_NCELL; i++, incDiff = incDiffTemp) {
					incDiffTemp = DominoYval[ch] - FitYval[i];

					if (incDiff*incDiffTemp < 0) {
						if(abs(incDiffTemp) < abs(incDiff))
							incCht = i;
						else
							incCht = i - 1;
						break;
					}
				}
				xflag = 1;
			}
			else if(abs(decDiffTemp) < abs(decDiff) || (decDiff*decDiffTemp < 0 && abs(incDiffTemp) > abs(incDiff))) {
				for (int j = ch; j >= 0 ; j--, decDiff = decDiffTemp) {
					decDiffTemp = DominoYval[ch] - FitYval[j];

					if (decDiff*decDiffTemp < 0) {
						if(abs(decDiffTemp) < abs(decDiff))
							decCht = j;
						else
							decCht = j + 1;
						break;
					}
				}
				xflag = -1;
			}

			if(xflag == 1)
				xdiff = FitXval[incCht] - DominoXval[ch];
			else
				xdiff = FitXval[decCht] - DominoXval[ch];

			((TH1F *) hCellTimeResList->At(ch))->Fill(xdiff);

			//cout << "ch: " << ch << " xdiff: " << xdiff << endl;
		}
		gProgress->Increment(1);
		gSystem->DispatchOneEvent(kTRUE);
		ievt++; // next event

	}

	//get RMS and mean for each cell and save them in grTimeRes graph

	TGraphErrors *grTimeRes = new TGraphErrors(DOMINO_NCELL);

	TH1 *hCellTmp;

	RMSsum = 0.;

	for(ch = 0; ch < DOMINO_NCELL; ch++) {
		hCellTmp = ((TH1F *) hCellTimeResList->At(ch));
		//mean = hCellTmp->GetMean();
		Double_t hentries = hCellTmp->GetEntries();
		if(hentries != 0.) {
			hCellTmp->Fit("gaus", "Q");
			TF1 *fgausFit = hCellTmp->GetFunction("gaus");
			Double_t chisquare = fgausFit->GetChisquare();
			Double_t mean = fgausFit->GetParameter(1);
			Double_t rms = fgausFit->GetParameter(2);
			//if(chisquare < 1.0e-11) {
			RMSsum += rms;
			grTimeRes->SetPoint(ch, (Double_t) ch, mean);
			grTimeRes->SetPointError(ch, 0 ,rms);
			//    				cout << "ch: "  << ch <<  " chisquare: " << chisquare << endl;
			//}
			/*else {
				mean = hCellTmp->GetMean();
				rms = hCellTmp->GetRMS();
				RMSsum += rms;
				grTimeRes->SetPoint(ch, (Double_t) ch, mean);
				grTimeRes->SetPointError(ch, 0 ,rms);
			}*/
		}
	}
	RMSmean = RMSsum / DOMINO_NCELL;

	// save root file with graph containing time residuals.
	TString OutFile = "TimeResidualsGraph";
	OutFile += nevt;
	OutFile += "events.root";
	TFile *fg = new TFile(OutFile, "RECREATE");
	TString key = "TimeResGraph";
	grTimeRes->Write(key);
	fg->Close();

	TString OutFile = "TimeResidualsHisto";
	OutFile += nevt;
	OutFile += "events.root";
	TFile *fh = new TFile(OutFile, "RECREATE");
	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		TString key = "TimeResHCell";
		key += ch;
		((TH1F *) hCellTimeResList->At(ch))->Write(key);
	}

	fh->Close();

	//plots


	cout << "fit scartati: " << nevt - fitusati << endl;

	TCanvas *c = new TCanvas("c", "ChannelTest", 1000, 1000);
	c->Divide(3, 3);
	for (int i = 1; i < 10; i++) {
		c->cd(i);
		hCellTmp = ((TH1 *) hCellTimeResList->At(i*100 + 22));
		mean = hCellTmp->GetMean();
		rms = hCellTmp->GetRMS();
		// gaussian fit for each plotted cell
		hCellTmp->Fit("gaus", "Q");
		hCellTmp->DrawCopy();

	}

	cout << "Draw grTimeRes" << endl;
	TString Title = "grTimeRes";
	TCanvas *grTimeResCanvas = new TCanvas("grTimeRes", Title, 1200, 780);
	grTimeRes->SetMarkerStyle(20);
	grTimeRes->SetMarkerSize(0.3);
	grTimeRes->GetYaxis()->SetLabelSize(0.12);
	grTimeRes->GetXaxis()->SetLabelSize(0.12);
	grTimeRes->Draw("APE");

	hCellTimeResList->Delete();

	((TGMainFrame *) gProgress->GetParent())->CloseWindow();
	fclose(fdata);

}

