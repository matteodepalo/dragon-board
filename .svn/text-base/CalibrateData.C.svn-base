#define anaChannel 5
#define trigChannel 0
#define DOMINO_NCELL 1024
#define DOMINO_NCH 8
#define DOMINO_DEPTH 12  // number of bit
#define NCALIBFILES 11
#include "datastructure.h"
#include "io.h"
#include "functions.h"

const Double_t Pi = TMath::Pi();

void CalibrateData(Int_t nevt,Int_t startEv = 1, char *PedFile = "drs4_20100311_t_ped.root",
		Int_t mVStart = -250 , Int_t mVStep = 50) {

	// create progress bar
	TGHProgressBar *gProgress = ProgressBar("Calibrazione dati");

	// Redefine DOMINO Depth in ADC counts
	const Float_t DominoDepthADC = pow(2, DOMINO_DEPTH);

	// create list of histograms for pedestals
	TList *grPedList = new TList();
	TGraphErrors *grPed;

	// create list of histograms for channels
	TList *hCellCalibList = new TList();
	TH1F *hCellCalib;

	TList *grCellCalibList = new TList();
	TGraphErrors *grCellCalib;

	int mV[NCALIBFILES];

	for (int iFile = 0; iFile < NCALIBFILES; iFile++) {
		mV[iFile] = mVStart;
		mVStart += mVStep;
	}

	char *calibrationFile;
	char *calibrationFileArray[NCALIBFILES];
	calibrationFileArray[0] = "drs4_20100311_m250.dat";
	calibrationFileArray[1] = "drs4_20100311_m200.dat";
	calibrationFileArray[2] = "drs4_20100311_m150.dat";
	calibrationFileArray[3] = "drs4_20100311_m100.dat";
	calibrationFileArray[4] = "drs4_20100311_m050.dat";
	calibrationFileArray[5] = "drs4_20100311_p000.dat";
	calibrationFileArray[6] = "drs4_20100311_p050.dat";
	calibrationFileArray[7] = "drs4_20100311_p100.dat";
	calibrationFileArray[8] = "drs4_20100311_p150.dat";
	calibrationFileArray[9] = "drs4_20100311_p200.dat";
	calibrationFileArray[10] = "drs4_20100311_p250.dat";

	for (int iFile = 0; iFile < NCALIBFILES; iFile++) {
		for (int ch = 0; ch < DOMINO_NCELL; ch++) {
			//
			TString title = "Calibration signal file:";
			title += iFile;
			title += " ch:";
			title += ch;
			hCellCalib = new TH1F(title,title, DominoDepthADC, -DominoDepthADC, DominoDepthADC);
			hCellCalibList->Add(hCellCalib);
		}
	}


	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		grCellCalib = new TGraphErrors(NCALIBFILES);
		grCellCalibList->Add(grCellCalib);
	}


	// calculate or read pedestals from file
	grPedList = OpenPedestals(PedFile);
	grPedData = (TGraphErrors *) grPedList->At(anaChannel);
	grPedTrig = (TGraphErrors *) grPedList->At(trigChannel);

	// create gauss function
	TF1 *fgauss = new TF1("fgauss", "TMath::Gaus(x,[0],[1],0)", -DOMINO_NCELL, DOMINO_NCELL);
	fgauss->SetParameter(0,0.);
	fgauss->SetParameter(1,1.);
	fgauss->SetParLimits(0, 0., DominoDepthADC);
	fgauss->SetParLimits(1, 0.1, 20.);

	TCanvas *ctest = new TCanvas("ChannelTest", "ChannelTest", 800, 600);
	ctest->Divide(3, 4);

	gProgress->Reset();
	gProgress->SetMax(nevt*NCALIBFILES);

	gSystem->ProcessEvents();


	for (int iFile = 0; iFile < NCALIBFILES; iFile++) {

		// open file

		calibrationFile = calibrationFileArray[iFile];
		FILE *fdata = OpenDataFile(calibrationFile);
		struct channel_struct *p;
		struct channel_struct *dep;

		Double_t refval=0, reftmp = 0;
		Double_t PedVal, itmp;

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
				<< calibrationFile << endl;

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
		Int_t iTrig = 0;
		Int_t flagEnd = 0;
		Double_t chtmp, chtrig;
		Double_t ratio;
		Double_t mean, rms;

		// loop on events

		while (ievt <= nevt && !flagEnd) {

			fread((void *) &event_data, 1, sizeof(event_data), fdata);
			if (feof(fdata))
				flagEnd = 1;

			p = (struct channel_struct *) &event_data.ch[0]; // read bunch of data
			dep = (struct channel_struct *) &event_data.ch[1]; // read bunch of data

			//now anaChannel analysis

			p += anaChannel;

			// create array of domino points

			if(iFile!=5) {

				Double_t DominoYval[DOMINO_NCELL];

				for(int ch = 0; ch < DOMINO_NCELL; ch++) {
					grPedData->GetPoint(ch, itmp, PedVal);
					chtmp = (Double_t)(p->data[ch]);
					chtmp = chtmp - PedVal;
					DominoYval[ch] = chtmp;
				}

				// find the domino point with max y-value

				if(iFile<5) {
					Double_t DominoMax = DominoYval[0];
					for(int ch = 1; ch < DOMINO_NCELL; ch++)
						if(DominoYval[ch] < DominoMax)
							DominoMax = DominoYval[ch];
				} else {
					Double_t DominoMax = DominoYval[0];
					for(int ch = 1; ch < DOMINO_NCELL; ch++)
						if(DominoYval[ch] > DominoMax)
							DominoMax = DominoYval[ch];
				}

				refval = DominoMax / 2;

				if(ievt == 1)
					cout << "refval="<< refval<<endl;

				vector<int> *indexes = GetIntersections(DominoYval, refval);
				vector<int> realindexes;
				realindexes.push_back(indexes->at(0));

				for(int i = 0; i < indexes->size() - 1; i++)
					if(indexes->at(i + 1) - indexes->at(i) > 20)
						realindexes.push_back(indexes->at(i + 1));

				if(ievt == 1) {
					for(i = 0; i < realindexes.size();i++) {
						cout << "Index " << i << ": " << realindexes.at(i) << endl;
					}
				}

				Int_t i;

				if(TMath::Abs(DominoYval[0]) > TMath::Abs(refval))
					i = 1;
				else
					i = 0;


				while(i + 1 < realindexes.size()) {
					for(int ch = realindexes.at(i) + 20; ch < realindexes.at(i + 1) - 20; ch++) {
						// Read pedestal value for this cell
						grPedData->GetPoint(ch, itmp, PedVal);
						chtmp = (Double_t)(p->data[ch]); // data value
						chtmp = chtmp - PedVal;
						if(TMath::Abs(chtmp) > TMath::Abs(refval))
							((TH1 *) hCellCalibList->At(iFile*DOMINO_NCELL+ch))->Fill(chtmp);
						//cout << ch << " " << iFile << " " << chtmp << endl;
					}
					i+=2;
				}
			}
			else {
				for(int ch = 0; ch < DOMINO_NCELL; ch++) {
					// Read pedestal value for this cell
					grPedData->GetPoint(ch, itmp, PedVal);
					chtmp = (Double_t)(p->data[ch]); // data value
					chtmp = chtmp - PedVal;
					((TH1 *) hCellCalibList->At(iFile*DOMINO_NCELL+ch))->Fill(chtmp);
					//cout << ch << " " << iFile << " " << chtmp << endl;
				}
			}

			gProgress->Increment(1);
			gSystem->DispatchOneEvent(kTRUE);
			ievt++; // next event
		}


		TH1 *hCellTmp;
		for(ch = 0; ch < DOMINO_NCELL;ch++) {
			hCellTmp = ((TH1 *) hCellCalibList->At(iFile*DOMINO_NCELL+ch));
			//hCellTmp->Fit("gaus", "Q");
			//mean = (hCellTmp->GetFunction("gaus"))->GetParameter(1);
			//rms = (hCellTmp->GetFunction("gaus"))->GetParameter(2);
			mean = hCellTmp->GetMean();
			rms = hCellTmp->GetRMS();
			((TGraphErrors *) (grCellCalibList->At(ch)))->SetPoint(iFile, (Double_t) mV[iFile], mean);
			((TGraphErrors *) (grCellCalibList->At(ch)))->SetPointError(iFile, 0., rms);
		}


		ctest->cd(iFile + 1);
		hCellTmp = ((TH1 *) hCellCalibList->At(567 + iFile*DOMINO_NCELL));
		hCellTmp->Fit("gaus","Q");
		hCellTmp->DrawCopy();
	}

	TString OutFile = "CalibrationDataNew";
	OutFile += nevt;
	OutFile += "events.root";
	TFile *f = new TFile(OutFile, "RECREATE");
	for (int ch = 0; ch < DOMINO_NCELL; ch++) {
		TString key = "CalibDataCell";
		key += ch;
		((TGraphErrors*) grCellCalibList->At(ch))->Write(key);
	}
	f->Close();

	hCellCalibList->Delete();

	((TGMainFrame *) gProgress->GetParent())->CloseWindow();
	fclose(fdata);
}

