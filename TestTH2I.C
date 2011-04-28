void TestTH2I() {
	TCanvas *c2 = new TCanvas("c2","c2",600,400);
	TH2F *hsurf1 = new TH2F("hsurf1","Option SURF1 example ",30,-4,4,30,-20,20);
	Float_t px, py;
	for (Int_t i = 0; i < 25000; i++) {
		gRandom->Rannor(px,py);
		hsurf1->Fill(px-1,5*py);
		hsurf1->Fill(2+0.5*px,2*py-10.,0.1);
	}
	hsurf1->Draw("SURF1");

	TCanvas *c3 = new TCanvas("c3", "c3", 1000,800);

	TH2F *hnoise = new TH2F("hnoise", "noise histo", 30, -10., 10., 30, -10., 10.);

	Float_t GHz4[8] = {1.23,1.20,1.26,1.29,1.23,1.23,1.20,1.10};
	Float_t GHz2[8] = {0.74,0.67,0.67,0.67,0.67,0.92,0.67,0.52};
	Float_t GHz1[8] = {0.64,0.67,0.61,0.61,0.67,0.86,0.64,0.49};
	Float_t GHz05[8] = {0.64,0.67,0.64,0.64,0.67,0.86,0.67,0.52};

	for (Int_t j = 0; i < 8; j++) {
		hnoise->Fill(4,j, GHz4[j]);
	}
	for (Int_t j = 0; i < 8; j++) {
		hnoise->Fill(2,j,GHz2[j]);
	}

	for (Int_t j = 0; i < 8; j++) {
		hnoise->Fill(1,j,GHz1[j]);
	}

	for (Int_t j = 0; i < 8; j++) {
		hnoise->Fill(0.5,j,GHz05[j]);
	}

	hnoise->Draw("SURF1");

}
