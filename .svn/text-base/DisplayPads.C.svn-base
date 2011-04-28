void DisplayPads(char *infile="output.dat") {

	#define DOMINO_NCH 8
	#define DOMINO_NCELL 1024

	struct channel_struct {
		unsigned short data[DOMINO_NCELL];
		} channel_data;

		struct event_struct {
			unsigned short tcell;
			unsigned short depth;
			unsigned short dummy1;
			unsigned short dummy2;
			unsigned short dummy3;
			unsigned short dummy4;
			unsigned short dummy5;
			unsigned short dummy6;
			unsigned short dummy7;
			unsigned short dummy8;
			unsigned short dummy9;
			unsigned short dummy10;
			unsigned short dummy11;
			unsigned short dummy12;
			unsigned short dummy13;
			unsigned short dummy14;
			struct channel_struct ch[DOMINO_NCH];
			} event_data;

			FILE *fdata;
			int i,j, n;
			struct channel_struct *p;

			TString Title, TitID, Titoli[8]={"ch00","ch01","ch02","ch03","ch04","ch05","ch06","ch07"};

			cout << " Input file: " << infile << endl;

			fdata=fopen(infile,"r");

			TCanvas *c1 = new TCanvas("c1",Title,1200,700);
			c1->Divide(1,8,0.001,0.001);
			ix=1;

			Title="Channel 1";
			TitID="ch "; 
			TH1F *ch0 = new TH1F(Titoli[0],Titoli[0],1024,0,1023);
			TH1F *ch1 = new TH1F(Titoli[1],Titoli[1],1024,0,1023);
			TH1F *ch2 = new TH1F(Titoli[2],Titoli[2],1024,0,1023);
			TH1F *ch3 = new TH1F(Titoli[3],Titoli[3],1024,0,1023);
			TH1F *ch4 = new TH1F(Titoli[4],Titoli[4],1024,0,1023);
			TH1F *ch5 = new TH1F(Titoli[5],Titoli[5],1024,0,1023);
			TH1F *ch6 = new TH1F(Titoli[6],Titoli[6],1024,0,1023);
			TH1F *ch7 = new TH1F(Titoli[7],Titoli[7],1024,0,1023);
			ch0->SetMinimum(1900);
			ch1->SetMinimum(1900);
			ch2->SetMinimum(1900);
			ch3->SetMinimum(1900);
			ch4->SetMinimum(1900);
			ch5->SetMinimum(1900);
			ch6->SetMinimum(1900);
			ch7->SetMinimum(1900);

			
// questo ciclo while riempe ogni canale con la media dei piedistalli per tutti gli eventi
				
			n=fread((void *)&event_data,1,sizeof(event_data),fdata);			

			while(n>0) {	
			
				p=(struct channel_struct *)&event_data.ch[0];
				
				for(i=0;i<1024;i++)
				{	
					ch0->Fill(i,p->data[i]);
				}
				p++;
				for(i=0;i<1024;i++)
				{	
					ch1->Fill(i,p->data[i]);
				}
				p++;
				for(i=0;i<1024;i++)
				{	
					ch2->Fill(i,p->data[i]);
				}
				p++;
				for(i=0;i<1024;i++)
				{	
					ch3->Fill(i,p->data[i]);
				}
				p++;
				for(i=0;i<1024;i++)
				{	
					ch4->Fill(i,p->data[i]);
				}
				p++;
				for(i=0;i<1024;i++)
				{	
					ch5->Fill(i,p->data[i]);
				}
				p++;
				for(i=0;i<1024;i++)
				{	
					ch6->Fill(i,p->data[i]);
				}
				p++;
				for(i=0;i<1024;i++)
				{	
					ch7->Fill(i,p->data[i]);
				}
				n=fread((void *)&event_data,1,sizeof(event_data),fdata);			
			}



			c1->cd( ix++); ch0->DrawCopy();
			c1->cd( ix++); ch1->Draw();
			c1->cd( ix++); ch2->Draw();
			c1->cd( ix++); ch3->Draw();
			c1->cd( ix++); ch4->Draw();
			c1->cd( ix++); ch5->Draw();
			c1->cd( ix++); ch6->Draw();
			c1->cd( ix++); ch7->Draw();
		
			TCanvas *c2 = new TCanvas("c2",Title,1200,700);
			c2->cd(1); ch0->Draw();

			fclose(fdata);
		}
	
