#include "ep_ana.cxx"

void ana0518(){
	
int opt=1;
/*
    if(opt1 ==1)return hpp;
    if(opt1 ==2)return hpe;
    if(opt1 ==4)return hmiss;
    if(opt1 ==5)return e1;
    if(opt1 ==6)return e2;
    if(opt1 ==7)return m1;
*/

// inclusive MC
    int n1 = 1;int n2 = 91;
	
    TH1F *hh = read(n1,n2,opt,0);

//con-data

	int n1=1,n2=18;

    TH1F *chh = read(n1,n2,opt,1);



//bhabha_mc

    int n1=62,n2=104;

    TH1F *bhh = read(n1,n2,opt,2);


//data,jpsi
    int n1=1,n2=350;

    TH1F *dhh = read(n1,n2,opt,3);


//signal

	TFile *f6 = new TFile("../root/ep0512.root");
	TTree *t6 = (TTree*)f6->Get("TreeAna");
	TH1F *h6 = ep_ana(t6,opt);

// s/sqrt(s+b) cut
/*
    int N = h6->GetNbinsX();
    double sig,bkg,ry,rx;
    double ratio = 0;
    TH1F *hr = ep_ana(t6,opt);
    hr->SetTitle("ratio_s_sqrt(b)");

    for(int j = 0;j<N;j++){
        sig = sig + h6->GetBinContent(j);
        bkg = bkg + hh->GetBinContent(j);
        if(bkg>=1){
            ry = sig/sqrt(bkg);
            if( ry>=ratio){ratio = ry; rx = 0.01*j;}
            hr->SetBinContent(j,ry);
            cout<<j<<" "<<ry<<endl;
        }
        else{hr->SetBinContent(j,0);}

    }
    cout<<"max "<<rx<<" "<<ratio<<endl;
*/    

    h6->SetLineColor(kRed);h6->SetStats(0);
    hh->SetLineColor(kGray);hh->SetStats(0); //h1->SetTitle();
	chh->SetLineColor(kBlue);chh->SetStats(0);
	bhh->SetLineColor(kGreen);bhh->SetStats(0);
    dhh->SetLineColor(kBlack);dhh->SetStats(0);
    dhh->SetMarkerStyle(20);

        
    



    //cocktail MC
    double norm0=1.0;
    double norm3=dhh->Integral();//1.0;

    double norm1=0.74;
    double norm2=0.30;
    
    Double_t scale0 = norm0/hh->Integral();
    hh->Scale(scale0);    
    Double_t scale3 = norm3/dhh->Integral();
    dhh->Scale(scale3);

    Double_t scale1 = norm1/chh->Integral();
    chh->Scale(norm3*scale1);
    Double_t scale2 = norm2/bhh->Integral();
    bhh->Scale(norm3*scale2);

    //ADD:chh,bhh

    TH1F *cbhh = chh->Clone();
    cbhh->Add(bhh);cbhh->SetLineColor(kCyan);
    
    TLegend *l = new TLegend(0.15,0.6,0.3,0.75);
    l->AddEntry(hh,"inclusive_mc");
    l->AddEntry(chh,"continuum");
    //l->AddEntry(h6,"signal");
    l->AddEntry(bhh,"jpsitoee");
    l->AddEntry(dhh,"data");
    l->AddEntry(cbhh,"total_c+b");

    //box for blind
    int Y =6500;
    TBox *b = new TBox(1.37,0,1.44,Y);
    b->SetFillColor(kRed-5);

	//TCanvas *c1 = new TCanvas("c1","",700,500);c1->Divide(1,1);c1->cd(1);//c1->SetTitle("p-neg");
    dhh->GetYaxis()->SetRangeUser(0,Y);
    dhh->Draw("PE");
    
    //h6->DrawNormalized("same");
    chh->Draw("same");
	//hh->Draw("same");
    bhh->Draw("same");
    cbhh->Draw("same");	
    l->Draw("same");
    
    b->Draw("same");
    
/*    
    TLine *t1 = new TLine(rx,0,rx,h6->GetBinContent(h6->GetMaximumBin()));
    t1->SetLineColor(kOrange);t1->SetLineStyle(2);t1->Draw("same");
*/    
   /*
	c1->cd(2);
    hr->Draw();
	TLine *t2 = new TLine(rx,0,rx,2*hr->GetBinContent(hr->GetMaximumBin()));
    t2->SetLineColor(kOrange);t2->SetLineStyle(2);t2->Draw("same");
    */
   
}

TH1F *read(int n1, int n2,  int opt,int opt2){
    int n = n2-n1+1;

    ostringstream s3[n];

    for(int i=n1; i <= n2;i++){
        
        int p = i-n1;        
        //s[i] << "../p3770Alg/dataroot/0414" << p << ".root"<<endl;
        
        if(opt2==3){s3[p]<<"../root/data/jpsidata_"<<i<<".root";}
        if(opt2==2){s3[p]<<"../root/bhabha_mc/bhabha_"<<i<<".root";}
        if(opt2==1){s3[p]<<"../root/con_mc/con_mc_"<<i<<".root";}
        if(opt2==0){s3[p]<<"../root/mc/mc_"<<i<<".root";}
//        cout << s3[p].str() << endl;

        const char* r = s3[p].str().c_str();//"test.root";
          
        cout<<r<<endl;
    
        TFile *f = new TFile(r);

        if (f->IsOpen()) {
            
            if ( (TTree*)f->Get("TreeAna") ){

                if(i==n1)TH1F *dhh = ep_ana((TTree*)f->Get("TreeAna"),opt); 
                
                TH1F *h = ep_ana((TTree*)f->Get("TreeAna"),opt); 
                //h->Draw();
                if(i!=n1)dhh->Add(h);       
            }
        }

        f->Close();

    }

    return dhh;

}



