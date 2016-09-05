#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>


//处理root 文件的底层函数
TH1F *ep_ana(TTree *t1,int opt1 ){

    const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};

    double cut1,cut2,cut3;

	// ---------------- settting -----------------

    TLorentzVector *km = new TLorentzVector();
    TLorentzVector *kp = new TLorentzVector();
    
    Double_t Chisq4c;

    double m_eppos_ratio;
    double m_epneg_ratio;

    double px,py,pz,p;
    double pxe,pye,pze,pe;

    double chiDeDx[4][5];
    double chiTof1[4][5];
    double chiTof2[4][5];
    double chisq_pid[4][5];
    
    double hits,miss;


    t1->SetBranchAddress("pneg", &pe);
    t1->SetBranchAddress("ppos", &p);
    //t1->SetBranchAddress("emchit_p", &hits);
    //t1->SetBranchAddress("emchit_p", &hits);
    t1->SetBranchAddress("EPRatioPlus",&m_eppos_ratio);
    t1->SetBranchAddress("EPRatioMinu",&m_epneg_ratio);
    t1->SetBranchAddress("Chisq_low",&Chisq4c);
    //t1->SetBranchAddress("chiDeDx",&chiDeDx);    
    //t1->SetBranchAddress("chiTof1",&chiTof1);
    //t1->SetBranchAddress("chiTof2",&chiTof2);
    //t1->SetBranchAddress("chisq_pid",&chisq_pid);
    t1->SetBranchAddress("km_unfitted",&km);
    t1->SetBranchAddress("kp_unfitted",&kp);
    

    TH1F *hpp = new TH1F("p-p","p-p",20,1.3,1.5);
    TH1F *hpe = new TH1F("p-e","p-e",20,1.3,1.5);

    //TH1F *hit = new TH1F("hit","EMC_hit",50,0.0,50);

//    TH1F *hmiss = new TH1F("u-miss","u-miss",20,0.0,0.2);
    TH1F *hmiss = new TH1F("u-miss","u-miss",70,0.0,0.07);

    TH1F *e1 = new TH1F("e/p-p","e/p-p",60,0,0.6);
    TH1F *e2 = new TH1F("e/p-e","e/p-e",60,0.6,1.2);

    TH1F *m1 = new TH1F("Chisq_4C","Chisq_4C",50,0.0,40);
/*
    TH1F *de1 = new TH1F("chi_dedx-p","chi_dedx-p",100,-5.0,10.0);
    TH1F *de2 = new TH1F("chi_dedx-e","chi_dedx-e",100,-5.0,10.0);

    TH1F *tof1p = new TH1F("chi_tof1-p","chi_tof1-p",100,-11.0,25.0);
    TH1F *tof1e = new TH1F("chi_tof1-e","chi_tof1-e",100,-11.0,25.0);

    TH1F *tof2p = new TH1F("chi_tof2-p","chi_tof2-p",100,-5.0,10.0);
    TH1F *tof2e = new TH1F("chi_tof2-e","chi_tof2-e",100,-5.0,5.0);

    TH1F *pid1 = new TH1F("chisq_pid-p","chisq_pid-p",100,0.0,200);
    TH1F *pid2 = new TH1F("chisq_pid-e","chisq_pid-e",100,0.0,40);

    TH1F *h2 = new TH1F("km","km",100,1.35,1.45);

    TH1F *h3 = new TH1F("kp","kp",100,1.65,1.75);
*/
   
    //-----------------------------------------------
	
	
	Long64_t nentries1 = t1->GetEntries();
    
    for (Long64_t i=0;i<nentries1;i++) {
        
        t1->GetEntry(i);
        
        TLorentzVector tot = *km + *kp;
        miss = 3.097 - (sqrt(pe*pe+xmass[0]*xmass[0])+sqrt(p*p+xmass[4]*xmass[4]) )+ tot.P();
        //if (Chisq4c<=40 && 1.3<=p && p<=1.5 && 1.3<= pe && pe<=1.5 ){
       /*
        if (Chisq4c <= 9 
            && 1.34 <= p 
            && p <= 1.46 
            && 1.34 <= pe 
            && pe <= 1.46 
            && miss <= 0.07 
            && 0.12 <= m_eppos_ratio 
            && m_eppos_ratio <= 0.3 
            && 0.9 <= m_epneg_ratio 
            && m_epneg_ratio <= 1.06){
        */    
            hpp->Fill(p);
            hpe->Fill(pe);

            //hit->Fill(hits);

            hmiss->Fill(miss);

            e1->Fill(m_eppos_ratio);
            e2->Fill(m_epneg_ratio);
            m1->Fill(Chisq4c);
/*
            tof1p->Fill(chiTof1[0][4]);
            tof1e->Fill(chiTof1[1][0]);

            tof2p->Fill(chiTof2[0][4]);
            tof2e->Fill(chiTof2[1][0]);

            de1->Fill(chiDeDx[0][4]);
            de2->Fill(chiDeDx[1][0]);

            pid1->Fill(chisq_pid[1][0]);
            pid2->Fill(chisq_pid[0][4]);
*/
        //}
               
    }
            //对于data 暂时没有the_decay,故自己写，否则要解除上面的注释h4
    if(opt1 ==1)return hpp;
    if(opt1 ==2)return hpe;
    //if(opt1 ==3)return hit;
    if(opt1 ==4)return hmiss;
    if(opt1 ==5)return e1;
    if(opt1 ==6)return e2;
    if(opt1 ==7)return m1;
    if(opt1 ==0)return nentries1;
    /*
    if(opt1 ==8)return tof1p;
    if(opt1 ==9)return tof1e;
    if(opt1 ==10)return tof2p;
    if(opt1 ==11)return tof2e;
    if(opt1 ==12)return de1;
    if(opt1 ==13)return de2;
    if(opt1 ==14)return pid1;
    if(opt1 ==15)return pid2;
    */
}