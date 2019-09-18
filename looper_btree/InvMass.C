
float InvMass(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2){
    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
    p2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
    return (p1+p2).M();
}
