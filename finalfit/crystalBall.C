double CrystalBall(const double *x0, const double *p){
    double norm = p[0];
    double mean = p[1];
    double sigma = p[2];
    double alpha = p[3];
    double n = p[4];
    double x = x0[0];

    double t = (x-mean)/sigma;
    if(alpha < 0)
        t = -t;

    double absAlpha = fabs(alpha);
    
    if(t >= -absAlpha)
        return norm * exp(-0.5*t*t);
    else{
        double a =  TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
        double b= n/absAlpha - absAlpha;

        return norm * a/TMath::Power(b - t, n);
    }
}

bool REJECT = false;
double Exp(const double *x0, const double *p){
    double x = x0[0];
    if(REJECT && x>120 && x<130){
        TF1::RejectPoint();
        return 0;
    }
    return p[0]/(p[1]*(1-exp(-80/p[1])))*exp(-(x-100) / p[1]);
}
double Quadratic(const double *x0, const double *p){
    double x = x0[0];
    if(REJECT && x>120 && x<130){
        TF1::RejectPoint();
        return 0;
    }
    return p[0] + p[1]*x + p[2]*x*x;
}
double Power(const double *x0, const double *p){
    double x = x0[0];
    if(REJECT && x>120 && x<130){
        TF1::RejectPoint();
        return 0;
    }
    return p[0] * TMath::Power(x, p[1]) + p[2] * TMath::Power(x, p[3]);
}
