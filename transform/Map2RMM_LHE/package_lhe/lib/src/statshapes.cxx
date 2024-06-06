// user calculations inside loop over MC events 
using namespace std;
#include<iostream>
#include "TMath.h"
#include "LParticle.h"
#include "CParticle.h"
#include "TLorentzVector.h"

// see implementation details in
// http://ask.metafilter.com/36213/Best-Fit-Ellipse


double* statshapes(LParticle jet) {


	int i, j;
	vector<CParticle> constit = jet.GetConstituents();
	double d_eta2; double d_etaphi;
	double wm_eta, wm_phi;
	double a, b, a0, b0;
	double d_eta;
	double d_phi;
	double eTot, Fmax;
        double *  stat = new double[23]; 


        // remove bad constituenst (with pt=0)
        unsigned int NTT=0;
        for (i = 0; i < constit.size(); i++){
            CParticle cp = constit.at(i);
            if (cp.calcPt()<=0.0) NTT++; 
        }
 
        int num =  constit.size()-NTT;
        double eta[num], phi[num], e[num], eta1[num], phi1[num];


	//get eta, phi, & e from CParticle/TLorentzVector and store in seperate vectors
	for (i = 0; i < num; i++){
		CParticle cp = constit.at(i);
                if (cp.calcPt()<=0.0) continue;
		  eta1[i] = cp.calcEta();
		  phi1[i] = cp.calcPhi();
		  e[i] = cp.e();
	}

		
	//invert jet if vertical for linear regression
	double mineta = 0; double maxeta = 0; double minphi = 0; double maxphi = 0;
	for (i = 0; i < num; i++){ 
		mineta = min(mineta, eta1[i]);
		maxeta = max(maxeta, eta1[i]);
		minphi = min(minphi, phi1[i]);
		maxphi = max(maxphi, phi1[i]);
	}
	if (maxphi - minphi > maxeta - mineta){
		for (i = 0; i < num; i ++){
			eta[i] = phi1[i];
			phi[i] = eta1[i];
		}
	}
	else {
		for (i = 0; i < num; i ++){
			eta[i] = eta1[i];
			phi[i] = phi1[i];
		}
	}
		
	
	//perform UNweighted linear regression
	eTot = 0;
	wm_eta = 0;
	wm_phi = 0;
	d_eta = 0; d_phi = 0; d_eta2 = 0; d_etaphi = 0;
	for (j = 0; j < num; j++) {
		double ei = 1;
		if (e[j] > 0) {
			eTot += ei;
			wm_eta += (eta[j] - wm_eta) * (ei / eTot);
			wm_phi += (phi[j] - wm_phi) * (ei / eTot);
		}
	}
	for (j = 0; j < num; j++) {
		if (e[j] > 0) {
			d_eta = eta[j] - wm_eta;
			d_phi = phi[j] - wm_phi;
			d_eta2 += (d_eta * d_eta) / eTot;
			d_etaphi += (d_eta * d_phi) / eTot;
		}
	}
	
	// recalculate eTot in case of unweighted linear regression 
	eTot = 0;
	for (i = 0; i < num; i++) {
		eTot += e[i];
	} 

	//linear regression (global avg. for major axis)
	//i.e. phi = a + eta*b
	b = d_etaphi / d_eta2;
	a = wm_phi - (wm_eta * b);

	//perpendicular line through mean
	b0 = -1 / b;
	a0 = wm_phi - (b0 * wm_eta);

	//rotate axes by 45 deg. to form quadrants
	double rotang = atan(b) + (TMath::Pi() / 4);
	double rotangp = atan(b0) + (TMath::Pi() / 4);
	double newb = TMath::Tan(rotang);
	double newb0 = TMath::Tan(rotangp);
	double newa = wm_phi - (newb * wm_eta);
	double newa0 = wm_phi - (newb0 * wm_eta);

	//find weighted center of each quadrant
	double dist[2][num];
        double points[2][num];
        double NewPoints[2][num];
        double R[2][2];
        for (i = 0; i < num; i++){
                points[0][i] = eta[i];
                points[1][i] = phi[i];
        }
        double theta = atan(b);
        R[0][0] = TMath::Cos(theta);
        R[0][1] = TMath::Sin(theta);
        R[1][0] = -TMath::Sin(theta);
        R[1][1] = TMath::Cos(theta);
        for (i = 0; i < 2; i++){
                for (j = 0; j < num; j++){
                        NewPoints[i][j] = (R[i][0] * points[0][j]) + (R[i][1] * points[1][j]);
                        dist[0][j] = NewPoints[0][j] - wm_eta;
                        dist[1][j] = NewPoints[1][j] - wm_phi;
                }
        }
	double wwm[4][2];        //coordiantes of weighted quadrant centers
	double Ep[4];            //partial energy for each region (quadrant or nq)
	double wwm1[4][2];       //coordinates of quadrant centers in adjusted basis
	double a_wwm[4][2];      //coordinates of nq centers
	double a1_wwm[4][2];     //coordinates of nq centers in adjusted basis
	for (i = 0; i < 4; i++){
		Ep[i] = 0;
		for (j = 0; j < 2; j++){
			wwm[i][j] = 0;
			wwm1[i][j] = 0;
			a_wwm[i][j] = 0;
			a1_wwm[i][j] = 0;
		}
	}
	for (i = 0; i < num; i++){
		double ei = e[i];
		if (ei > 0){
			//[0] above major axis
			if (phi[i] > a + b * eta[i]){
				Ep[0] += ei;
				a_wwm[0][0] += (eta[i] - a_wwm[0][0]) * (ei / Ep[0]);
				a_wwm[0][1] += (phi[i] - a_wwm[0][1]) * (ei / Ep[0]);
				a1_wwm[0][0] += (dist[0][i] - a1_wwm[0][0]) * (ei / Ep[0]);
				a1_wwm[0][1] += (dist[1][i] - a1_wwm[0][1]) * (ei / Ep[0]);
			}
			//[1] below major axis
			if (phi[i] < a + b * eta[i]){
				Ep[1] += ei;
				a_wwm[1][0] += (eta[i] - a_wwm[1][0]) * (ei / Ep[1]);
				a_wwm[1][1] += (phi[i] - a_wwm[1][1]) * (ei / Ep[1]);
				a1_wwm[1][0] += (dist[0][i] - a1_wwm[1][0]) * (ei / Ep[1]);
				a1_wwm[1][1] += (dist[1][i] - a1_wwm[1][1]) * (ei / Ep[1]);
			}
			//[2] above minor axis
			if (phi[i] > a0 + b0 * eta[i]){
				Ep[2] += ei;
				a_wwm[2][0] += (eta[i] - a_wwm[2][0]) * (ei / Ep[2]);
				a_wwm[2][1] += (phi[i] - a_wwm[2][1]) * (ei / Ep[2]);
				a1_wwm[2][0] += (dist[0][i] - a1_wwm[2][0]) * (ei / Ep[2]);
				a1_wwm[2][1] += (dist[1][i] - a1_wwm[2][1]) * (ei / Ep[2]);
			}
			//[3] below minor axis
			if (phi[i] < a0 + b0 * eta[i]){
				Ep[3] += ei;
				a_wwm[3][0] += (eta[i] - a_wwm[3][0]) * (ei / Ep[3]);
				a_wwm[3][1] += (phi[i] - a_wwm[3][1]) * (ei / Ep[3]);
				a1_wwm[3][0] += (dist[0][i] - a1_wwm[3][0]) * (ei / Ep[3]);
				a1_wwm[3][1] += (dist[1][i] - a1_wwm[3][1]) * (ei / Ep[3]);
			}
		}
	}
	for (i = 0; i < 4; i++){
                Ep[i] = 0;
	}
		for (i = 0; i < num; i++){
        	        double ei = e[i];
                	if (ei > 0){
                        	//[0] in top quadrant
                        	if (phi[i] > newa + (newb * eta[i]) && phi[i] > newa0 + (newb0 * eta[i])){
                                	Ep[0] += ei;
                               		wwm[0][0] += (eta[i] - wwm[0][0]) * (ei / Ep[0]);
                               		wwm[0][1] += (phi[i] - wwm[0][1]) * (ei / Ep[0]);
                                	wwm1[0][0] += (dist[0][i] - wwm1[0][0]) * (ei / Ep[0]);
                                	wwm1[0][1] += (dist[1][i] - wwm1[0][1]) * (ei / Ep[0]);
                        	}
                        	//[1] in left quadrant
                        	if (phi[i] > newa + (newb * eta[i]) && phi[i] < newa0 + (newb0 * eta[i])){
                                	Ep[1] += ei;
                                	wwm[1][0] += (eta[i] - wwm[1][0]) * (ei / Ep[1]);
                                	wwm[1][1] += (phi[i] - wwm[1][1]) * (ei / Ep[1]);
                                	wwm1[1][0] += (dist[0][i] - wwm1[1][0]) * (ei / Ep[1]);
                              		wwm1[1][1] += (dist[1][i] - wwm1[1][1]) * (ei / Ep[1]);
                        	}
                        	//[2] in bottom quadrant
                        	if (phi[i] < newa + (newb * eta[i]) && phi[i] < newa0 + (newb0 * eta[i])){
                                	Ep[2] += ei;
                                	wwm[2][0] += (eta[i] - wwm[2][0]) * (ei / Ep[2]);
                                	wwm[2][1] += (phi[i] - wwm[2][1]) * (ei / Ep[2]);
                                	wwm1[2][0] += (dist[0][i] - wwm1[2][0]) * (ei / Ep[2]);
                                	wwm1[2][1] += (dist[1][i] - wwm1[2][1]) * (ei / Ep[2]);
                        	}
                        	//[3] in right quadrant
                        	if (phi[i] < newa + (newb * eta[i]) && phi[i] > newa0 + (newb0 * eta[i])){
                                	Ep[3] += ei;
                                	wwm[3][0] += (eta[i] - wwm[3][0]) * (ei / Ep[3]);
                                	wwm[3][1] += (phi[i] - wwm[3][1]) * (ei / Ep[3]);
                                	wwm1[3][0] += (dist[0][i] - wwm1[3][0]) * (ei / Ep[3]);
                                	wwm1[3][1] += (dist[1][i] - wwm1[3][1]) * (ei / Ep[3]);
                        	}
                	}
        	}

	//if quadrant empty, place quadrant-center at jet-center
	for (i = 0; i < 4; i++){
		if (Ep[i] == 0){
        		wwm[i][0] = wm_eta;
                	wwm[i][1] = wm_phi;
                	wwm1[i][0] = wm_eta;
                	wwm1[i][1] = wm_phi;
                }
       }

	//get Fmax parameter
	double emax = 0;
	for (i = 0; i < num; i++) {
		if (e[i] > emax) {
			emax = e[i];
		}
	}
	Fmax = emax / eTot;


	//distances between opposing quadrant centers (length of axes)
	double a1 = TMath::Sqrt((wwm[0][0] - wwm[2][0]) * (wwm[0][0] - wwm[2][0]) 
			+ (wwm[0][1] - wwm[2][1]) * (wwm[0][1] - wwm[2][1]));
	double a2 = TMath::Sqrt((wwm[3][0] - wwm[1][0]) * (wwm[3][0] - wwm[1][0]) 
			+ (wwm[3][1] - wwm[1][1]) * (wwm[3][1] - wwm[1][1]));
	//length of axes method 2
	double a7 = TMath::Sqrt((wwm1[0][1] - wwm1[2][1]) * (wwm1[0][1] - wwm1[2][1]));
        double a8 = TMath::Sqrt((wwm1[3][0] - wwm1[1][0]) * (wwm1[3][0] - wwm1[1][0]));

	//distances between jet center and non-quadrant centers
	//top semiaxis
	double a3 = TMath::Sqrt((wm_eta - a_wwm[0][0]) * (wm_eta - a_wwm[0][0]) 
			+ (wm_phi - a_wwm[0][1]) * (wm_phi - a_wwm[0][1]));
	//bottom semiaxis
	double a4 = TMath::Sqrt((wm_eta - a_wwm[1][0]) * (wm_eta - a_wwm[1][0]) 
			+ (wm_phi - a_wwm[1][1]) * (wm_phi - a_wwm[1][1]));
	//left semiaxis
	double a5 = TMath::Sqrt((wm_eta - a_wwm[2][0]) * (wm_eta - a_wwm[2][0]) 
			+ (wm_phi - a_wwm[2][1]) * (wm_phi - a_wwm[2][1]));
	//right semiaxis
	double a6 = TMath::Sqrt((wm_eta - a_wwm[3][0]) * (wm_eta - a_wwm[3][0]) 
			+ (wm_phi - a_wwm[3][1]) * (wm_phi - a_wwm[3][1]));

	//length of non-quandrant axes
	double a9 = TMath::Sqrt((a_wwm[0][0] - a_wwm[1][0]) * (a_wwm[0][0] - a_wwm[1][0])
                        + (a_wwm[0][1] - a_wwm[1][1]) * (a_wwm[0][1] - a_wwm[1][1]));
	double a10 = TMath::Sqrt((a_wwm[2][0] - a_wwm[3][0]) * (a_wwm[2][0] - a_wwm[3][0])
                        + (a_wwm[2][1] - a_wwm[3][1]) * (a_wwm[2][1] - a_wwm[3][1]));

	//distance between jet center and quadrant centers
	//top semiaxis
	double a11 = TMath::Sqrt((wm_eta - wwm[0][0]) * (wm_eta - wwm[0][0])
			+ (wm_phi - wwm[0][1]) * (wm_phi - wwm[0][1]));
	//left semiaxis
	double a12 = TMath::Sqrt((wm_eta - wwm[1][0]) * (wm_eta - wwm[1][0])
			+ (wm_phi - wwm[1][1]) * (wm_phi - wwm[1][1]));
	//bottom semiaxis
	double a13 = TMath::Sqrt((wm_eta - wwm[2][0]) * (wm_eta - wwm[2][0])
			+ (wm_phi - wwm[2][1]) * (wm_phi - wwm[2][1]));
	//right semiaxis
	double a14 = TMath::Sqrt((wm_eta - wwm[3][0]) * (wm_eta - wwm[3][0])
			+ (wm_phi - wwm[3][1]) * (wm_phi - wwm[3][1]));

	//length of non-quadrant axis method 2
	double a15 = TMath::Sqrt((a1_wwm[0][1] - a1_wwm[1][1]) * (a1_wwm[0][1] - a1_wwm[1][1]));
	double a16 = TMath::Sqrt((a1_wwm[2][0] - a1_wwm[3][0]) * (a1_wwm[2][0] - a1_wwm[3][0]));

	//get eccentricity
	double minorLength = a1;
	double majorLength = a2;
	if (a1 > a2){
		minorLength = a2;
		majorLength = a1;
	}


	//get eccentricity method 2
	double minorLength_meth2 = a7;
	double majorLength_meth2 = a8;
	if (a7 > a8){
		minorLength_meth2 = a8;
		majorLength_meth2 = a7;
	}


	//get major eccentricity
	double majorLength1 = a12;
	double majorLength2 = a14;
	if (a12 > a14){
		majorLength1 = a14;
		majorLength2 = a12;
	}

	//get minor eccentricity
	double minorLength1 = a11;
	double minorLength2 = a13;
	if (a11 > a13){
		minorLength1 = a13;
		minorLength2 = a11;
	}

	//get non-quadrant major eccentricity
	double nq_majorLength1 = a5;
	double nq_majorLength2 = a6;
	if (a5 > a6){
		nq_majorLength1 = a6;
		nq_majorLength2 = a5;
	}


	//get non-quadrant minor eccentricity
	double nq_minorLength1 = a3;
	double nq_minorLength2 = a4;
	if (a3 > a4){
		nq_minorLength1 = a4;
		nq_minorLength2 = a3;
	}


	//get non-quadrant eccentricity
	double nq_minorLength = a9;
	double nq_majorLength = a10;
	if (a9 > a10){
		nq_minorLength = a10;
		nq_majorLength = a9;
	}

	//get non-quadrant eccentricity method 2
	double nq_minorLength_meth2 = a15;
	double nq_majorLength_meth2 = a16;
	if (a15 > a16){
		nq_minorLength_meth2 = a16;
		nq_majorLength_meth2 = a15;
	}


       //get absolute jet size.
       //find dist from jet center in new basis, and find max-dist constituents
	
	double min_eta = dist[0][0];
	double max_eta = dist[0][0];
	double min_phi = dist[1][0];
	double max_phi = dist[1][0];
	int l = 0;
	int q = 0;
	int r = 0;
	int s = 0;
	//find eta' and phi' points furthest from jet center
	for (i = 1; i < num; i++){
		if (dist[0][i] > max_eta){
			max_eta = dist[0][i];
			l = i;
		}
		if (dist[0][i] < min_eta){
			min_eta = dist[0][i];
			q = i;
		}
		if (dist[1][i] > max_phi){
			max_phi = dist[1][i];
			r = i;
		}
		if (dist[1][i] < min_phi){
			min_phi = dist[1][i];
			s = i;
		}
	}
	double GlobalMajor = dist[0][l] - dist[0][q];    //global major length
	double GlobalMinor = dist[1][r] - dist[1][s];    //global minor length

	stat[0] = majorLength;
	stat[1] = minorLength;
	stat[2] = 1 - minorLength / majorLength;
	stat[3] = majorLength1;
	stat[4] = majorLength2;
	stat[5] = 1 - majorLength1 / majorLength2;
	stat[6] = minorLength1;
	stat[7] = minorLength2;
	stat[8] = 1 - minorLength1 / minorLength2;
	stat[9] = GlobalMajor;
	stat[10] = GlobalMinor;
	stat[11] = majorLength_meth2;
	stat[12] = minorLength_meth2;
	stat[13] = 1 - minorLength_meth2 / majorLength_meth2;
	stat[14] = nq_majorLength;
	stat[15] = nq_minorLength; 
	stat[16] = 1 - nq_minorLength / nq_majorLength;
	stat[17] = nq_majorLength_meth2;
	stat[18] = nq_minorLength_meth2;
	stat[19] = 1 - nq_minorLength_meth2 / nq_majorLength_meth2;
	stat[20] = 1 - nq_majorLength1 / nq_majorLength2;
	stat[21] = 1 - nq_minorLength1 / nq_minorLength2;
	stat[22] = Fmax;
	//if minorECC = NaN, replace with 0
	if (stat[8] != stat[8]) {stat[8] = 0;}	
    
	return stat;

} 

