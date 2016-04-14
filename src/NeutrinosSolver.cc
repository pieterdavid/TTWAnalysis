#include <cp3_llbb/TTAnalysis/interface/NeutrinosSolver.h>

#include <Math/Vector3D.h>

std::vector<std::pair<NeutrinosSolver::LorentzVector, NeutrinosSolver::LorentzVector>> NeutrinosSolver::getNeutrinos(const LorentzVector& lepton1_p4, 
        const LorentzVector& lepton2_p4, 
        const LorentzVector& bjet1_p4, 
        const LorentzVector& bjet2_p4,
        const LorentzVector& met) {


    // pT = transverse total momentum of the visible particles
    // It will be used to reconstruct neutrinos, but we want to take into account the measured ISR (pt_isr = - pt_met - pt_vis),
    // so we add pt_isr to pt_vis in order to have pt_vis + pt_nu + pt_isr = 0 as it should be.

    const LorentzVector& p3 = lepton1_p4;
    const LorentzVector& p4 = bjet1_p4;
    const LorentzVector& p5 = lepton2_p4;
    const LorentzVector& p6 = bjet2_p4;

    double s13 = w_mass * w_mass;
    double s134 = t_mass * t_mass;
    double s25 = w_mass * w_mass;
    double s256 = t_mass * t_mass;

    LorentzVector ISR = -(lepton1_p4 + lepton2_p4 + bjet1_p4 + bjet2_p4 + met);
    LorentzVector pT = lepton1_p4 + lepton2_p4 + bjet1_p4 + bjet2_p4 + ISR;

    const double p34 = p3.Dot(p4);
    const double p56 = p5.Dot(p6);
    const double p33 = p3.M2();
    const double p44 = p4.M2();
    const double p55 = p5.M2();
    const double p66 = p6.M2();

    // A1 p1x + B1 p1y + C1 = 0, with C1(E1,E2)
    // A2 p1y + B2 p2y + C2 = 0, with C2(E1,E2)
    // ==> express p1x and p1y as functions of E1, E2

    const double A1 = 2.*( -p3.Px() + p3.Pz()*p4.Px()/p4.Pz() );
    const double A2 = 2.*( p5.Px() - p5.Pz()*p6.Px()/p6.Pz() );

    const double B1 = 2.*( -p3.Py() + p3.Pz()*p4.Py()/p4.Pz() );
    const double B2 = 2.*( p5.Py() - p5.Pz()*p6.Py()/p6.Pz() );

    const double Dx = B2*A1 - B1*A2;
    const double Dy = A2*B1 - A1*B2;

    const double X = 2*( pT.Px()*p5.Px() + pT.Py()*p5.Py() - p5.Pz()/p6.Pz()*( 0.5*(s25 - s256 + p66) + p56 + pT.Px()*p6.Px() + pT.Py()*p6.Py() ) ) + p55 - s25;
    const double Y = p3.Pz()/p4.Pz()*( s13 - s134 + 2*p34 + p44 ) - p33 + s13;

    // p1x = alpha1 E1 + beta1 E2 + gamma1
    // p1y = ...(2)
    // p1z = ...(3)
    // p2z = ...(4)
    // p2x = ...(5)
    // p2y = ...(6)

    const double alpha1 = -2*B2*(p3.E() - p4.E()*p3.Pz()/p4.Pz())/Dx;
    const double beta1 = 2*B1*(p5.E() - p6.E()*p5.Pz()/p6.Pz())/Dx;
    const double gamma1 = B1*X/Dx + B2*Y/Dx;

    const double alpha2 = -2*A2*(p3.E() - p4.E()*p3.Pz()/p4.Pz())/Dy;
    const double beta2 = 2*A1*(p5.E() - p6.E()*p5.Pz()/p6.Pz())/Dy;
    const double gamma2 = A1*X/Dy + A2*Y/Dy;

    const double alpha3 = (p4.E() - alpha1*p4.Px() - alpha2*p4.Py())/p4.Pz();
    const double beta3 = -(beta1*p4.Px() + beta2*p4.Py())/p4.Pz();
    const double gamma3 = ( 0.5*(s13 - s134 + p44) + p34 - gamma1*p4.Px() - gamma2*p4.Py() )/p4.Pz();

    const double alpha4 = (alpha1*p6.Px() + alpha2*p6.Py())/p6.Pz();
    const double beta4 = (p6.E() + beta1*p6.Px() + beta2*p6.Py())/p6.Pz();
    const double gamma4 = ( 0.5*(s25 - s256 + p66) + p56 + (gamma1 + pT.Px())*p6.Px() + (gamma2 + pT.Py())*p6.Py() )/p6.Pz();

    const double alpha5 = -alpha1;
    const double beta5 = -beta1;
    const double gamma5 = -pT.Px() - gamma1;

    const double alpha6 = -alpha2;
    const double beta6 = -beta2;
    const double gamma6 = -pT.Py() - gamma2;

    // a11 E1^2 + a22 E2^2 + a12 E1E2 + a10 E1 + a01 E2 + a00 = 0
    // id. with bij

    const double a11 = -1 + ( SQ(alpha1) + SQ(alpha2) + SQ(alpha3) );
    const double a22 = SQ(beta1) + SQ(beta2) + SQ(beta3);
    const double a12 = 2.*( alpha1*beta1 + alpha2*beta2 + alpha3*beta3 );
    const double a10 = 2.*( alpha1*gamma1 + alpha2*gamma2 + alpha3*gamma3 );
    const double a01 = 2.*( beta1*gamma1 + beta2*gamma2 + beta3*gamma3 );
    const double a00 = SQ(gamma1) + SQ(gamma2) + SQ(gamma3);

    const double b11 = SQ(alpha5) + SQ(alpha6) + SQ(alpha4);
    const double b22 = -1 + ( SQ(beta5) + SQ(beta6) + SQ(beta4) );
    const double b12 = 2.*( alpha5*beta5 + alpha6*beta6 + alpha4*beta4 );
    const double b10 = 2.*( alpha5*gamma5 + alpha6*gamma6 + alpha4*gamma4 );
    const double b01 = 2.*( beta5*gamma5 + beta6*gamma6 + beta4*gamma4 );
    const double b00 = SQ(gamma5) + SQ(gamma6) + SQ(gamma4);

    // Find the intersection of the 2 conics (at most 4 real solutions for (E1,E2))
    std::vector<double> E1, E2;
    solve2Quads(a11, a22, a12, a10, a01, a00, b11, b22, b12, b10, b01, b00, E1, E2);

    // For each solution (E1,E2), find the neutrino 4-momenta p1,p2
    std::vector<std::pair<LorentzVector, LorentzVector>> neutrinos;

    for (size_t i = 0; i < E1.size(); i++){
        const double e1 = E1.at(i);
        const double e2 = E2.at(i);

        if (e1 < 0. || e2 < 0.)
            continue;

        LorentzVector p1(
                alpha1*e1 + beta1*e2 + gamma1,
                alpha2*e1 + beta2*e2 + gamma2,
                alpha3*e1 + beta3*e2 + gamma3,
                e1);

        LorentzVector p2(
                alpha5*e1 + beta5*e2 + gamma5,
                alpha6*e1 + beta6*e2 + gamma6,
                alpha4*e1 + beta4*e2 + gamma4,
                e2);

        neutrinos.push_back(std::make_pair(p1, p2));
    }

    return neutrinos;
}

bool solveQuadratic(const double a, const double b, const double c, std::vector<double>& roots) {

    if(!a){
        if(!b){
            return false;
        }
        roots.push_back(-c/b);
        return true;
    }

    const double rho = SQ(b) - 4.*a*c;

    if(rho >= 0.){
        if(b == 0.){
            roots.push_back( std::sqrt(rho)/(2.*a) );
            roots.push_back( -std::sqrt(rho)/(2.*a) );
        }else{
            const double x = -0.5*(b + std::copysign(std::sqrt(rho), b));
            roots.push_back(x/a);
            roots.push_back(c/x);
        }
        return true;
    }else{
        return false;
    }
}

bool solveCubic(const double a, const double b, const double c, const double d, std::vector<double>& roots) {

    if(a == 0)
        return solveQuadratic(b, c, d, roots);

    const double an = b/a;
    const double bn = c/a;
    const double cn = d/a;

    const double Q = SQ(an)/9. - bn/3.;
    const double R = CB(an)/27. - an*bn/6. + cn/2.;

    if( SQ(R) < CB(Q) ){
        const double theta = std::acos( R/std::sqrt(CB(Q)) )/3.;

        roots.push_back( -2. * std::sqrt(Q) * std::cos(theta) - an/3. );
        roots.push_back( -2. * std::sqrt(Q) * cosXpm2PI3(theta, 1.) - an/3. );
        roots.push_back( -2. * std::sqrt(Q) * cosXpm2PI3(theta, -1.) - an/3. );
    }else{
        const double A = - std::copysign(std::cbrt(std::abs(R) + std::sqrt( SQ(R) - CB(Q))), R);

        double B;

        if(A == 0.)
            B = 0.;
        else
            B = Q/A;

        const double x = A + B - an/3.;

        roots.push_back(x);
        roots.push_back(x);
        roots.push_back(x);
    }

    return true;
}

bool solveQuartic(const double a, const double b, const double c, const double d, const double e, std::vector<double>& roots) {

    if(!a)
        return solveCubic(b, c, d, e, roots);

    if(!b && !c && !d){
        roots.push_back(0.);
        roots.push_back(0.);
        roots.push_back(0.);
        roots.push_back(0.);
    }else{
        const double an = b/a;
        const double bn = c/a - (3./8.) * SQ(b/a);
        const double cn = CB(0.5*b/a) - 0.5*b*c/SQ(a) + d/a;
        const double dn = -3.*QU(0.25*b/a) + e/a - 0.25*b*d/SQ(a) + c*SQ(b/4.)/CB(a);

        std::vector<double> res;
        solveCubic(1., 2.*bn, SQ(bn) - 4.*dn, -SQ(cn), res);
        short pChoice = -1;

        for(unsigned short i = 0; i<res.size(); ++i){
            if(res[i] > 0){
                pChoice = i;
                break;
            }
        }

        if(pChoice < 0){
            return false;
        }

        const double p = std::sqrt(res[pChoice]);
        solveQuadratic(p, SQ(p), 0.5*( p*(bn + res[pChoice]) - cn ), roots);
        solveQuadratic(p, -SQ(p), 0.5*( p*(bn + res[pChoice]) + cn ), roots);

        for(unsigned short i = 0; i<roots.size(); ++i)
            roots[i] -= an/4.;
    }

    size_t nRoots = roots.size();

    return nRoots > 0;
}

bool solve2Quads(const double a20, const double a02, const double a11, const double a10, const double a01, const double a00, const double b20, const double b02, const double b11, const double b10, const double b01, const double b00, std::vector<double>& E1, std::vector<double>& E2){

    // The procedure used in this function relies on a20 != 0 or b20 != 0
    if(a20 == 0. && b20 == 0.){

        if(a02 != 0. || b02 != 0.){
            // Swapping E1 <-> E2 should suffice!
            return solve2Quads(a02, a20, a11, a01, a10, a00,
                    b02, b20, b11, b01, b10, b00,
                    E2, E1);
        }else{
            return solve2QuadsDeg(a11, a10, a01, a00,
                    b11, b10, b01, b00,
                    E1, E2);
        }

    }

    const double alpha = b20*a02-a20*b02;
    const double beta = b20*a11-a20*b11;
    const double gamma = b20*a10-a20*b10;
    const double delta = b20*a01-a20*b01;
    const double omega = b20*a00-a20*b00;

    const double a = a20*SQ(alpha) + a02*SQ(beta) - a11*alpha*beta;
    const double b = 2.*a20*alpha*delta - a11*( alpha*gamma + delta*beta ) - a10*alpha*beta + 2.*a02*beta*gamma + a01*SQ(beta);
    const double c = a20*SQ(delta) + 2.*a20*alpha*omega - a11*( delta*gamma + omega*beta ) - a10*( alpha*gamma + delta*beta )
        + a02*SQ(gamma) + 2.*a01*beta*gamma + a00*SQ(beta);
    const double d = 2.*a20*delta*omega - a11*omega*gamma - a10*( delta*gamma + omega*beta ) + a01*SQ(gamma) + 2.*a00*beta*gamma;
    const double e = a20*SQ(omega) - a10*omega*gamma + a00*SQ(gamma);

    solveQuartic(a, b, c, d, e, E2);

    for(unsigned short i = 0; i < E2.size(); ++i){

        const double e2 = E2[i];

        if(beta*e2 + gamma != 0.){
            // Everything OK

            const double e1 = -(alpha * SQ(e2) + delta*e2 + omega)/(beta*e2 + gamma);
            E1.push_back(e1);

        }else if(alpha*SQ(e2) + delta*e2 + omega == 0.){
            // Up to two solutions for e1

            std::vector<double> e1;

            if( !solveQuadratic(a20, a11*e2 + a10, a02*SQ(e2) + a01*e2 + a00, e1) ){

                if( !solveQuadratic(b20, b11*e2 + b10, b02*SQ(e2) + b01*e2 + b00, e1) ){
                    E1.clear();
                    E2.clear();
                    return false;
                }

            }else{
                // We have either a double, or two roots for e1
                // In this case, e2 must be twice degenerate!
                // Since in E2 degenerate roots are grouped, E2[i+1] shoud exist and be equal to e2
                // We then go straight for i+2.

                if(i < E2.size() - 1){

                    if(e2 != E2[i+1]){
                        E1.clear();
                        E2.clear();
                        return false;
                    }

                    E1.push_back(e1[0]);
                    E1.push_back(e1[1]);
                    ++i;
                    continue;

                }else{
                    E1.clear();
                    E2.clear();
                    return false;
                }

            }

        }else{
            // There is no solution given this e2
            E2.erase(E2.begin() + i);
            --i;
        }

    }

    return true;
}

bool solve2QuadsDeg(const double a11, const double a10, const double a01, const double a00, const double b11, const double b10, const double b01, const double b00, std::vector<double>& E1, std::vector<double>& E2) {

    if(a11 == 0. && b11 == 0.)
        return solve2Linear(a10, a01, a00, b10, b01, b00, E1, E2);

    bool result = solveQuadratic(a11*(b11*a10-a11*b10),
            a01*(b11*a10-a11*b10) - a01*(b11*a01-a11*b01) + a11*(b11*a00-a11*b00),
            a01*(b11*a00-a11*b00),
            E1);

    if(!result){
        return false;
    }

    for(unsigned short i=0; i<E1.size(); ++i){

        double denom = a11*E1[i] + a01;

        if(denom != 0){
            E2.push_back( -(a10*E1[i] + a00)/denom );

        }else{
            denom = b11*a01 - a11*b01;
            if(denom != 0.){
                E2.push_back( -( (b11*a10 - a11*b10)*E1[i] + b11*a00 - a11*b00 )/denom );
            }else{
                denom = b11*E1[i] + b01;
                if(denom != 0.){
                    E2.push_back( -(b10*E1[i] + b00)/denom );
                }else{
                    E1.erase(E1.begin() + i);
                    --i;
                }
            }
        }

    }

    return E1.size();
}

bool solve2Linear(const double a10, const double a01, const double a00, const double b10, const double b01, const double b00, std::vector<double>& E1, std::vector<double>& E2) {

    const double det = a10*b01 - b10*a01;

    if(det == 0.){
        if(a00 != 0. || b00 != 0.){
            return false;
        }else{
            return false;
        }
    }

    const double e2 = (b10*a00-a10*b00)/det;
    E2.push_back(e2);
    double e1;
    if(a10 == 0.)
        e1 = -(b00 + b01*e2)/b10;
    else
        e1 = -(a00 + a01*e2)/a10;
    E1.push_back(e1);

    return true;
}
