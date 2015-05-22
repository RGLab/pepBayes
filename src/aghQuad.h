#ifndef AGHQUAD_H_
#define AGHQUAD_H_

class densityAlphaGamma0 {
public:
    densityAlphaGamma0(const double & nu, const double & n_rep) :
        log_gamma_ratio_data(R::lgammafn((n_rep + nu) / 2.0) -
                R::lgammafn(nu / 2.0)),
        log_gamma_ratio_ranef(R::lgammafn((1 + nu) / 2.0) -
                R::lgammafn(nu / 2.0)) {} ;

    void calc(arma::vec & out, const arma::vec & alpha,
            const double & offset0, const double & delta0,
            const double & offset1, const double & delta1,
            const int & nu, const double & nu_eps, const double & nu_alpha,
            const double & n_rep);
private:
    const double log_gamma_ratio_data;
    const double log_gamma_ratio_ranef;
};

class densityAlphaGamma1 {
public:
    densityAlphaGamma1(const double & nu, const double & n_rep) :
        log_gamma_ratio_data(R::lgammafn((n_rep + nu) / 2.0) -
                R::lgammafn(nu / 2.0)),
        log_gamma_ratio_ranef(R::lgammafn((1 + nu) / 2.0) -
                R::lgammafn(nu / 2.0)) {} ;

    void calc(arma::vec & out, const arma::vec & alpha,
            const double & offset, const double & delta,
            const int & nu, const double & nu_eps, const double & nu_alpha,
            const double & n_rep);
private:
    const double log_gamma_ratio_data;
    const double log_gamma_ratio_ranef;
};

double aghQuad(const arma::vec & gauss_rule_w,
        const arma::vec & gauss_rule_x,
        const arma::vec & gauss_work_gx,
        const double & sigma);

void densityGamma1(double * out, const double & offset0, const double & delta0,
        const double & offset1, const double & delta1, const int & nu,
        const double & nu_eps, const double & nu_alpha0, const double & nu_alpha1,
        const double & n_rep,
        const double & m0, const double & m1,
        const double & sigma0, const double & sigma1,
        const arma::vec & gauss_rule_x, const arma::vec & gauss_rule_w,
        arma::vec & gauss_work_x, arma::vec & gauss_work_gx);

double densityGamma0(const double & offset0, const double & delta0,
        const double & offset1, const double & delta1, const int & nu,
        const double & nu_eps, const double & nu_alpha, const double & n_rep,
        const arma::vec & gauss_rule_x, const arma::vec & gauss_rule_w,
        arma::vec & gauss_work_x, arma::vec & gauss_work_gx);

#endif
